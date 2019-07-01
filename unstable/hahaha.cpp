//
// hahaha
//
// Process HAHA data to find genome phase
//
// Copyright 2018-2019 Peter Andrews @ CSHL
//

// ~/mumdex/hahaha /data/safe/paa/analysis/mums/hg38/hg38.fa 500 160 chr1

#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>

#include <algorithm>
#include <exception>
#include <iostream>
#include <fstream>
#include <functional>
#include <iomanip>
#include <memory>
#include <random>
#include <sstream>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include "cn.h"
#include "error.h"
#include "logs.h"
#include "matrix.h"
#include "mumdex.h"
#include "psplot.h"
#include "stats.h"
#include "threads.h"

using std::cout;
using std::cerr;
using std::endl;
using std::exception;
using std::fixed;
using std::future;
using std::lock_guard;
using std::make_unique;
using std::max;
using std::min;
using std::mutex;
using std::ofstream;
using std::ostringstream;
using std::pair;
using std::string;
using std::to_string;
using std::unique_lock;
using std::unique_ptr;
using std::vector;

// Main algorithm control points
const bool use_old_coverage{1};
const uint64_t max_em_iter{200};
const unsigned int n_em_trials{3};
constexpr bool cut_hets{false};
const size_t afterburner_range{15};
const uint64_t sample_coverage_low_cutoff{1};
const uint64_t sample_coverage_high_cutoff{1000000000};  // basically unused
const uint64_t het_coverage_low_cutoff{1};
const uint64_t het_coverage_high_cutoff{1000000000};  // basically unused
const double max_alleles{50};

namespace paa {

template <class Val, uint64_t I>
using Array = std::array<Val, I>;

using uArray4 = Array<unsigned int, 4>;

constexpr uArray4 flip_haplo{{0, 2, 1, 3}};
constexpr uArray4 flip_obs{{0, 2, 1, 3}};

class HetInfo {
 public:
  HetInfo(const unsigned int chromosome__,
          const unsigned int position__,
          const int known_phase_info) :
      position_{position__},
    chromosome_{static_cast<unsigned char>(chromosome__)},
      is_phased_{known_phase_info >= 0},
      mom_is_a_{known_phase_info == 1} {}

  unsigned int chromosome() const { return chromosome_; }
  unsigned int position() const { return position_; }
  bool is_phased() const { return is_phased_; }
  bool mom_is_a() const { return mom_is_a_; }

 private:
  unsigned int position_;
  unsigned char chromosome_;
  bool is_phased_{false};
  bool mom_is_a_{false};
};

class HetTable {
 public:
  enum {
    Obs0 = 0,
    ObsA = 1,
    ObsB = 2,
    ObsAB = 3
  };

  // Load from cache
  HetTable(const std::string & het_table_name,
           const std::string & het_pos_file_name,
           const std::string & het_data_file_name) :
      hets_{het_pos_file_name},
    data_{het_data_file_name},
    n_samples_{data_.size() / n_hets()},
    flip_(n_hets()) {
      std::random_device rd{};
      std::mt19937_64 mersenne{rd()};
      std::function<uint64_t()> gen{std::bind(
          std::uniform_int_distribution<uint64_t>(0, 1), mersenne)};
      for (uint64_t h{0}; h != flip_.size(); ++h) flip_[h] = gen();
      if (1) std::cout << "Loaded " << n_hets() << " hets from "
                       << n_samples() << " samples from file "
                       << het_table_name << std::endl;
    }

  // Load from text input
  explicit HetTable(const std::string & het_table_file_name,
                    const ChromosomeIndexLookup & lookup) {
    const std::string het_base_name{
      remove_substring(het_table_file_name, ".txt")};
    const std::string het_pos_file_name{het_base_name + ".pos.bin"};
    const std::string het_data_file_name{het_base_name + ".data.bin"};
    if (!readable(het_pos_file_name)) {
      std::cout << "Constructing the het info binary cache for "
                << het_table_file_name << std::endl;
      // Input file
      std::ifstream in_stream{het_table_file_name.c_str()};
      if (!in_stream)
        throw Error("Could not open input het file") << het_table_file_name;
      std::vector<HetInfo> temp_hets;
      std::vector<std::vector<unsigned char>> temp_data;
      std::vector<unsigned int> obs_counts;
      std::vector<unsigned int> A_counts;
      std::vector<unsigned int> B_counts;
      std::string line;
      std::string chr_name;
      unsigned int pos;
      int known_phase_value;
      std::string state_string;
      unsigned int n_rows{0};
      unsigned int last_n_cols{0};
      in_stream.ignore(1000000, '\n');
      while (getline(in_stream, line)) {
        ++n_rows;
        std::istringstream line_stream{line.c_str()};
        line_stream >> chr_name >> pos >> known_phase_value;
        if (!line_stream) throw Error("Het info parse error in HetTable");
        const HetInfo het_info{lookup[chr_name], pos, known_phase_value};
        temp_hets.push_back(het_info);
        unsigned int n_cols{0};
        obs_counts.push_back(0);
        A_counts.push_back(0);
        B_counts.push_back(0);
        unsigned int state;
        while (line_stream >> state) {
          if (!known_phase_value) {
            if (state == 1) {
              state = 2;
            } else if (state == 2) {
              state = 1;
            }
          }
          if (temp_data.size() == n_cols) temp_data.resize(n_cols + 1);
          temp_data[n_cols++].push_back(state);
          obs_counts.back() += static_cast<bool>(state);
          A_counts.back() += (state == 1 || state == 3);
          B_counts.back() += (state == 2 || state == 3);
        }
        if (last_n_cols && last_n_cols != n_cols)
          throw Error("Column number discrepancy in HetTable");
        last_n_cols = n_cols;
      }

      // Filter out unusual loci - high counts across samples
      // or very unbalanced in A/B counts
      std::vector<unsigned char> filter_locus(n_rows);
      if (cut_hets) {
        const double max_stdev{4};
        const double p_val_cutoff{1.0 / temp_hets.size()};
        const NormalParams norm{obs_counts};
        unsigned int n_p_cut{0};
        unsigned int n_s_cut{0};
        unsigned int n_cut{0};
        for (uint64_t l{0}; l != temp_hets.size(); ++l) {
          // Note -1 to make less extreme
          const double p_val{gsl_cdf_binomial_Q(
              std::max(A_counts[l], B_counts[l]) - 1, 0.5,
              A_counts[l] + B_counts[l])};
          if (0)
            std::cout << obs_counts[l]
                      << " " << A_counts[l]
                      << " " << B_counts[l]
                      << " " << p_val
                      << std::endl;
          if (p_val < p_val_cutoff) {
            filter_locus[l] = true;
            ++n_p_cut;
          }
          if (obs_counts[l] > norm.mean + max_stdev * norm.stdev) {
            filter_locus[l] = true;
            ++n_s_cut;
          }
          if (filter_locus[l]) ++n_cut;
        }
        std::cout << "Cut a total of " << n_cut
                  << " of " << temp_hets.size() << " het loci: "
                  << n_p_cut << " by binomial test, "
                  << n_s_cut << " by stdev" << std::endl;
      }

      // Het data output
      std::ofstream pos_stream{het_pos_file_name.c_str(), std::ios::binary};
      if (!pos_stream) throw Error("Could not open output het file")
                           << het_pos_file_name;

      for (uint64_t l{0}; l != temp_hets.size(); ++l)
        if (!filter_locus[l])
          pos_stream.write(reinterpret_cast<const char *>(&temp_hets[l]),
                           sizeof(HetInfo));

      std::ofstream data_stream{het_data_file_name.c_str(), std::ios::binary};
      if (!data_stream) throw Error("Could not open output het data file")
                            << het_data_file_name;
      for (uint64_t c{0}; c != last_n_cols; ++c) {
        for (uint64_t l{0}; l != temp_data[c].size(); ++l) {
          if (!filter_locus[l]) {
            data_stream.write(reinterpret_cast<const char *>(&temp_data[c][l]),
                              sizeof(unsigned char));
          }
        }
      }
    }
    new (this) HetTable{het_table_file_name,
          het_pos_file_name, het_data_file_name};
  }

  uint64_t n_hets() const { return hets_.size(); }
  uint64_t n_samples() const { return n_samples_; }
  bool flip(const uint64_t het_index) const { return flip_[het_index]; }
  unsigned char operator()(const uint64_t sample_index,
                           const uint64_t het_index) const {
    if (flip_[het_index]) {
      return flip_obs[data_[sample_index * n_hets() + het_index]];
    } else {
      return data_[sample_index * n_hets() + het_index];
    }
  }
  unsigned char unflipped(const uint64_t sample_index,
                          const uint64_t het_index) const {
    return data_[sample_index * n_hets() + het_index];
  }
  HetInfo operator[](const uint64_t index) const { return hets_[index]; }
  unsigned int max_position() const { return hets_[n_hets() - 1].position(); }

 private:
  PreMappedVector<HetInfo> hets_{"/dev/null", false};
  PreMappedVector<unsigned char> data_{"/dev/null", false};
  uint64_t n_samples_{0};
  std::vector<unsigned char> flip_{};
};

}  // namespace paa

using paa::flip_haplo;
using paa::get_next_file;
using paa::lnsum;
using paa::Bounds;
using paa::ChromosomeIndexLookup;
using paa::CN_abspos;
using paa::Error;
using paa::HetTable;
using paa::Marker;
using paa::Matrix;
using paa::PSDoc;
using paa::PSGraph;
using paa::PSHSeries;
using paa::PSPage;
using paa::PSShade;
using paa::PSXYMSeries;
using paa::PSXYSeries;
using paa::Reference;
using paa::ThreadPool;

// Flip this to test effect of simple afterburner
const bool do_simple_afterburner{true};
// const bool do_afterburner{true};

// Defines data to be returned from the EM over a bin
using ObsType = vector<vector<unsigned char>>;
struct EMResult {
  EMResult() : het_start{0}, haplo_probs{1}, member{1}, failure{true} { }
  explicit EMResult(const uint64_t het_start_,
                    const uint64_t n_loci, const uint64_t n_samples,
                    const unsigned int bin_size_) :
      het_start{het_start_},
      haplo_probs{n_loci},
      member{n_samples},
      bin_size{bin_size_} {
        for (uint64_t s{0}; s != n_samples; ++s) {
          member[s][0] = member[s][1] = member[s][2] = member[s][3] = log(0.25);
        }
        for (uint64_t h{0}; h != n_loci; ++h)
          haplo_probs[h][0] = haplo_probs[h][1] = log(0.5);
      }
  void flip() {
    if (do_simple_afterburner) {
      for (uint64_t sample{0}; sample != haplo_probs.I; ++sample)
        std::swap(haplo_probs[sample][0], haplo_probs[sample][1]);
      for (uint64_t het{0}; het != member.I; ++het)
        std::swap(member[het][1], member[het][2]);
    }
  }
  bool operator<(const EMResult & rhs) const {
    return lnlike < rhs.lnlike;
  }
  const uint64_t het_start;
  Matrix<double, 0, 2> haplo_probs;
  Matrix<double, 0, 4> member;
  double coverage{0.0};
  double background{0.0};
  double lnlike{0.0};
  uint64_t n_iter{0};
  unsigned int bin_size{0};
  bool failure{false};
};

// EM function to run on obs data over a bin defined by start and stop
EMResult em(const ObsType & obs, const uint64_t start, const uint64_t stop,
            const unsigned int bin_size, PSDoc & movie, mutex & movie_mutex,
            const bool do_movie) {
  const uint64_t n_samples{obs.size()};
  const uint64_t n_hets{stop - start};

  // Set up random number generator
  std::random_device rd{};
  std::mt19937_64 mersenne{rd()};

  const uint64_t n_observations{
    [&obs, n_samples, start, stop]() {
      uint64_t result{0};
      for (uint64_t sample{0}; sample != n_samples; ++sample)
        for (uint64_t het{start}; het != stop; ++het)
          result += obs[sample][het] == 0 ? 0 :
              (obs[sample][het] == 3 ? 2 : 1);
      return result;
    }()};

  vector<EMResult> results;
  for (unsigned int trial{0}; trial != n_em_trials; ++trial) {
    // Info to return from function
    EMResult result(start, n_hets, n_samples, bin_size);

    // Set initial coverage and background
    uint64_t allele_seen{0};
    for (uint64_t sample{0}; sample != n_samples; ++sample)
      for (uint64_t het{start}; het != stop; ++het)
        allele_seen += obs[sample][het] > 0;
    result.coverage = 1.0 * allele_seen / obs.size() / n_hets;
    result.background = 0.05;

    // Get push point with good coverage
    vector<uint64_t> good_push_points;
    uint64_t push_thresh{5};
    while (good_push_points.empty()) {
      for (uint64_t h{0}; h != n_hets; ++h) {
        const uint64_t het{h + start};
        uint64_t as{0};
        uint64_t bs{0};
        for (uint64_t sample{0}; sample != n_samples; ++sample) {
          const unsigned char val{obs[sample][het]};
          if (val == 1) ++as;
          if (val == 2) ++bs;
        }
        if (as > push_thresh && bs > push_thresh) good_push_points.push_back(h);
      }
      if (push_thresh == 0) break;
      --push_thresh;
    }
    if (good_push_points.empty()) {
      push_thresh = 5;
      while (good_push_points.empty()) {
        for (uint64_t h{0}; h != n_hets; ++h) {
          const uint64_t het{h + start};
          uint64_t as{0};
          uint64_t bs{0};
          for (uint64_t sample{0}; sample != n_samples; ++sample) {
            const unsigned char val{obs[sample][het]};
            if (val == 1) ++as;
            if (val == 2) ++bs;
          }
          if (as > push_thresh || bs > push_thresh)
            good_push_points.push_back(h);
        }
        if (push_thresh == 0) break;
        --push_thresh;
      }
    }
    if (good_push_points.empty())
      throw Error("No good push points found");

    // Push at random push point, in random direction
    std::function<uint64_t()> push_gen{
      std::bind(std::uniform_int_distribution<uint64_t>(
          0, good_push_points.size() - 1), mersenne)};
    const uint64_t push_point{good_push_points[push_gen()]};
    if (push_gen() % 2) {
      result.haplo_probs[push_point][0] = log(0.95);
      result.haplo_probs[push_point][1] = log(0.05);
    } else {
      result.haplo_probs[push_point][0] = log(0.05);
      result.haplo_probs[push_point][1] = log(0.95);
    }

    // Do EM
    for (; result.n_iter != max_em_iter; ++result.n_iter) {
      // prepare movie
      unique_lock<mutex> lock{movie_mutex};
      static vector<unique_ptr<PSXYMSeries>> frames;
      frames.reserve(max_em_iter);
      static mutex frame_mutexes[max_em_iter];
      if (result.n_iter == frames.size() && do_movie) {
        frames.emplace_back(make_unique<PSXYMSeries>(
            movie, "Haplotypes for iteration " + to_string(result.n_iter) +
            ";H1 alleles;H2 alleles", Bounds{0, max_alleles, 0, max_alleles}));
      }
      lock.unlock();

      // populate movie
      if (do_movie && trial == 0) {
        static thread_local Matrix<double, 0, 2> h_counts(n_samples);
        static thread_local vector<string> colors;
        colors.clear();
        for (uint64_t sample{0}; sample != n_samples; ++sample) {
          for (const bool m2 : {false, true}) h_counts[sample][m2] = 0;
          for (uint64_t h{0}; h != n_hets; ++h) {
            const uint64_t het{h + start};
            for (const bool A : {false, true}) {
              for (const bool m2 : {false, true}) {
                const double val{exp(result.haplo_probs[h][m2])};
                h_counts[sample][m2] += (A ? (1 - val) : val) *
                    (obs[sample][het] == A + 1 || obs[sample][het] == 3);
              }
            }
          }
          ostringstream color;
          for (uint64_t component{1}; component != 4; ++component)
            color << " " << exp(result.member[sample][component]);
          colors.push_back(color.str());
        }
        lock_guard<mutex> frame_lock{frame_mutexes[result.n_iter]};
        PSXYMSeries * frame{frames[result.n_iter].get()};
        for (uint64_t sample{0}; sample != n_samples; ++sample) {
          const Marker marker{
            paa::circle(), 0.5, colors[sample], 0.1, true, colors[sample]};
          frame->add_point(min(max_alleles, h_counts[sample][0]),
                           min(max_alleles, h_counts[sample][1]), marker);
        }
      }

      // simple emission
      Matrix<double, 4, 4> m1_emit{};
      const double c{result.coverage};
      const double cq{1 - c};
      const double cb{c * result.background};
      const double cbq{1 - cb};
      m1_emit = {
        {log(cbq * cbq), log(cb * cbq), log(cbq * cb), log(cb * cb)},
        {log(cq  * cbq), log(c  * cbq), log(cq * cb),  log(c  * cb)},
        {log(cbq * cq),  log(cb * cq),  log(cbq * c),  log(cb * c)},
        {log(cq  * cq),  log(c  * cq),  log(cq  * c),  log(c  * c)}};

      // get member(sample, state) matrix with state probabilities
      double total_lnlike{0};
      for (uint64_t sample{0}; sample != n_samples; ++sample) {
        for (uint64_t state{0}; state != 4; ++state)
          result.member[sample][state] = 0;
        double * mem{result.member[sample]};
        double norm{0};
        for (uint64_t state{0}; state != 4; ++state) {
          for (uint64_t h{0}; h != n_hets; ++h) {
            const uint64_t het{h + start};
            mem[state] += lnsum(
                result.haplo_probs[h][0] +
                m1_emit[state][obs[sample][het]],
                result.haplo_probs[h][1] +
                m1_emit[flip_haplo[state]][obs[sample][het]]);
          }
          norm = state ? lnsum(norm, mem[state]) : mem[state];
        }
        total_lnlike += norm;
        for (uint64_t state{0}; state != 4; ++state) mem[state] -= norm;
      }
      result.lnlike = total_lnlike;

      // Get new haplo_probs(position, model) matrix: prob of model M1 vs M2
      bool all_close{true};
      for (uint64_t h{0}; h != n_hets; ++h) {
        const uint64_t het{h + start};
        double p_obs_model[2]{0, 0};
        double norm{0};
        for (const bool m2 : {false, true}) {
          for (uint64_t sample{0}; sample != n_samples; ++sample) {
            double term{0};
            for (uint64_t state{0}; state != 4; ++state) {
              const uint64_t fstate{m2 ? flip_haplo[state] : state};
              const double value{result.member[sample][state] +
                    m1_emit[fstate][obs[sample][het]]};
              term = state ? lnsum(term, value) : value;
            }
            p_obs_model[m2] += term;
          }
          norm = m2 ? lnsum(norm, p_obs_model[m2]) : p_obs_model[m2];
        }
        for (const bool m2 : {false, true}) {
          p_obs_model[m2] -= norm;
          if (fabs(exp(result.haplo_probs[h][m2]) -
                   exp(p_obs_model[m2])) > 0.000001)
            all_close = false;
        }
        for (const bool m2 : {false, true})
          result.haplo_probs[h][m2] = p_obs_model[m2];
      }

      // Terminate?
      if (all_close) {
        break;
      }
      if (use_old_coverage) {
        // Get new coverage and background
        double opportunity[2][2]{{0, 0}, {0, 0}};
        double observed[2][2]{{0, 0}, {0, 0}};
        const bool COV{0};
        const bool BKG{1};
        for (uint64_t h{0}; h != n_hets; ++h) {
          const uint64_t het{h + start};
          for (uint64_t sample{0}; sample != n_samples; ++sample) {
            for (const bool B : {false, true}) {
              const double cov_opp{exp(lnsum(lnsum(
                  result.member[sample][1 + B] + result.haplo_probs[h][0],
                  result.member[sample][2 - B] + result.haplo_probs[h][1]),
                                             result.member[sample][3]))};
              opportunity[COV][B] += cov_opp;
              const double bkg_opp{exp(lnsum(lnsum(
                  result.member[sample][2 - B] + result.haplo_probs[h][0],
                  result.member[sample][1 + B] + result.haplo_probs[h][1]),
                                             result.member[sample][0]))};
              opportunity[BKG][B] += bkg_opp;
              if (obs[sample][het] == 3 || obs[sample][het] == 1 + B) {
                observed[COV][B] += cov_opp;
                observed[BKG][B] += bkg_opp;
              }
            }
          }
        }
        result.coverage = (observed[COV][0] + observed[COV][1]) /
            (opportunity[COV][0] + opportunity[COV][1]);
        result.background = (observed[BKG][0] + observed[BKG][1]) /
            (opportunity[BKG][0] + opportunity[BKG][1]) / result.coverage;
      } else {
        double ploidy{0.0};
        for (uint64_t sample{0}; sample != n_samples; ++sample)
          ploidy += exp(lnsum(lnsum(lnsum(
              result.member[sample][1], result.member[sample][2]),
                                    result.member[sample][3]),
                              result.member[sample][3]));
        ploidy *= n_hets;
        result.coverage = n_observations / ploidy;
      }
      if (0) cerr << start << " " << result.n_iter
                  << " " << result.coverage
                  << " " << result.background << "\n";
    }
    results.push_back(result);
  }
  if (results.size()) {
    unsigned int max_like{0};
    for (unsigned int r{1}; r != results.size(); ++r) {
      if (results[max_like] < results[r]) max_like = r;
    }
    return results[max_like];
  } else {
    return EMResult();
  }
}

int main(int argc, char * argv[]) try {
  // Process command line arguments
  if (--argc < 5) throw Error("hahaha ref n_loci n_threads down chr ...");
  const Reference ref{argv[1]};
  const ChromosomeIndexLookup chr_lookup{ref};
  const CN_abspos cn_abspos{ref};
  const unsigned int n_loci{static_cast<unsigned int>(atoi(argv[2]))};
  const unsigned int n_threads{static_cast<unsigned int>(atoi(argv[3]))};
  const string downsample{argv[4]};
  ThreadPool pool{n_threads};

  const string bad_loci_name{get_next_file("bad_loci", "txt")};
  ofstream bad_loci{bad_loci_name.c_str()};
  bad_loci << "chr pos bin_size cover" << endl;
  const string all_loci_name{"all_loci.txt"};
  ofstream all_loci{all_loci_name.c_str()};
  all_loci << "chr pos bin_size cover" << endl;

  std::random_device rd{};
  std::mt19937_64 mersenne{rd()};
  std::function<double()> gen{
    std::bind(std::uniform_real_distribution<double>(-0.1, 0.1), mersenne)};
  std::function<double()> gen2{
    std::bind(std::uniform_real_distribution<double>(0, 1), mersenne)};

  PSDoc plots{"hahaha"};
  PSDoc rocs{"roc"};
  PSDoc phasing{"phasing"};
  PSDoc proportions{"proportion"};
  PSPage all_phasing{phasing, "All phasing", "2 11"};
  const Marker black_marker{paa::circle(), 0.1, "0 0 0", 0.1, true};
  PSXYSeries proportion_series{proportions,
        "Chromosome Slicing;% Covering Mom;% Covering Dad",
        Bounds{-0.1, 100.1, -0.1, 100.1}, black_marker};

  PSDoc movie{"movie"};
  mutex movie_mutex;
  PSXYMSeries alleles_plot{movie, "Haplotypes;H1 alleles;H2 alleles",
        Bounds{0, max_alleles, 0, max_alleles}};
  PSGraph alleles_hist_graph{movie, "Haplotypes;H1 alleles;H2 alleles",
        Bounds{0, max_alleles, 0, max_alleles}};
  PSShade alleles_hist{alleles_hist_graph, "1 0 0"};
  ofstream alleles_out{"alleles.txt"};

  using uHist = PSHSeries<uint64_t, uint64_t>;
  using dHist = PSHSeries<double, uint64_t>;
  dHist sample_coverage_plot{plots, "Sample Coverage;Coverage;N",
        Bounds{0.0, 0.2}, 100};
  uHist het_coverage_plot{
    plots, "Het Coverage;Coverage;N", Bounds{0.0, 50.0}, 50};
  dHist cover_hist{
    plots, "EM Coverage per Bin;Coverage;N", Bounds{0, 0.05}, 50};
  dHist background_hist{
    plots, "EM Background per Bin;Background;N", Bounds{0, 0.2}, 50};
  uHist n_iter_hist{plots, "EM iterations;N iterations;N bins",
        100, Bounds{0, 200}};

  using ROC_Point = pair<double, bool>;
  using ROC = vector<ROC_Point>;
  ROC all_roc;
  auto do_roc = [n_loci, &plots, &rocs, &mersenne]
      (ROC & roc, const string & name) {
    shuffle(roc.begin(), roc.end(), mersenne);
    stable_sort(roc.begin(), roc.end());
    PSPage & roc_page{*PSPage::create(plots)};
    rocs.add(&roc_page);
    PSGraph & roc_graph{*PSGraph::create(
        roc_page, "ROC-like curve - " + name + ";% observed;% incorrect",
        Bounds{-0.1, 100.1, 0.001, 20.0})};
    roc_graph.log_y(true);
    PSXYSeries & roc_series{*PSXYSeries::create(roc_graph)};
    double last_confidence{100000000};
    double last_confidence_point{last_confidence};
    uint64_t total_correct{0};
    ostringstream confidence_ps;
    const vector<double> thresholds{0.01, 0.03, 0.1, 0.3, 1.0, 3, 10, 30};
    vector<double> last_incorrect(thresholds.size(), 0.0);

    confidence_ps << "12 sf\n";
    for (uint64_t rr{0}; rr != roc.size(); ++rr) {
      const uint64_t r{roc.size() - rr - 1};
      const ROC_Point & roc_point{roc[r]};
      const double confidence{roc_point.first};
      const bool correct{roc_point.second};
      if (correct) ++total_correct;
      const double incorrect{100 - 100.0 * total_correct / (rr + 1)};
      for (unsigned int nn{0}; nn != thresholds.size(); ++nn) {
        const double threshold{thresholds[nn]};
        if (last_incorrect[nn] > 0 && last_incorrect[nn] <= threshold &&
            incorrect > last_incorrect[nn] && incorrect > threshold) {
          cerr << name << " threshold " << threshold
               << " between last " << last_incorrect[nn] << " and "
               << incorrect
               << " at " << (rr + 1) << " of " << roc.size()
               << " or " << 100.0 * (rr + 1) / roc.size() << "% for "
               << n_loci << " loci bins at confidence " << confidence << endl;
          last_incorrect[nn] = -1;
        }
        if (last_incorrect[nn] >= 0) last_incorrect[nn] = incorrect;
      }
      if (last_confidence > confidence)
        roc_series.add_point(100.0 * rr / roc.size(), incorrect);
      if (fabs(last_confidence_point - confidence) > 1) {
        last_confidence_point = confidence;
        if (incorrect > 0)
          confidence_ps << 100.0 * rr / roc.size()
            << " " << incorrect
            << " gc m ("
            << static_cast<unsigned int>(100 * confidence) / 100
                        << ") show\n";
      }
      last_confidence = confidence;
    }
    roc_graph.ps(confidence_ps.str());
  };

  // Output files
  ofstream em_parameters_file{"em_params.txt"};
  if (!em_parameters_file) throw Error("Problem opening EM parameters");
  ofstream phasing_results_file{"phasing_results.txt"};
  if (!phasing_results_file) throw Error("Problem opening phasing results");
  ofstream well_states_file{"well_states.txt"};
  if (!well_states_file) throw Error("Problem opening well states");

  argc -= 4;
  argv += 4;
  while (argc--) {
    // Get chromosome
    const string chr_name{argv++[1]};
    const unsigned int chr{chr_lookup[chr_name]};

    PSXYSeries & chr_proportion_series{*PSXYSeries::create(
        proportions, "Chromosome Slicing for " + chr_name +
        ";% Covering Mom;% Covering Dad",
        Bounds{-0.1, 100.1, -0.1, 100.1}, black_marker)};

    // Het input data - Note HetTable does locus cuts on load
    const string data_dir{"/data/safe/paa/analysis/haha/ssc_new/data"};
    const string table_name{data_dir + "/" + chr_name + downsample +
          ".obs.txt"};
    const HetTable hets{table_name, chr_lookup};

    // Coverage and background plots, from data and then from EM
    PSPage & phase_page{*PSPage::create(
        plots, "Phasing for " + chr_name, "1 3 (0.3 0.4)")};
    phasing.add(&phase_page);

    // Plot sample and het coverage
    vector<uint64_t> sample_coverage(hets.n_samples());
    vector<uint64_t> het_coverage(hets.n_hets());
    vector<uint64_t> single_het_coverage(hets.n_hets());
    vector<uint64_t> het_good(hets.n_hets());  // UNUSED !!!
    for (uint64_t sample{0}; sample != hets.n_samples(); ++sample) {
      uint64_t sample_cover{0};
      for (uint64_t h{0}; h != hets.n_hets(); ++h) {
        const uint64_t obs{hets(sample, h)};
        const uint64_t count{obs > 0 ? 1U : 0U};
        sample_cover += count;
        het_coverage[h] += count;
        if (obs == 1 || obs == 2) ++single_het_coverage[h];
      }
      sample_coverage[sample] = sample_cover;
      sample_coverage_plot.add_point(1.0 * sample_cover / hets.n_hets());
    }
    for (uint64_t h{0}; h != hets.n_hets(); ++h)
      het_coverage_plot.add_point(het_coverage[h]);

    // Make cuts and put selected data into obs
    ObsType obs;
    for (uint64_t h{0}; h != hets.n_hets(); ++h) {
      het_good[h] = het_coverage[h] >= het_coverage_low_cutoff &&
          het_coverage[h] < het_coverage_high_cutoff;
    }
    for (uint64_t sample{0}; sample != hets.n_samples(); ++sample) {
      if (sample_coverage[sample] >= sample_coverage_low_cutoff &&
          sample_coverage[sample] < sample_coverage_high_cutoff) {
        obs.emplace_back(0);
        for (uint64_t h{0}; h != hets.n_hets(); ++h)
          if (het_good[h]) obs.back().push_back(hets.unflipped(sample, h));
      } else {
        if (false) cerr << "Sample cut " << sample
                        << " " << sample_coverage[sample] << endl;
      }
    }
    const uint64_t n_samples{obs.size()};
    const uint64_t n_hets{obs.front().size()};
    if (hets.n_samples() != n_samples || hets.n_hets() != n_hets)
      cerr << "Cut samples " << hets.n_samples() << " to " << n_samples
           << " and hets " << hets.n_hets() << " to " << n_hets << endl;
    if (hets.n_hets() != n_hets)
      throw Error("Cutting het loci may lead to unexpected bugs");

    // run EM jobs in parallel
    vector<future<EMResult>> em_futures;
    for (uint64_t het_start{0}; het_start < n_hets; het_start += n_loci) {
      const uint64_t het_stop{
        het_start + n_loci <= n_hets ? het_start + n_loci : n_hets};
      em_futures.push_back(pool.run(em, std::ref(obs), het_start, het_stop,
                                    hets[het_stop].position() -
                                    hets[het_start].position(),
                                    std::ref(movie), std::ref(movie_mutex),
                                    het_start == 10 * n_loci &&
                                    chr_name == "chr17"));
    }

    // Get results
    vector<EMResult> results;
    vector<unsigned int> bin_sizes;
    vector<vector<double>> haplo_probs(2);
    vector<vector<double>> member(4);
    vector<uint64_t> n_iter;
    double total_like{0.0};
    uint64_t bin_number{0};
    vector<uint64_t> bin_numbers;
    em_parameters_file
        << "bin\thet_start\tchr\tpos\tcover\tbackground"
        << "\tlnlike\tn_iter\tfailure" << endl;
    for (future<EMResult> & future : em_futures) {
      EMResult result{future.get()};
      if (result.failure) continue;
      total_like += result.lnlike;
      cover_hist.add_point(result.coverage);
      background_hist.add_point(result.background);
      em_parameters_file << bin_number
                         << '\t' << result.het_start
                         << '\t'
                         << ref.name(hets[result.het_start].chromosome())
                         << '\t' << hets[result.het_start].position()
                         << '\t' << result.coverage
                         << '\t' << result.background
                         << '\t' << result.lnlike
                         << '\t' << result.n_iter
                         << '\t' << result.failure
                         << endl;

      n_iter.push_back(result.n_iter);
      if (results.size()) {
        int do_flip{0};
        // See if we should flip this result to match orientation of last few
        for (unsigned int i{0};
             i != min(afterburner_range, results.size()); ++i) {
          const EMResult & last{results[results.size() - 1 - i]};
          double A{0};
          double B{0};
          const uint64_t h1{1};
          const uint64_t h2{2};
          for (uint64_t sample{0}; sample != n_samples; ++sample) {
            A += exp(result.member[sample][h1] + last.member[sample][h1]) +
                exp(result.member[sample][h2] + last.member[sample][h2]);
            B += exp(result.member[sample][h1] + last.member[sample][h2]) +
                exp(result.member[sample][h2] + last.member[sample][h1]);
          }
          do_flip += B > A ? 1 : -1;
        }
        if (do_flip > 0) result.flip();
      }
      // plot alleles
      if (result.failure == false) {
        for (uint64_t sample{0}; sample != n_samples; ++sample) {
          double h_counts[2]{0, 0};
          for (uint64_t h{0}; h != result.haplo_probs.I; ++h) {
            const uint64_t het{h + result.het_start};
            for (const bool A : {false, true}) {
              for (const bool m2 : {false, true}) {
                const double val{exp(result.haplo_probs[h][m2])};
                h_counts[m2] += (A ? (1 - val) : val) *
                    (obs[sample][het] == A + 1 || obs[sample][het] == 3);
              }
            }
          }
          ostringstream color;
          color.precision(6);
          color << fixed;
          for (uint64_t component{1}; component != 4; ++component) {
            if (component > 1) color << " ";
            color << max(0.000001, exp(result.member[sample][component]));
          }
          const Marker marker{paa::circle(), 0.2, color.str(), 0.1,
                true, color.str()};
          const double x{min(max_alleles, h_counts[0] + gen2())};
          const double y{min(max_alleles, h_counts[1] + gen2())};
          alleles_out << ref.name(hets[result.het_start].chromosome())
                      << '\t' << hets[result.het_start].position()
                      << '\t' << sample
                      << '\t' << x << '\t' << y << '\t' << color.str() << '\n';
          alleles_plot.add_point(x, y, marker);
          alleles_hist.add_point(x, y);
        }
      }

      // save results for chromosome in single data structure
      for (uint64_t h{0}; h != result.haplo_probs.I; ++h) {
        for (const bool m2 : {false, true}) {
          haplo_probs[m2].push_back(1 - exp(result.haplo_probs[h][m2]));
        }
        bin_sizes.push_back(result.bin_size);
        bin_numbers.push_back(bin_number);
      }
      for (uint64_t sample{0}; sample != result.member.I; ++sample)
        for (const uint64_t state : {0, 1, 2, 3})
          member[state].push_back(1 - exp(result.member[sample][state]));
      results.push_back(std::move(result));
      ++bin_number;
    }
    for (uint64_t sample{0}; sample != n_samples; ++sample) {
      double mom_total{0};
      double dad_total{0};
      for (const uint64_t state : {1, 2, 3}) {
        for (uint64_t bin{0}; bin != results.size(); ++bin) {
          if (state == 1 || state == 3)
            mom_total += exp(results[bin].member[sample][state]);
          if (state == 2 || state == 3)
            dad_total += exp(results[bin].member[sample][state]);
        }
      }
      chr_proportion_series.add_point(100.0 * mom_total / results.size(),
                                      100.0 * dad_total / results.size());
      proportion_series.add_point(100.0 * mom_total / results.size(),
                                  100.0 * dad_total / results.size());
    }

    sort(n_iter.begin(), n_iter.end());
    for (const uint64_t iter : n_iter)
      n_iter_hist.add_point(std::min(iter, static_cast<uint64_t>(99)));

    cerr << "Total log likelihood = " << total_like << ", prob per obs = "
         << exp(total_like / hets.n_hets() / hets.n_samples()) << endl;
    // Flip results globally?
    {
      uint64_t n{0};
      double sum{0.0};
      for (uint64_t h{0}; h != haplo_probs[0].size(); ++h) {
        if (-log10(0.5 - fabs(haplo_probs[0][h] - 0.5)) >= 5) {
          ++n;
          sum += haplo_probs[0][h];
        }
      }
      if (sum / n < 0.5)
        haplo_probs[0].swap(haplo_probs[1]);
    }
    phasing_results_file << "het\tbin\tchr\tpos\tp_mom_a\tis_phased\tphase"
                         << endl;
    for (uint64_t h{0}; h != haplo_probs[0].size(); ++h) {
      phasing_results_file
          << h
          << '\t' << bin_numbers[h]
          << '\t' << ref.name(hets[h].chromosome())
          << '\t' << hets[h].position()
          << '\t' << haplo_probs[0][h]
          << '\t' << hets[h].is_phased()
          << '\t' << (hets[h].is_phased() ?
                      (hets[h].mom_is_a() ? 'M' : 'F') : '-')
          << endl;
    }
    bin_number = 0;
    well_states_file << "bin\thet_start\tchr\tpos";
    const string state_names{"0MFB"};
    for (uint64_t well{0}; well != results.front().member.I; ++well) {
        for (uint64_t state{0}; state != 4; ++state) {
          well_states_file << '\t' << to_string(well) << state_names[state];
        }
    }
    well_states_file << endl;
    for (const EMResult & result : results) {
      well_states_file
          << bin_number++
          << '\t' << result.het_start
          << '\t'
          << ref.name(hets[result.het_start].chromosome())
          << '\t' << hets[result.het_start].position();
      for (uint64_t well{0}; well != result.member.I; ++well) {
        for (uint64_t state{0}; state != 4; ++state) {
          well_states_file << '\t' << exp(result.member[well][state]);
        }
      }
      well_states_file << endl;
    }

    // Roc-like plot
    ROC roc;
    for (uint64_t h{0}; h != haplo_probs[0].size(); ++h) {
      const double confidence{
        -log10(0.5 - fabs(haplo_probs[0][h] - 0.5))};
      if (hets[h].is_phased() && confidence > 0.000001) {
        if (haplo_probs[0][h] < 0.01) {
          bad_loci << chr_name
                   << " " << hets[h].position()
                   << " " << bin_sizes[h]
                   << " " << single_het_coverage[h] << endl;
        }
        if (haplo_probs[0][h] > 0.99) {
          all_loci << chr_name
                   << " " << hets[h].position()
                   << " " << bin_sizes[h]
                   << " " << single_het_coverage[h] << endl;
        }
        roc.emplace_back(confidence, !(haplo_probs[0][h] <= 0.5));
        all_roc.emplace_back(confidence, !(haplo_probs[0][h] <= 0.5));
      }
    }
    do_roc(roc, chr_name);

    string colors[4]{"0 0 0", "1 0 0", "0 0 1", "0 1 0"};
    PSGraph & state_graph{*PSGraph::create(
        plots, "Well-Bin State Probabilities;Well Rank;1 - prob")};
    state_graph.log_y(true);
    for (const uint64_t state : {0, 1, 2, 3}) {
      sort(member[state].begin(), member[state].end());
      const Marker marker{paa::circle(), 0.1, colors[state], 0.1, true};
      PSXYSeries * series{PSXYSeries::create(state_graph, marker)};
      for (uint64_t sample{0}; sample != member[state].size(); ++sample)
        series->add_point(sample, member[state][sample]);
    }

    PSGraph & all_prob_graph{*PSGraph::create(
        phase_page, ";Position;Prob Mom is A",
        Bounds{0.0, 1.0 * hets[hets.n_hets() - 1].position(), -0.1, 1.1})};
    PSXYSeries & all_prob_series{*PSXYSeries::create(all_prob_graph)};

    dHist & prob_hist{*dHist::create(
        phase_page, ";Probability Mom is A;N", Bounds{-0.1, 1.1}, 100u)};

    PSGraph & prob_graph{*PSGraph::create(
        phase_page, ";Position;Prob Mom is A",
        Bounds{0.0, 1.0 * hets[hets.n_hets() - 1].position(), -0.1, 1.1})};
    all_phasing.add(&prob_graph);
    PSXYSeries & prob_series{*PSXYSeries::create(prob_graph)};
    ofstream prob_output{"probs.txt"};
    for (uint64_t h{0}; h != haplo_probs[0].size(); ++h) {
      if (hets[h].is_phased()) {
        prob_hist.add_point(haplo_probs[0][h]);
        all_prob_series.add_point(hets[h].position(),
                                  haplo_probs[0][h] + gen());
        if (-log10(0.5 - fabs(haplo_probs[0][h] - 0.5)) >= 5)
          prob_series.add_point(hets[h].position(), haplo_probs[0][h] + gen());
      }
      prob_output << chr_name << " " << hets[h].position()
                  << " " << cn_abspos(chr, hets[h].position())
                  << " " << haplo_probs[0][h] << '\n';
    }

    PSGraph & phase_graph{*PSGraph::create(
        plots, "Phasing confidence;Locus Rank;1 - prob")};
    phase_graph.log_y(true);
    for (const bool m2 : {false, true}) {
      sort(haplo_probs[m2].begin(), haplo_probs[m2].end());
      const Marker marker{paa::circle(),
            0.1, m2 ? "1 0 0" : "0 0 1", 0.1, true};
      PSXYSeries * series{PSXYSeries::create(phase_graph, marker)};
      for (uint64_t h{0}; h != haplo_probs[m2].size(); ++h)
        series->add_point(m2 ? h : haplo_probs[m2].size() - h - 1,
                          haplo_probs[m2][h]);
    }
  }
  do_roc(all_roc, "All Chromosomes");
  cerr << "Exit" << endl;

  return 0;
} catch (Error & e) {
  cerr << "paa::Error:" << endl;
  cerr << e.what() << endl;
  return 1;
} catch (exception & e) {
  cerr << "std::exception" << endl;
  cerr << e.what() << endl;
  return 1;
} catch (...) {
  cerr << "unknown exception was caught" << endl;
  return 1;
}
