//
// haha.h
//
// Code used to process HAHA data
//
// Copyright Peter Andrews 2018 CSHL
//

#ifndef PAA_HAHA_H_
#define PAA_HAHA_H_

#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>

#include <algorithm>
#include <array>
#include <fstream>
#include <functional>
#include <initializer_list>
#include <limits>
#include <map>
#include <mutex>
#include <random>
#include <string>
#include <utility>
#include <vector>

#include "files.h"
#include "lowess.h"
#include "mumdex.h"
#include "psplot.h"
#include "stats.h"
#include "threads.h"

namespace paa {

constexpr bool log_paranoia{false};
constexpr bool cut_hets{true};
constexpr bool cut_loci{true};
constexpr bool cut_samples{false};

template <class FLOAT>
inline FLOAT lnproduct(const FLOAT lhs, const FLOAT rhs) {
  if (!log_paranoia) return lhs + rhs;
  if (lhs > std::numeric_limits<FLOAT>::lowest() &&
      rhs > std::numeric_limits<FLOAT>::lowest()) {
    return lhs + rhs;
  } else {
    return std::numeric_limits<FLOAT>::lowest();
  }
}

template <class FLOAT>
inline FLOAT lnproduct(const FLOAT a1, const FLOAT a2, const FLOAT a3) {
  return lnproduct(a1, lnproduct(a2, a3));
}

template <class FLOAT>
inline FLOAT lnsum(const FLOAT lhs, const FLOAT rhs) {
  if (log_paranoia) {
    if (lhs <= std::numeric_limits<FLOAT>::lowest())
      return rhs;
    if (rhs <= std::numeric_limits<FLOAT>::lowest())
      return lhs;
  }
  if (lhs > rhs)
    return lhs + log(1 + exp(rhs - lhs));
  return rhs + log(1 + exp(lhs - rhs));
}

template <class Vec>
inline double lnsum(const Vec & vec, const uint64_t n) {
  double result{vec[0]};
  for (uint64_t i{1}; i != n; ++i) result = lnsum(vec[i], result);
  return result;
}

template <class Val, uint64_t I>
using Array = std::array<Val, I>;

// Matrix with two compile-time dimensions
template <class Val, uint64_t II = 0, uint64_t JJ = 0>
class Matrix {
 public:
  static constexpr uint64_t I{II};
  static constexpr uint64_t J{JJ};
  Matrix() {}  // Uninitialized!
  Matrix(const std::initializer_list<std::initializer_list<Val> > & input) {
    auto rows = input.begin();
    for (uint64_t i{0}; i != I; ++i, ++rows) {
      auto val = rows->begin();
      for (uint64_t j{0}; j != J; ++j, ++val) (*this)[i][j] = *val;
    }
  }
  const Val * operator[](const uint64_t i) const { return data + i * J; }
  Val * operator[](const uint64_t i) { return data + i * J; }

 private:
  Val data[I * J];
};

// Base class for dynamic matrices
template <class Val> struct MBase {
  explicit MBase(const uint64_t n_elem) : data(n_elem) { }
  std::vector<Val> data;
};
// Matrix with two run-time dimensions
template <class Val>
class Matrix<Val, 0, 0> : private MBase<Val> {
 public:
  const uint64_t I;
  const uint64_t J;
  Matrix(const uint64_t I_, const uint64_t J_) :
      MBase<Val>{I_ * J_}, I{I_}, J{J_} { }
  const Val * operator[](const uint64_t i) const { return &this->data[i * J]; }
  Val * operator[](const uint64_t i) { return &this->data[i * J]; }
};
// Matrix with run-time first dimension and compile-time second dimension
template <class Val, uint64_t JJ>
class Matrix<Val, 0, JJ> : private MBase<Val> {
 public:
  const uint64_t I;
  static constexpr uint64_t J{JJ};
  explicit Matrix(const uint64_t I_) : MBase<Val>{I_ * J}, I{I_} { }
  const Val * operator[](const uint64_t i) const { return &this->data[i * J]; }
  Val * operator[](const uint64_t i) { return &this->data[i * J]; }
  template <class Array> void assign(const Array & a) {
    for (uint64_t i{0}; i != I; ++i)
      for (uint64_t j{0}; j != J; ++j)
        this->data[i * J + j] = a[j];
  }
};
// Matrix with compile-time first dimension and run-time second dimension
template <class Val, uint64_t II>
class Matrix<Val, II, 0> : private MBase<Val> {
 public:
  static constexpr uint64_t I{II};
  const uint64_t J;
  explicit Matrix(const uint64_t J_) : MBase<Val>{I * J_}, J{J_} { }
  const Val * operator[](const uint64_t i) const { return &this->data[i * J]; }
  Val * operator[](const uint64_t i) { return &this->data[i * J]; }
};
// Output for matrices
template <class Val, uint64_t I, uint64_t J>
std::ostream & operator<<(std::ostream & stream,
                        const Matrix<Val, I, J> & matrix) {
  for (uint64_t i{0}; i != matrix.I; ++i) {
    for (uint64_t j{0}; j != matrix.J; ++j) {
      if (j) stream << '\t';
      stream << matrix[i][j];
    }
    if (i + 1 != matrix.I) stream << '\n';
  }
  return stream;
}

// Agrees with eigenvalue-based general method
inline Matrix<double, 2, 2> matrix_power(const Matrix<double, 2, 2> & unit,
                                         const uint64_t power) {
  const double a{unit[0][1]};
  const double b{unit[1][0]};
  const double a_plus_b{a + b};
  const double gamma{pow(1 - a - b, power)};
  return Matrix<double, 2, 2>{
    {(a * gamma + b) / a_plus_b, (a - a * gamma) / a_plus_b},
    {(b - b * gamma) / a_plus_b, (b * gamma + a) / a_plus_b}};
}

#if 0
#include <Eigen/Eigenvalues>
template <int NRows>
class MatrixPowers {
 public:
  using Internal = Eigen::Matrix<double, NRows, NRows>;
  MatrixPowers() {}

  template <class Matrix>
  explicit MatrixPowers(const Matrix & original) {
    initialize(original);
  }

  template <class Matrix>
  void initialize(const Matrix & original) {
    if (original.size() != NRows)  // || original[0].size() != NRows)
      throw Error("Bad Matrix size in MatrixPowers");
    for (unsigned int i{0}; i != NRows; ++i) {
      for (unsigned int j{0}; j != NRows; ++j) {
        input(i, j) = original[i][j];
      }
    }
    std::cerr << "Input matrix:" << std::endl
              << input << std::endl;
    solver.compute(input);
    std::cerr << "Eigenvalues:" << std::endl
              << (eigenvalues = solver.eigenvalues().real()) << std::endl;
    std::cerr << "Eigenvector matrix:" << std::endl
              << (eigenvectors = solver.eigenvectors().real()) << std::endl;
    eigenvalues = solver.eigenvalues().real();
    inverse_eigenvectors = solver.eigenvectors().real().inverse();
  }

  Internal operator()(const double d) const {
    Internal result{Internal::Zero()};
    for (unsigned int r{0}; r != eigenvalues.size(); ++r) {
      result(r, r) = pow(eigenvalues(r), d);
    }
    return eigenvectors * result * inverse_eigenvectors;
  }

  static void test() {
    // Test Matrix powers
    const double a{0.00001};
    const double b{0.00005};
    const Matrix2d<2, 2> trans{{1 - a, a}, {b, 1 - b}};
    const paa::MatrixPowers<2> powers{trans};
    for (const double power : {0, 1, 10, 100, 10000, 100000, 1000000}) {
      std::cout << "Power is " << power << ":" << std::endl;
      std::cout << powers(power) << std::endl;
      std::cout << "new:" << std::endl;
      std::cout << matrix_power(trans, power) << std::endl;
      std::cout << std::endl;
    }
  }

 private:
  Eigen::Matrix<double, NRows, 1> eigenvalues{};
  Internal input{};
  Internal eigenvectors{};
  Internal inverse_eigenvectors{};
  Eigen::EigenSolver<Internal> solver{};
};
#endif

class FragmentInfo {
 public:
  FragmentInfo(const unsigned int chromosome__,
               const unsigned int start__,
               const unsigned int stop__) :
      chromosome_{chromosome__},
    start_{start__},
    stop_{stop__} { }

  unsigned int chromosome() const { return chromosome_; }
  unsigned int start() const { return start_; }
  unsigned int stop() const { return stop_; }
  unsigned int length() const { return stop_ - start_; }

 private:
  uint64_t chromosome_ : 8;
  uint64_t start_ : 28;
  uint64_t stop_ : 28;
};

class FragmentTable {
 public:
  using CountType = uint8_t;
  // Load from cache
  FragmentTable(const std::string & fragment_table_file_name,
                const std::string & all_fragment_pos_file_name,
                const std::string & fragment_pos_file_name,
                const std::string & fragment_data_file_name) :
      all_fragments_{all_fragment_pos_file_name},
      fragments_{fragment_pos_file_name},
    data_{fragment_data_file_name} {
      n_samples_ = data_.size() / n_fragments();
      std::cout << "Loaded " << n_fragments() << " fragments from "
                << n_samples() << " samples from file "
                << fragment_table_file_name << std::endl;
    }

  // Load from text input
  explicit FragmentTable(const std::string & fragment_table_file_name,
                         const ChromosomeIndexLookup & lookup,
                         const unsigned int max_count = 2) {
    const std::string fragment_base_name{
      remove_substring(fragment_table_file_name, ".txt")};
    const std::string all_fragment_pos_file_name{
      fragment_base_name + ".allpos.bin"};
    const std::string fragment_pos_file_name{fragment_base_name + ".pos.bin"};
    const std::string fragment_data_file_name{fragment_base_name + ".data.bin"};
    if (!readable(fragment_pos_file_name)) {
      std::cout << "Constructing fragment info binary cache for "
                << fragment_table_file_name << std::endl;
      // Input file
      std::ifstream in_stream{fragment_table_file_name.c_str()};
      if (!in_stream) throw Error("Could not open input fragment file")
                          << fragment_table_file_name;

      std::vector<FragmentInfo> temp_fragments;
      std::vector<std::vector<CountType>> temp_data;
      std::vector<unsigned int> locus_count_sums;
      std::vector<unsigned int> sample_count_sums;
      std::string line;
      std::string chr_name;
      unsigned int id;
      unsigned int pos_start;
      unsigned int pos_stop;
      unsigned int count;
      unsigned int n_rows{0};
      unsigned int last_n_cols{0};
      while (getline(in_stream, line)) {
        ++n_rows;
        std::istringstream line_stream{line.c_str()};
        line_stream >> chr_name >> id >> pos_start >> pos_stop;
        if (!line_stream)
          throw Error("Fragment info parse error in FragmentTable");
        const FragmentInfo fragment_info{lookup[chr_name], pos_start, pos_stop};
        temp_fragments.push_back(fragment_info);
        unsigned int n_cols{0};
        unsigned int locus_count{0};
        while (line_stream >> count) {
          if (temp_data.size() == n_cols) {
            temp_data.resize(n_cols + 1);
            sample_count_sums.resize(n_cols + 1);
          }
          const unsigned int mcount{std::min(count, max_count)};
          sample_count_sums[n_cols] += mcount;
          temp_data[n_cols++].push_back(mcount);
          locus_count += mcount;
        }
        locus_count_sums.push_back(locus_count);
        if (last_n_cols && last_n_cols != n_cols)
          throw Error("Column number discrepancy in FragmentTable");
        last_n_cols = n_cols;
      }

      // Filter out unusual loci - low or high counts across samples
      std::vector<unsigned char> filter_locus(n_rows);
      if (cut_loci) {
        const double max_stdev{4};
        const unsigned int min_count{4};
        const NormalParams norm{locus_count_sums};
        unsigned int n_cut{0};
        unsigned int n_s_cut{0};
        unsigned int n_m_cut{0};
        for (uint64_t f{0}; f != temp_fragments.size(); ++f) {
          if (locus_count_sums[f] < min_count) {
            filter_locus[f] = true;
            ++n_m_cut;
          }
          if (locus_count_sums[f] > norm.mean + max_stdev * norm.stdev) {
            filter_locus[f] = true;
            ++n_s_cut;
          }
          if (filter_locus[f]) ++n_cut;
        }
        std::cout << "Cut a total of " << n_cut
                  << " of " << temp_fragments.size() << " fragment loci: "
                  << n_m_cut << " for low count, "
                  << n_s_cut << " by stdev" << std::endl;
      }

      // Filter out bad samples
      std::vector<unsigned char> filter_sample(sample_count_sums.size());
      if (cut_samples) {
        const double max_stdev{1.5};
        const NormalParams norm{sample_count_sums};
        unsigned int n_samples_cut{0};
        unsigned int n_samples_low_cut{0};
        unsigned int n_samples_high_cut{0};
        for (uint64_t sample{0}; sample != filter_sample.size(); ++sample) {
          if (sample_count_sums[sample] < norm.mean - max_stdev * norm.stdev) {
            ++n_samples_low_cut;
            filter_sample[sample] = true;
          }
          if (sample_count_sums[sample] > norm.mean + max_stdev * norm.stdev) {
            ++n_samples_high_cut;
            filter_sample[sample] = true;
          }
          if (filter_sample[sample]) ++n_samples_cut;
        }
        std::cout << "Cut a total of " << n_samples_cut
                  << " of " << sample_count_sums.size() << " samples: "
                  << n_samples_low_cut << " for low count, "
                  << n_samples_high_cut << " for high count" << std::endl;
      }

      // Fragment data output
      std::ofstream all_pos_stream{all_fragment_pos_file_name.c_str(),
            std::ios::binary};
      if (!all_pos_stream) throw Error("Could not open output fragment file")
                               << all_fragment_pos_file_name;
      std::ofstream pos_stream{fragment_pos_file_name.c_str(),
            std::ios::binary};
      if (!pos_stream) throw Error("Could not open output fragment file")
                           << fragment_pos_file_name;
      for (uint64_t f{0}; f != temp_fragments.size(); ++f) {
        all_pos_stream.write(reinterpret_cast<const char *>(&temp_fragments[f]),
                         sizeof(FragmentInfo));
        if (!filter_locus[f]) {
          pos_stream.write(reinterpret_cast<const char *>(&temp_fragments[f]),
                           sizeof(FragmentInfo));
        }
      }

      // Count data output
      std::ofstream data_stream{fragment_data_file_name.c_str(),
            std::ios::binary};
      if (!data_stream) throw Error("Could not open output fragment data file")
                            << fragment_data_file_name;
      for (uint64_t c{0}; c != last_n_cols; ++c) {
        if (filter_sample[c]) continue;
        for (uint64_t f{0}; f != temp_data[c].size(); ++f) {
          if (!filter_locus[f]) {
            data_stream.write(reinterpret_cast<const char *>(&temp_data[c][f]),
                              sizeof(CountType));
          }
        }
      }
    }
    new (this) FragmentTable{fragment_table_file_name,
          all_fragment_pos_file_name, fragment_pos_file_name,
          fragment_data_file_name};
  }

  // uint64_t size() const { return fragments_.size(); }
  uint64_t n_fragments() const { return fragments_.size(); }
  uint64_t n_samples() const { return n_samples_; }
  unsigned char operator()(const uint64_t sample_index,
                           const uint64_t fragment_index) const {
    return data_[sample_index * n_fragments() + fragment_index];
  }
  FragmentInfo operator[](const uint64_t fragment_index) const {
    return fragments_[fragment_index];
  }
  unsigned int max_position() const {
    return fragments_[n_fragments() - 1].start();
  }

  std::vector<unsigned int> get_distances() const {
     // Get locus distances
    std::vector<unsigned int> result;
    result.reserve(n_fragments());
    result.push_back(0);
    unsigned int last_position{fragments_[0].start()};
    for (uint64_t f{1}; f != n_fragments(); ++f) {
      // Get distances between loci
      const unsigned int position{fragments_[f].start()};
      result.push_back(position - last_position);
      if (position <= last_position) throw Error("Fragment ordering problem");
      last_position = position;
    }
    return result;
  }

  const PreMappedVector<FragmentInfo> & all_fragments() const {
    return all_fragments_;
  }

 private:
  PreMappedVector<FragmentInfo> all_fragments_{"/dev/null", false};
  PreMappedVector<FragmentInfo> fragments_{"/dev/null", false};
  PreMappedVector<CountType> data_{"/dev/null", false};
  uint64_t n_samples_{0};
};

std::vector<double> coverage_lowess(const std::string & title,
                                    const FragmentTable & fragment_table,
                                    PSDoc & plots) {
  // Assemble data for coverage LOWESS
  std::vector<unsigned int> combined_coverage(fragment_table.n_fragments());
  std::vector<unsigned int> fragment_lengths(fragment_table.n_fragments());
  unsigned int min_coverage{std::numeric_limits<unsigned int>::max()};
  for (uint64_t f{0}; f != fragment_table.n_fragments(); ++f) {
    for (uint64_t s{0} ; s != fragment_table.n_samples(); ++s) {
      combined_coverage[f] += fragment_table(s, f);
    }
    if (combined_coverage[f] > 0 && combined_coverage[f] < min_coverage)
      min_coverage = combined_coverage[f];
    fragment_lengths[f] = fragment_table[f].length();
  }

  // Random numbers for fuzzing plots
  std::random_device rd{};
  std::mt19937_64 mersenne{rd()};
  std::function<double()> gen{std::bind(
      std::uniform_real_distribution<double>(0, 1), mersenne)};

  // Compute LOWESS correction factors
  std::vector<double> smoothed_coverage{
    lowess_correction(fragment_lengths, combined_coverage)};
  for (double & coverage : smoothed_coverage)
    if (coverage < min_coverage) coverage = min_coverage;

  // Create and fill plots
  using Hist = PSHSeries<unsigned int, unsigned int>;
  Hist * const fragments_hist{new Hist{plots,
          title + " fragment coverage;Position;N",
          Bounds{0, 1.0 * fragment_table.max_position()}}};
  ownp(fragments_hist);
  Hist * const samples_hist{new Hist{plots,
          title + " sample counts;Sample index;N",
          Bounds{0, 1.0 * fragment_table.n_samples()},
          static_cast<unsigned int>(fragment_table.n_samples())}};
  ownp(samples_hist);
  const Marker black_marker{paa::circle(), 0.1, "0 0 0", 0.1, true};
  const Marker red_marker{paa::circle(), 0.1, "1 0 0", 0.1, true};
  std::unique_ptr<PSXYSeries> coverage_series{std::make_unique<PSXYSeries>(
      plots, title + " fragment Coverage;Position;Combined Coverage",
      black_marker)};
  PSXYSeries::PSPart & coverage_graph{coverage_series->graph()};
  coverage_graph.log_y(true);
  std::unique_ptr<PSXYSeries> lowess_series{std::make_unique<PSXYSeries>(
      plots, title + " fragment Coverage;Fragment Length;Combined Coverage",
      black_marker)};
  PSXYSeries::PSPart & lowess_graph{lowess_series->graph()};
  lowess_graph.log_y(true);
  std::unique_ptr<PSXYSeries> smoothed_series{std::make_unique<PSXYSeries>(
      lowess_graph, red_marker)};
  std::unique_ptr<PSXYSeries> smoothed_coverage_series{
    std::make_unique<PSXYSeries>(
      plots, title + " smoothed Fragment Coverage;Position;Smoothed Coverage",
      black_marker)};
  PSXYSeries::PSPart & smoothed_coverage_graph{
    smoothed_coverage_series->graph()};
  smoothed_coverage_graph.log_y(true);

  for (unsigned int sample{0}; sample != fragment_table.n_samples(); ++sample) {
    for (uint64_t frag{0}; frag != fragment_table.n_fragments(); ++frag) {
      const uint64_t count{fragment_table(sample, frag)};
      for (uint64_t c{0}; c != count; ++c) {
        samples_hist->add_point(sample);
      }
    }
  }
  for (uint64_t f{0}; f != fragment_table.n_fragments(); ++f) {
    fragments_hist->add_point(fragment_table[f].start());
    coverage_series->add_point(fragment_table[f].start(),
                               combined_coverage[f] + gen());
    lowess_series->add_point(fragment_table[f].length(),
                             combined_coverage[f] + gen());
    smoothed_series->add_point(fragment_table[f].length(),
                               smoothed_coverage[f]);
    smoothed_coverage_series->add_point(fragment_table[f].start(),
                                        combined_coverage[f] /
                                        smoothed_coverage[f]);
  }
  coverage_graph.own(std::move(coverage_series));
  smoothed_coverage_graph.own(std::move(smoothed_coverage_series));
  lowess_graph.own(std::move(lowess_series));
  lowess_graph.own(std::move(smoothed_series));
  return smoothed_coverage;
}

using uArray3 = Array<unsigned int, 3>;
constexpr const uArray3 copies{{0, 1, 2}};
constexpr const uArray3 counts{{0, 1, 2}};

class CoverageHMM {
 public:
  static constexpr double SMALL{std::numeric_limits<double>::lowest()};
  static constexpr uint64_t nCopies{copies.size()};
  static constexpr uint64_t nObservations{counts.size()};
  static_assert(nCopies == 3, "nCopies check");
  static_assert(nObservations == 3, "nObservations check");

  using Trans = Matrix<double, nCopies, nCopies>;
  using Emit = Matrix<double, nCopies, nObservations>;
  using StateProbs = Array<double, nCopies>;
  // using StateChoices = unsigned char[nCopies];  // fails to compile on mac
  using StateChoices = Array<unsigned char, nCopies>;

  CoverageHMM(const std::string & title_,
              const FragmentTable & fragments_,
              ThreadPool & pool_,
              PSDoc & plots_,
              const double init0_,
              const double background_,
              const double fragment_length_) :
      title{title_},
      fragments{fragments_},
    nSamples{fragments.n_samples()},
    nFragments{fragments.n_fragments()},
    pool{pool_},
    locus_distances{fragments.get_distances()},
    init0{init0_},
    init1{1 - init0},
    background{background_},
    fragment_length{fragment_length_},
    plots{plots_},
    coverages{coverage_lowess(title, fragments, plots)},
    transitions(nFragments),
    emissions(nFragments),
    viterbis(nSamples, nFragments),
    viterbi_probs(nSamples),
    segment_starts(nSamples) {
      initialize();
    }

  void initialize() {
    for (double & coverage : coverages) coverage /= (2 * nSamples);
    update_model();
    show_hmm();
  }
  void update_model() {
    initial_state.assign(nSamples, calculate_initial_state());
    auto fut = pool.run([this]() { calculate_transitions(); });
    calculate_emissions();
    fut.get();
  }
  StateProbs calculate_initial_state() const {
    return {{log(init0 * init0),
            log(2 * init0 * init1), log(init1 * init1)}};
  }
  void calculate_transitions() {
    const double rate10{1.0 / fragment_length};
    const double rate01{rate10 * init0 / init1};
    const Matrix<double, 2, 2> strandT{
      {1 - rate01, rate01}, {rate10, 1 - rate10}};
    for (uint64_t f{1}; f != nFragments; ++f) {
      const Matrix<double, 2, 2> strandP{
        matrix_power(strandT, locus_distances[f])};
      transitions[f] = {
        {log(strandP[0][0] * strandP[0][0]),
         log(2 * strandP[0][0] * strandP[0][1]),
         log(strandP[0][1] * strandP[0][1])},
        {log(strandP[1][0] * strandP[0][0]),
         log(strandP[1][1] * strandP[0][0] + strandP[1][0] * strandP[0][1]),
         log(strandP[1][1] * strandP[0][1])},
        {log(strandP[1][0] * strandP[1][0]),
         log(2 * strandP[1][0] * strandP[1][1]),
         log(strandP[1][1] * strandP[1][1])}};
    }
    if (0) show_transitions();
  }
  void calculate_emissions() {
    for (uint64_t f{0}; f != nFragments; ++f) {
      const double est_rate{coverages[f] / init1};
      for (const uint64_t c : copies) {
        double total{0.0};
        for (const uint64_t o : counts) {
          if (o + 1 == nObservations) {
            emissions[f][c][o] = log(1 - total);
          } else {
            const double val{gsl_ran_poisson_pdf(
                static_cast<unsigned int>(o), (c + background) * est_rate)};
            total += val;
            emissions[f][c][o] = log(val);
          }
        }
      }
    }
    if (0) show_emissions();
  }

  double viterbi(const uint64_t s) {
    std::vector<StateProbs> T1(nFragments);
    std::vector<StateChoices> T2(nFragments);

    for (uint64_t i{0}; i != nCopies; ++i) {
      T1[0][i] = initial_state[s][i] + emissions[0][i][fragments(s, 0)];
      T2[0][i] = 0;
    }
    for (uint64_t f{1}; f != nFragments; ++f) {
      for (uint64_t j{0}; j != nCopies; ++j) {
        double max{std::numeric_limits<double>::lowest()};
        uint64_t maxk{0};
        for (uint64_t k{0}; k != nCopies; ++k) {
          const double val{T1[f - 1][k] + transitions[f][k][j] +
                emissions[f][j][fragments(s, f)]};
          if (max < val) {
            max = val;
            maxk = k;
          }
        }
        T1[f][j] = max;
        T2[f][j] = maxk;
      }
    }
    unsigned int * const v{viterbis[s]};
    unsigned int z(static_cast<unsigned int>(std::max_element(
        &T1[nFragments - 1][0], &T1[nFragments - 1][0] + nCopies) -
                                             &T1[nFragments - 1][0]));
    viterbi_probs[s] = T1[nFragments - 1][z];
    if (v[nFragments - 1] != z) converged = false;
    v[nFragments - 1] = z;
    for (uint64_t f{nFragments - 1}; f != 0; --f) {
      z = T2[f][v[f]];
      if (v[f - 1] != z) converged = false;
      v[f - 1] = z;
    }
    return viterbi_probs[s];
  }
  double viterbi() {
    ++iteration;
    converged = true;
    static thread_local std::vector<std::future<double> > futures;
    futures.clear();
    for (uint64_t s{0}; s != nSamples; ++s) futures.push_back(pool.run(
             [this](const uint64_t sample) { return viterbi(sample); }, s));
    viterbi_prob = 0;
    for (uint64_t s{0}; s != nSamples; ++s) {
      const double val{futures[s].get()};
      if (0) std::cout << s << " " << val << std::endl;
      viterbi_prob = lnproduct(viterbi_prob, val);
    }
    return viterbi_prob;
  }
  struct update_counts {
    uint64_t lengths[nCopies]{0, 0, 0};
    uint64_t n_state[nCopies]{0, 0, 0};
    uint64_t count{0};
    double expected{0};
    update_counts & operator+=(const update_counts & rhs) {
      for (uint64_t s{0}; s != nCopies; ++s) {
        lengths[s] += rhs.lengths[s];
        n_state[s] += rhs.n_state[s];
      }
      count += rhs.count;
      expected += rhs.expected;
      return *this;
    }
  };
  void recalculate_parameters(const update_counts & nc) {
    background = 1.0 * nc.count / nc.expected;
    fragment_length = 2.0 * nc.lengths[0] / nc.n_state[0];
    init0 = sqrt(1.0 * nc.lengths[0] /
                 (nc.lengths[0] + nc.lengths[1] + nc.lengths[2]));
    init1 = 1 - init0;
    update_model();
  }
  update_counts update_parameters() {
    update_counts uc{};
    static thread_local std::vector<std::future<void> > futures;
    futures.clear();
    std::mutex update_mutex;
    for (uint64_t sample{0}; sample != nSamples; ++sample) {
      futures.push_back(pool.run(
          [this, &uc, &update_mutex]
          (const uint64_t s) {
            update_counts luc;
            uint64_t lastf{0};
            uint64_t last{0};
            segment_starts[s].clear();
            for (uint64_t f{0}; f != nFragments; ++f) {
              const uint64_t state{viterbis[s][f]};
              if (state == 0) {
                luc.count += fragments(s, f);
                luc.expected += coverages[f] / init1;
              }
              if (f + 1 == nFragments || (state != last && f)) {
                segment_starts[s].push_back(static_cast<unsigned int>(lastf));
                const uint64_t length{
                  fragments[f].start() - fragments[lastf].start()};
                luc.lengths[state] += length;
                ++luc.n_state[state];
                lastf = f;
              }
              last = state;
            }
            std::lock_guard<std::mutex> update_lock(update_mutex);
            for (const uint64_t c : copies) {
              uc.lengths[c] += luc.lengths[c];
              uc.n_state[c] += luc.n_state[c];
              uc.count += luc.count;
              uc.expected += luc.expected;
            }
          }, sample));
    }
    for (uint64_t s{0}; s != nSamples; ++s) futures[s].get();
    return uc;
  }

  void show_hmm() const {
    std::cout << "title = " << title << ' '
              << "nSamples = " << nSamples << ' '
              << "nFragments = " << nFragments << '\n';
  }
  void show_all(std::ostream & out) const {
    out << "frag\tpos\tlength\tcover\tsmoothed\tdistance";
    for (const uint64_t i : copies)
      for (const uint64_t o : counts)
        out << "\temit" << i << o;
    for (const uint64_t i : copies)
      for (const uint64_t j : copies)
        out << "\ttrans" << i << j;
    out << '\n';
    for (uint64_t f{0}; f != nFragments; ++f) {
      const uint64_t tot_count{[this, f]() {
          uint64_t result{0};
          for (uint64_t s{0}; s != nSamples; ++s) result += fragments(s, f);
          return result;
        }()};
      out << f
          << '\t' << fragments[f].start()
          << '\t' << fragments[f].length()
          << '\t' << tot_count
          << '\t' << coverages[f] * nSamples * 2
          << '\t' << locus_distances[f];
      for (const uint64_t i : copies)
        for (const uint64_t o : counts)
          out << '\t' << exp(emissions[f][i][o]);
      for (const uint64_t i : copies)
        for (const uint64_t j : copies)
          out << '\t' << exp(transitions[f][i][j]);
      out << '\n';
    }
  }
  void show_iteration() const {
    std::cout
        << title << ",\t"
        << "Iter = " << iteration << ",\t"
        << "init0 = " << init0 << ",\t"
        << "background = " << background << ",\t"
        << "fragment_length = " << fragment_length << ",\t"
        << "per obs prob = " << exp(viterbi_prob / nSamples / nFragments);
    if (converged) std::cout << "\tConverged";
    std::cout << '\n';
  }
  void show_summary() const {
    if (0) show_emissions();
    if (0) show_transitions();
    if (0) show_viterbi();
    if (0) show_path(0);
    if (0) show_sample(std::cout, 0);
    plot_states();
  }
  void show_emissions() const {
    for (uint64_t f{0}; f != nFragments; ++f) {
      std::cout << f << '\t' << coverages[f];
      for (const uint64_t i : copies)
        for (const uint64_t o : counts)
          std::cout << '\t' << exp(emissions[f][i][o]);
      std::cout << '\n';
    }
    std::cout << std::flush;
  }
  void show_transitions() const {
    for (uint64_t f{0}; f != nFragments; ++f) {
      std::cout << f
          // << '\t' << fragments[f].chromosome()
          // << '\t' << fragments[f].start()
                << '\t' << locus_distances[f];
      for (const uint64_t i : copies)
        for (const uint64_t j : copies)
          std::cout << '\t' << exp(transitions[f][i][j]);
      std::cout << '\n';
    }
    std::cout << std::flush;
  }
  void show_viterbi() const {
    std::vector<unsigned int> state_counts(nCopies);
    // std::map<unsigned int, unsigned int> state_counts;
    uint64_t total_copy{0};
    for (uint64_t s{0}; s != nSamples; ++s) {
      for (uint64_t f{0}; f != nFragments; ++f) {
        const uint64_t copy{viterbis[s][f]};
        ++state_counts[copy];
        total_copy += copy;
      }
    }
    std::cout << "Average copy = " << 1.0 * total_copy / nSamples / nFragments
              << "\n";
    for (uint64_t s{0}; s != nCopies; ++s) {
      std::cout << "State " << s
                << " seen " << state_counts[s] << " times" << std::endl;
    }
    if (1) show_path(0);
  }
  void show_path(const uint64_t sample) const {
    uint64_t lastf{0};
    uint64_t last{0};
    uint64_t total{0};
    std::cout << "state\tnobs\tfrags\tlength\n";
    for (uint64_t f{0}; f != nFragments; ++f) {
      const uint64_t state{viterbis[sample][f]};
      if (f + 1 == nFragments || (state != last && f)) {
        std::cout
            << last
            << '\t' << total
            << '\t' << f - lastf
            << '\t' << fragments[f].start() - fragments[lastf].start()
            << '\n';
        total = 0;
        lastf = f;
      }
      total += fragments(sample, f);
      last = state;
    }
  }

  void show_sample(std::ostream & stream, const uint64_t sample) const {
    for (uint64_t f{0}; f != nFragments; ++f) {
      stream << viterbis[sample][f]
             << '\t' << counts[fragments(sample, f)]
             << '\t' << fragments[f].chromosome()
             << '\t' << fragments[f].start()
             << '\t' << locus_distances[f]
             << "\n";
    }
    stream << std::flush;
  }
  void plot_states() const {
    std::cout << "Making plots for coverage HMM on " << title << std::endl;
    PSGraph * const states_graph{new PSGraph{plots,
            title + " states " +
            std::to_string(iteration) + ";Fragment Index;Sample Index",
            Bounds{0.0, 1.0 * nFragments, 0.0, 1.0 * nSamples}}};
    ownp(states_graph);
    const Marker marker{paa::circle(), 0.5, "0 0 0", 0.5, true};
    PSXYSeries * const ploidy_series{new PSXYSeries{
        plots, title + " Viterbi ploidy for iteration " +
            std::to_string(iteration) + ";Sample Rank;Ploidy", marker}};
    ploidy_series->parents()[0]->own(
        std::unique_ptr<PSXYSeries>(ploidy_series));
    PSXYSeries * const copy_series{new PSXYSeries{
        plots, title + " Viterbi average copy number for iteration " +
            std::to_string(iteration) + ";Position;Average Ploidy", marker}};
    ownp(copy_series);
    PSXYSeries * const prob_series{new PSXYSeries{
        plots, title + " Viterbi observation probabilities for iteration " +
            std::to_string(iteration) +
            ";Sample Rank;Observation probability",
            marker}};
    ownp(prob_series);
    states_graph->ps(R"xxx(0 setlinewidth
/box {
newpath /y exch def /x2 exch def /x1 exch def
x1 y gc moveto setrgbcolor
x2 y gc lineto
x2 y 0.5 add gc lineto
x1 y 0.5 add gc lineto
closepath
fill 
} def
)xxx");
    const char * colors[]{"1 1 1", "1 0 0", "0 0 1"};
    thread_local std::vector<double> ploidies;
    ploidies.clear();
    thread_local std::vector<double> copy_no;
    copy_no.assign(nFragments, 0);
    thread_local std::vector<double> probs;
    probs.clear();
    for (uint64_t sample{0}; sample != nSamples; ++sample) {
      uint64_t ploidy{0};
      const std::vector<unsigned int> & starts{segment_starts[sample]};
      std::ostringstream ps;
      for (uint64_t segment{0}; segment != starts.size(); ++segment) {
        const uint64_t start{starts[segment]};
        const uint64_t stop{segment + 1 == starts.size() ?
              nFragments : starts[segment + 1]};
        const uint64_t copy{viterbis[sample][start]};
        ploidy += (stop - start) * copy;
        for (uint64_t f{start}; f != stop; ++f) {
          copy_no[f] += copy;
        }
        if (copy) {
          ps << colors[copy]
             << " " << start << " " << stop << " " << sample << " box\n";
        }
      }
      ploidies.push_back(ploidy);
      probs.push_back(viterbi_probs[sample]);
      states_graph->ps(ps.str());
    }
    for (double & pos : copy_no) pos /= nSamples;
    sort(ploidies.begin(), ploidies.end());
    sort(probs.begin(), probs.end());
    for (uint64_t p{0}; p != ploidies.size(); ++p) {
      const double ploidy{ploidies[p]};
      ploidy_series->add_point(p, 1.0 * ploidy / nFragments);
    }
    for (uint64_t p{0}; p != copy_no.size(); ++p) {
      copy_series->add_point(fragments[p].start(), copy_no[p]);
    }
    for (uint64_t p{0}; p != probs.size(); ++p) {
      const double prob{probs[p]};
      if (0) std::cout << p << " " << prob << " " << exp(prob / nFragments)
                       << std::endl;
      prob_series->add_point(p, exp(prob / nFragments));
    }
  }

  const std::string title;
  const FragmentTable & fragments;
  const uint64_t nSamples;
  const uint64_t nFragments;
  ThreadPool & pool;
  const std::vector<unsigned int> locus_distances;

  double init0;
  double init1;
  double background;
  double fragment_length;

  PSDoc & plots;

  std::vector<double> coverages;

  std::vector<Trans> transitions;
  std::vector<Emit> emissions;
  std::vector<StateProbs> initial_state{};

 private:
  Matrix<unsigned int> viterbis;
  std::vector<double> viterbi_probs;
  double viterbi_prob{SMALL};
  std::vector<std::vector<unsigned int> > segment_starts;
  uint64_t iteration{0};
  static std::mutex cout_mutex;

 public:
  bool converged{false};
};

std::mutex CoverageHMM::cout_mutex{};

inline unsigned char encode_base(const char base) {
  switch (base) {
    case 'A':
      return 0;
    case 'C':
      return 1;
    case 'G':
      return 2;
    case 'T':
      return 3;
    default:
      throw Error("Unknown base in encode_base:") << base;
  }
}
constexpr char bases[]{"ACGT"};

class HetInfo {
 public:
  HetInfo(const unsigned int chromosome__,
          const unsigned int position__,
          const unsigned char phase1__,
          const unsigned char phase2__,
          const bool phase_known = true) :
      position_{position__},
    chromosome_{static_cast<unsigned char>(chromosome__)},
    phase1_{phase_known ? phase1__ :
          (phase1__ < phase2__ ? phase1__ : phase2__)},
    phase2_{phase_known ? phase2__ :
          (phase1__ < phase2__ ? phase2__ : phase1__)} { }

  unsigned int chromosome() const { return chromosome_; }
  unsigned int position() const { return position_; }
  char phase1() const { return bases[phase1_]; }
  char phase2() const { return bases[phase2_]; }

 private:
  unsigned int position_;
  unsigned char chromosome_;
  unsigned char phase1_;
  unsigned char phase2_;
};

using uArray4 = Array<unsigned int, 4>;
constexpr const char * haplo_symbols[4]{"0", "M", "F", "MF"};
constexpr const char * alleles_symbols[4]{"_", "A", "B", "AB"};
constexpr uArray4 haplos{{0, 1, 2, 3}};
constexpr uArray4 flip_haplo{{0, 2, 1, 3}};
constexpr uArray4 haplo_copy{{0, 1, 1, 2}};
constexpr uArray4 het_obs{{0, 1, 2, 3}};
constexpr uArray4 flip_obs{{0, 2, 1, 3}};

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
      std::cout << "Loaded " << n_hets() << " hets from "
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
      std::cout << "Constructing het info binary cache for "
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
      unsigned int id;
      char phase1;
      char phase2;
      std::string state_string;
      unsigned int n_rows{0};
      unsigned int last_n_cols{0};
      while (getline(in_stream, line)) {
        ++n_rows;
        std::istringstream line_stream{line.c_str()};
        line_stream >> chr_name >> pos >> id >> phase1 >> phase2;
        if (!line_stream) throw Error("Het info parse error in HetTable");
        const HetInfo het_info{lookup[chr_name], pos,
              encode_base(phase1), encode_base(phase2)};
        temp_hets.push_back(het_info);
        unsigned int n_cols{0};
        obs_counts.push_back(0);
        A_counts.push_back(0);
        B_counts.push_back(0);
        while (line_stream >> state_string) {
          const unsigned int state{[&state_string]() {
              if (state_string == "_") {
                return 0U;
              } else if (state_string == "A") {
                return 1U;
              } else if (state_string == "B") {
                return 2U;
              } else if (state_string == "AB") {
                return 3U;
              } else {
                throw Error("Bad state string in HetTable:") << state_string;
              }
            }()};
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

class LowessCoverageLookup {
 public:
  using Lookup = std::pair<unsigned int, double>;
  LowessCoverageLookup(const FragmentTable & fragment_table,
                       const std::vector<double> & smoothed_coverage) :
      direct(smoothed_coverage.size()),
      lookup(smoothed_coverage.size()) {
    if (smoothed_coverage.size() != fragment_table.n_fragments())
      throw Error("LCL error");
    for (uint64_t f{0}; f != fragment_table.n_fragments(); ++f) {
      const double coverage{smoothed_coverage[f]};
      direct[f] = coverage;
      lookup[f] = {fragment_table[f].length(), coverage};
    }
    sort(lookup.begin(), lookup.end(), sort_by_length);
  }
  double operator[](const unsigned int fragment) const {
    return direct[fragment];
  }
  double operator()(const unsigned int fragment_length) const {
    auto found = equal_range(lookup.begin(), lookup.end(),
                             Lookup{fragment_length, 0.0}, sort_by_length);
    // Note no interpolation is done
    if (found.first == lookup.end()) return lookup.back().second;
    return (found.first + (found.second - found.first) / 2)->second;
  }
  uint64_t size() const { return direct.size(); }
  static bool sort_by_length(const Lookup & lhs, const Lookup & rhs) {
    return lhs.first < rhs.first;
  }

 private:
  std::vector<double> direct;
  std::vector<Lookup> lookup;
};


std::vector<double> simple_push(const HetTable & hets,
                                const bool set_to_truth = false) {
  std::vector<double> p_mom_as(hets.n_hets(), log(0.5));
  // special code to set to truth
  if (set_to_truth) {
    const double push{0.999};
    for (uint64_t h{0}; h != hets.n_hets(); ++h)
      p_mom_as[h] = hets.flip(h) ? log(1 - push) : log(push);
  } else {
    std::vector<uint64_t> well_covered_loci;
    const unsigned int min_singles{5};
    for (uint64_t h{0}; h != hets.n_hets(); ++h) {
      unsigned int n_singles{0};
      for (uint64_t s{0}; s != hets.n_samples(); ++s) {
        const unsigned char obs{hets(s, h)};
        if (obs == HetTable::ObsA || obs == HetTable::ObsB) ++n_singles;
      }
      if (n_singles >= min_singles) well_covered_loci.push_back(h);
    }
    if (well_covered_loci.empty())
      throw Error("Could not find well covered locus for push");
    std::random_device rd{};
    std::mt19937_64 mersenne{rd()};
    std::function<uint64_t()> gen{std::bind(
        std::uniform_int_distribution<uint64_t>(
            0, well_covered_loci.size() - 1), mersenne)};
    for (unsigned int i{0}; i != 1; ++i) {
      const uint64_t het{well_covered_loci[gen()]};
      // const double push{1 - 2 * std::numeric_limits<double>::epsilon()};
      const double push{0.999};
      p_mom_as[het] = hets.flip(het) ? log(1 - push) : log(push);
    }
  }
  return p_mom_as;
}

const unsigned int digits_precision{14};
const double min_confidence{1 / pow(10, digits_precision)};
struct Call {
  Call(const unsigned int id_,
       const unsigned int chr_, const unsigned int pos_,
       const double p_mom_a_) :
      id{id_}, chr{chr_}, pos{pos_}, p_mom_a{p_mom_a_} {
        if (p_mom_a <= 0.5 && p_mom_a >= 0.5)
          throw Error("No 0.5 p_mom_a allowed");
      }
  unsigned int id : 25;
  unsigned int chr :7;
  unsigned int pos;
  double p_mom_a;
  bool mom_a() const { return p_mom_a > 0.5; }
  bool agree(const Call & other) const { return mom_a() == other.mom_a(); }
  double confidence() const {
    return 0.5 - fabs(p_mom_a - 0.5) + min_confidence;
  }
  double confidence_score() const { return -log10(confidence()); }
};

const std::vector<std::string> colors{
  "0.898039 0 0", "0.145098 0 0.619608",
      "0 0.717647 0", "0.898039 0.745098 0",
      "0.0235294 0.337255 0.576471", "0.717647 0.866667 0",
      "0.898039 0.513725 0", "0.584314 0 0.584314",
      "0.972549 0.470588 0.972549", "0 0.0941176 0",
      "0 0.941176 0.533333", "0.564706 0.627451 0.533333",
      "0.972549 0.972549 0.627451", "0 0.658824 0.972549",
      "0.439216 0.313725 0.972549", "0.972549 0.0313725 0.972549",
      "0.470588 0.282353 0.188235", "0.972549 0.25098 0.470588",
      "0.470588 0.972549 0.376471", "0 0.156863 0.972549",
      "0.439216 0.596078 0", "0.12549 0.627451 0.376471",
      "0.972549 0.596078 0.470588", "0.627451 0.658824 0.972549",
      "0.282353 0.972549 0", "0.0627451 0.407843 0.0941176",
      "0.439216 0 0", "0.0313725 0.972549 0.909804",
      "0.376471 0 0.941176", "0.596078 0.313725 0.596078",
      "0.972549 0.784314 0.972549", "0.282353 0.784314 0.721569",
      "0.12549 0.156863 0.313725", "0.972549 0.972549 0.219608",
      "0.282353 0.501961 0.752941", "0.658824 0.878431 0.690196",
      "0.815686 0.25098 0.12549", "0.784314 0.25098 0.909804",
      "0 0.972549 0.219608", "0.878431 0 0.376471",
      "0.690196 0.470588 0.282353", "0.784314 0.784314 0.345098",
      "0 0.407843 0.909804", "0.345098 0.439216 0.439216",
      "0.282353 0.219608 0.721569", "0 0.627451 0.658824",
      "0.407843 0 0.313725", "0.376471 0.752941 0.25098",
      "0.784314 0.533333 0.721569", "0.784314 0.972549 0.972549",
      "0.690196 0 0.909804", "0.658824 0.156863 0.376471",
      "0.627451 0.407843 0", "0.313725 0.658824 0.972549",
      "0.847059 0.0941176 0.690196", "0.25098 0.25098 0",
      "0.878431 0.752941 0.658824", "0.376471 0.219608 0.439216",
      "0.972549 0.407843 0.25098", "0.407843 0.972549 0.658824",
      "0.658824 0.439216 0.941176", "0.690196 0.658824 0.12549",
      "0.658824 0 0.156863", "0.533333 0.156863 0.815686",
      "0.0627451 0.407843 0.345098", "0.25098 0.847059 0.470588",
      "0.188235 0.596078 0.12549", "0.188235 0 0.156863",
      "0.721569 0.972549 0.25098", "0 0.156863 0.721569",
      "0.972549 0.407843 0.658824", "0.470588 0.815686 0",
      "0 0.815686 0.752941", "0.596078 0.188235 0",
      "0.470588 0.564706 0.25098", "0.0941176 0.784314 0.219608",
      "0 0 0.407843", "0.313725 0.627451 0.564706",
      "0.533333 0.533333 0.752941", "0.0941176 0 0.847059",
      "0.25098 0.847059 0.972549", "0.0313725 0.909804 0",
      "0.784314 0.407843 0.501961", "0.25098 0.439216 0.972549",
      "0.909804 0.627451 0.219608", "0.470588 0.784314 0.878431",
      "0.313725 0.439216 0.0627451", "0.564706 0.815686 0.470588",
      "0.972549 0.12549 0.188235", "0.752941 0.972549 0.501961",
      "0.25098 0.941176 0.25098", "0.501961 0.972549 0.12549",
      "0.972549 0.627451 0.815686", "0.376471 0.0313725 0.690196",
      "0.25098 0.188235 0.941176", "0 0.752941 0.501961",
      "0.156863 0.972549 0.721569", "0.596078 0.815686 0.219608",
      "0.972549 0.847059 0.439216", "0.470588 0.972549 0.972549"
      };

class HaHaHMM {
 public:
  static constexpr double SMALL{std::numeric_limits<double>::lowest()};
  static constexpr uint64_t nHaplos{haplos.size()};
  static_assert(nHaplos == 4, "nHaplos check");

  static constexpr uint64_t nObservations{het_obs.size()};
  static_assert(nObservations == 4, "nObservations check");

  using Trans = Matrix<double, nHaplos, nHaplos>;
  using Emit = Matrix<double, nHaplos, nObservations>;
  using HaploProbs = Array<double, nHaplos>;

  struct LocusInfo {
    LocusInfo(const unsigned int index_,
              const bool is_het_) :
        index{index_}, is_het{is_het_} { }
    unsigned int index;
    bool is_het;
  };

  HaHaHMM(const std::string & title_,
          ThreadPool & pool_,
          PSDoc & plots_,
          const HetTable & hets_,
          const CoverageHMM & coverage_hmm_,
          const std::vector<double> initial_phase_ = std::vector<double>()) :
      title{title_}, pool{pool_}, plots{plots_},
    hets{hets_}, coverage_hmm{coverage_hmm_},
    p_mom_as{initial_phase_.size() ? initial_phase_ : simple_push(hets)} {
      initialize();
      make_initial_plots();
      run();
    }

  void initialize() {
    auto fut = pool.run([this]() { calculate_transitions(); });
    calculate_mom_a_emissions();
    fut.get();
    update_model();
    show_hmm();
  }

  void make_initial_plots() {
    using Hist = PSHSeries<unsigned int, unsigned int>;
    Hist * const hets_hist{new Hist{plots,
            "Het distribution;Position;N",
            Bounds{0, 1.0 * hets.max_position()}}};
    ownp(hets_hist);
    for (uint64_t h{0}; h != nHets; ++h)
      hets_hist->add_point(hets[h].position());
  }
  void update_model() {
    initial_haplo.assign(calculate_initial_haplo());
    calculate_het_sample_emissions();
  }
  HaploProbs calculate_initial_haplo() const {
    return {{log(init0 * init0), log(init0 * init1),
            log(init0 * init1), log(init1 * init1)}};
  }
  void calculate_transitions() {
    const double rate10{1 / fragment_length};
    const double rate01{rate10 * init0 / (1 - init0)};
    const Matrix<double, 2, 2> strandT{
      {1 - rate01, rate01}, {rate10, 1 - rate10}};
    for (uint64_t d{1}; d != nLoci; ++d) {
      const Matrix<double, 2, 2> strandP{
        matrix_power(strandT, locus_distances[d])};
      transitions[d] = {
        {log(strandP[0][0] * strandP[0][0]),
         log(strandP[0][1] * strandP[0][0]),
         log(strandP[0][0] * strandP[0][1]),
         log(strandP[0][1] * strandP[0][1])},
        {log(strandP[1][0] * strandP[0][0]),
         log(strandP[1][1] * strandP[0][0]),
         log(strandP[1][0] * strandP[0][1]),
         log(strandP[1][1] * strandP[0][1])},
        {log(strandP[0][0] * strandP[1][0]),
         log(strandP[0][1] * strandP[1][0]),
         log(strandP[0][0] * strandP[1][1]),
         log(strandP[0][1] * strandP[1][1])},
        {log(strandP[1][0] * strandP[1][0]),
         log(strandP[1][1] * strandP[1][0]),
         log(strandP[1][0] * strandP[1][1]),
         log(strandP[1][1] * strandP[1][1])}};
    }
  }
  void show_transitions(std::ostream & out) const {
    for (uint64_t l{0}; l != nLoci; ++l) {
      out << l << '\t' << locus_distances[l];
      for (uint64_t i{0}; i != nHaplos; ++i) {
        for (uint64_t j{0}; j != nHaplos; ++j) {
          out << '\t' << exp(transitions[l][i][j]);
        }
      }
      out << '\n';
    }
  }

  // loci x haplo x obs assuming mom is A allele
  void calculate_mom_a_emissions() {
    // Assumes a locus model of M = A!
    for (uint64_t h{0}; h != nHets; ++h) {
      const double c{coverages[h] / init1};
      const double cq{1 - c};
      const double cb{c * background};
      const double cbq{1 - cb};
      mom_a_emissions[h] = {
        {log(cbq * cbq), log(cb * cbq), log(cbq * cb), log(cb * cb)},
        {log(cq  * cbq), log(c  * cbq), log(cq * cb),  log(c  * cb)},
        {log(cbq * cq),  log(cb * cq),  log(cbq * c),  log(cb * c)},
        {log(cq  * cq),  log(c  * cq),  log(cq  * c),  log(c  * c)}};
    }
  }
  // Return emission prob for obs in haplo at het, based on mom_is_a
  double get_het_emission(const uint64_t het,
                          const uint64_t haplo,
                          const uint64_t obs,
                          const bool mom_is_a) const {
    return mom_is_a ? mom_a_emissions[het][haplo][obs] :
        mom_a_emissions[het][flip_haplo[haplo]][obs];
  }
  // Return emission prob for obs in haplo at het, based on p_mom_a
  double get_het_emission(const uint64_t het,
                          const uint64_t haplo,
                          const uint64_t obs,
                          const double log_p_mom_a,
                          const double log_q_mom_a) const {
    const uint64_t flipped_haplo{flip_haplo[haplo]};
    const bool symmetric_obs{obs == HetTable::Obs0 || obs == HetTable::ObsAB};
    return (flipped_haplo != haplo && !symmetric_obs) ?
        lnsum(log_p_mom_a + mom_a_emissions[het][haplo][obs],
              log_q_mom_a + mom_a_emissions[het][flipped_haplo][obs]) :
        mom_a_emissions[het][haplo][obs];
  }
  // locus * haplo * obs standard emission matrices for each het locus
  void calculate_het_emissions() {
    for (uint64_t h{0}; h != nHets; ++h) {
      // parallelize here?
      const double log_p_mom_a{p_mom_as[h]};
      const double log_q_mom_a{log(1 - exp(log_p_mom_a))};
      if (log_p_mom_a > 0 || log_q_mom_a > 0) throw Error("prob prob");
      for (uint64_t i{0}; i != nHaplos; ++i)
        for (uint64_t o{0}; o != nObservations; ++o)
          het_emissions[h][i][o] =
              get_het_emission(h, i, o, log_p_mom_a, log_q_mom_a);
    }
  }
  // sample * locus * haplo pre-calculated emissions for obs
  void calculate_het_sample_emissions() {
    calculate_het_emissions();
    static thread_local std::vector<std::future<void> > futures;
    futures.clear();
    for (uint64_t s{0}; s != nSamples; ++s)
      futures.push_back(pool.run([this](const uint64_t sample) {
            for (uint64_t h{0}; h != nHets; ++h)
              for (uint64_t i{0}; i != nHaplos; ++i)
                het_sample_emissions[sample][h][i] =
                    het_emissions[h][i][hets(sample, h)];
          }, s));
    for (uint64_t s{0}; s != nSamples; ++s) futures[s].get();
  }
  // Returns emission probabilities over fragments and hets
  double get_emission(const uint64_t sample,
                      const uint64_t locus,
                      const uint64_t haplo) const {
    const LocusInfo info{loci[locus]};
    if (info.is_het) {
      return het_sample_emissions[sample][info.index][haplo];
    } else {
      const uint64_t obs{fragments(sample, info.index)};
      return coverage_hmm.emissions[info.index][haplo_copy[haplo]][obs];
    }
  }
  void show_emissions(std::ostream & out) const {
    for (uint64_t l{0}; l != nLoci; ++l) {
      const LocusInfo info{loci[l]};
      out << l << "\t" << info.is_het;
      for (uint64_t i{0}; i != nHaplos; ++i) {
        for (uint64_t o{0}; o != nObservations; ++o) {
          out << '\t';
          if (info.is_het) {
            out << exp(het_emissions[info.index][i][o]);
          } else {
            out << exp(coverage_hmm.emissions
                       [info.index][haplo_copy[i]][haplo_copy[o]]);
          }
        }
      }
      out << '\n';
    }
  }

  void run() {
    iteration = 0;
    converged = false;
    while (!converged && iteration != max_iterations) {
      const double p{forward_backward()};
      if (p_mom_a_diffs < p_mom_a_tol) converged = true;
      show_iteration(p);
      if (false) {
        if (converged) {
          plot_iteration(25, true, 1);
        } else {
          plot_iteration();
        }
      }
    }
    std::ofstream phase_out{title + ".phase.txt"};
    show_loci(phase_out);
    if (converged) std::cout << "Converged!" << std::endl;
  }
  double forward(const uint64_t s) {
    HaploProbs * const f{forwards[s]};
    for (uint64_t i{0}; i != nHaplos; ++i)
      f[0][i] = initial_haplo[s][i] + get_emission(s, 0, i);
    for (uint64_t l{1}; l != nLoci; ++l) {
      for (uint64_t j{0}; j != nHaplos; ++j) {
        double logalpha{0};
        for (uint64_t i{0}; i != nHaplos; ++i) {
          const double product{
            lnproduct(f[l - 1][i], transitions[l - 1][i][j])};
          logalpha = i ? lnsum(logalpha, product) : product;
        }
        f[l][j] = lnproduct(logalpha, get_emission(s, l, j));
      }
    }
    return lnsum(f[nLoci - 1], nHaplos);
  }
  double backward(const uint64_t s) {
    HaploProbs * const b{backwards[s]};
    for (uint64_t i{0}; i != nHaplos; ++i) b[nLoci - 1][i] = 0;
    for (uint64_t ll{1}; ll != nLoci; ++ll) {
      const uint64_t l{nLoci - ll - 1};
      for (uint64_t i{0}; i != nHaplos; ++i) {
        double logbeta{0};
        for (uint64_t j{0}; j != nHaplos; ++j) {
          const double product{lnproduct(
              transitions[l][i][j], lnproduct(get_emission(s, l + 1, j),
                                              b[l + 1][j]))};
          logbeta = j ? lnsum(logbeta, product) : product;
        }
        b[l][i] = logbeta;
      }
    }
    double result{0};
    for (uint64_t i{0}; i != nHaplos; ++i) {
      const double product{
        b[0][i] + get_emission(s, 0, i) + initial_haplo[s][i]};
      result = i ? lnsum(product, result) : product;
    }
    return result;
  }
  void forward_backward(const uint64_t s) {
    const double fp{forward(s)};
    const double bp{backward(s)};
    if (0) std::cout << "f b " << fp << " " << bp
                     << " " << exp(bp) / exp(fp) << std::endl;
    sample_probs[s] = fp;
    for (uint64_t l{0}; l != nLoci; ++l) {
      for (uint64_t i{0}; i != nHaplos; ++i) {
        posterior[s][l][i] = lnproduct(forwards[s][l][i], backwards[s][l][i],
                                       -fp);
        const double fval{lnproduct(forwards[s][l][i],
                                    backwards[s][l][flip_haplo[i]])};
        const LocusInfo & locus{loci[l]};
        if (locus.is_het) {
          const uint64_t h{locus.index};
          p_flip[s][h] = i ? lnsum(fval, p_flip[s][h]) : fval;
        }
      }
    }
  }
  double forward_backward() {
    // Run forward-backward on all samples in parallel
    static thread_local std::vector<std::future<void> > futures;
    futures.clear();
    for (uint64_t s{0}; s != nSamples; ++s)
      futures.push_back(pool.run([this](const uint64_t sample) {
            return forward_backward(sample); }, s));
    for (uint64_t s{0}; s != nSamples; ++s) futures[s].get();

    // Get overall probability
    const double p_all{[this]() {
        double result{0};
        for (const double prob : sample_probs)
          result = lnproduct(prob, result);
        return result;
      }()};

    flip_diff = 0;
    flip_index = 0;
    if (++iteration >= first_flip_iteration) {
      // Determine if a flip is needed and do flip if so
      for (uint64_t h{0}; h != nHets; ++h) {
        flip_diffs[h] = -p_all;
        for (uint64_t s{0}; s != nSamples; ++s)
          flip_diffs[h] = lnproduct(flip_diffs[h], p_flip[s][h]);
      }
      flip_index = static_cast<uint64_t>(
          max_element(flip_diffs.begin(), flip_diffs.end() - 1) -
          flip_diffs.begin());
      if ((flip_diff = flip_diffs[flip_index]) >= flip_threshold) {
        ++n_flips;
        for (uint64_t h{flip_index + 1}; h != nHets; ++h)
          p_mom_as[h] = log(1 - exp(p_mom_as[h]));
        update_model();
        return p_all;
      }
    }
    calculate_p_mom_as();
    // if (iteration % 20 == 0) p_mom_a_sweep();
    update_model();
    return p_all;
  }

  void show_iteration(const double p) {
    const int64_t old_precision{std::cout.precision()};
    std::cout.precision(10);
    std::cout << "Iteration " << iteration
              << " prob " << exp(p / nLoci / nSamples)
              << " flip " << n_flips << " " << flip_index
              << " " << exp(flip_diff)
              << " diff " << p_mom_a_diffs
              << std::endl;
    std::cout.precision(old_precision);
    return;
    const std::string iter_id{std::to_string(iteration) +
          "." + std::to_string(n_flips)};
    std::ofstream emissions_stream{
      (std::string("emit.") + iter_id + ".txt").c_str()};
    show_emissions(emissions_stream);
    std::ofstream transitions_stream{
      (std::string("trans.") + iter_id + ".txt").c_str()};
    show_transitions(transitions_stream);
    std::ofstream p_mom_as_stream{
      (std::string("p_mom_a.") + iter_id + ".txt").c_str()};
    show_p_mom_as(p_mom_as_stream);
  }

  void show_loci(std::ostream & out) const {
    const int64_t old_precision{out.precision()};
    out.precision(20);
    for (uint64_t h{0}; h != nHets; ++h) {
      const double prob{hets.flip(h) ? 1 - exp(p_mom_as[h]) : exp(p_mom_as[h])};
      out << h << " " << hets[h].position()
          << " " << hets[h].phase1()
          << " " << hets[h].phase2()
          << " " << prob
          << " " << 0.5 - fabs(0.5 - prob)
          << " " << flip_diffs[h];
      if (false) {
        for (uint64_t s{0}; s != nSamples; ++s)
          out << std::max_element(&posterior[s][h][0],
                                  &posterior[s][h][0] + nHaplos) -
              &posterior[s][h][0];
      }
      out << '\n';
    }
    out.precision(old_precision);
  }

  void show_hmm() const {
    std::cout << "title = " << title << ' '
              << "nSamples = " << nSamples << ' '
              << "nHets = " << nHets << ' '
              << "nFragments = " << nFragments << ' '
              << "nLoci = " << nLoci << '\n';
  }

  void calculate_p_mom_as() {
    const double min_prob{1.0 / nSamples};
    const double log_min_p{log(min_prob)};
    const double log_max_p{log(1 - min_prob)};
    const bool limit_ps{false};
    p_mom_a_diffs = 0;
    for (uint64_t l{0}; l != nLoci; ++l) {
      if (loci[l].is_het) {
        const uint64_t h{loci[l].index};
        double model_probs[2]{0.0, 0.0};
        for (uint64_t s{0}; s != nSamples; ++s) {
          for (const bool mom_is_a : {false, true}) {
            double sample_prob{0.0};
            for (uint64_t i{0}; i != nHaplos; ++i) {
              const double obs_prob{posterior[s][l][i] +
                    get_het_emission(h, i, hets(s, h), mom_is_a)};
              sample_prob = i ? lnsum(sample_prob, obs_prob) : obs_prob;
            }
            model_probs[mom_is_a] += sample_prob;
          }
        }
        const double last_p_mom_a{exp(p_mom_as[h])};
        p_mom_as[h] = model_probs[1] - lnsum(model_probs, 2);
        const double p_mom_a_diff{last_p_mom_a - exp(p_mom_as[h])};
        p_mom_a_diffs += fabs(p_mom_a_diff);
        if (limit_ps) {
          if (p_mom_as[h] < log_min_p) p_mom_as[h] = log_min_p;
          if (p_mom_as[h] > log_max_p) p_mom_as[h] = log_max_p;
        }
      }
    }
  }

  void p_mom_a_sweep() {
    // Does not work as hoped....
    std::cerr << "do sweep..." << std::flush;
    const uint64_t block_size{500};
    bool last_flip{false};
    for (uint64_t r{1}; r != nHets; ++r) {
      if (!(p_mom_as[r] < log(0.5) || p_mom_as[r] > log(0.5))) continue;
      const bool r_mom_a{p_mom_as[r] > log(0.5)};
      int64_t val{0};
      const uint64_t start{block_size > r ? 0 : r - block_size};
      for (uint64_t s{0}; s != nSamples; ++s) {
        const unsigned int r_obs{hets(s, r)};
        if (r_obs == 0 || r_obs == 3) continue;
        for (uint64_t l{start}; l != r; ++l) {
          if (!(p_mom_as[l] < log(0.5) || p_mom_as[l] > log(0.5))) continue;
          const unsigned int l_obs{hets(s, l)};
          if (l_obs == 0 || l_obs == 3) continue;
          const bool l_mom_a{p_mom_as[l] > log(0.5)};
          if ((l_mom_a == r_mom_a && l_obs == r_obs) ||
              (l_mom_a != r_mom_a && l_obs != r_obs)) {
            ++val;
          } else {
            --val;
          }
        }
      }
      const bool do_flip{val != 0 ? val < 0 : last_flip};
      if (do_flip) p_mom_as[r] = log(1 - exp(p_mom_as[r]));
      last_flip = do_flip;
    }
    std::cerr << " done" << std::endl;
  }

  void show_p_mom_as(std::ostream & out) {
    for (uint64_t h{0}; h != nHets; ++h)
      out << h << " " << exp(p_mom_as[h]) << "\n";
  }

  void plot_iteration(const uint64_t n_per_page = 5,
                      const bool show_all = false,
                      const unsigned int sample = 10) {
    const Marker flip_marker{paa::circle(), 0.1, "0 0 1", 0.1, true, "0 0 1"};
    const std::vector<std::string> state_colors{
      "0 0 0", "1 0 0", "0 1 0", "0 0 1"};
    const std::vector<Marker> state_markers{[&state_colors]() {
        std::vector<Marker> result;
        for (const std::string color : state_colors)
          result.emplace_back(paa::circle(), 0.1, color, 0.1, true, color);
        return result;
      }()};

    doc_defaults.ticks(false);
    for (uint64_t start_sample{0}; start_sample < nSamples;
         start_sample += n_per_page) {
      const uint64_t stop_sample{
        std::min(start_sample + n_per_page, nSamples)};
      const uint64_t n_this_page{stop_sample - start_sample};
      std::ostringstream layout_stream;
      layout_stream << "3 (0.8 0.1) 1 ^0 0 " << n_this_page << " 1^";
      std::ostringstream title_stream;
      title_stream << title
                   << ", samples " << start_sample + 1 << " - " << stop_sample
                   << ", iteration " << iteration;
      PSPage * const phase_page{
        new PSPage{plots, title_stream.str(), layout_stream.str()}};
      for (uint64_t s{start_sample}; s != stop_sample; ++s) {
        PSGraph * const state_graph{new PSGraph{*phase_page, ""}};
        state_graph->do_border(false);
        for (const unsigned int i : haplos) {
          PSXYSeries * state_series{new PSXYSeries{
              *state_graph, state_markers[i]}};
          for (unsigned int l{0}; l < nLoci; l += sample)
            state_series->add_point(exp(posterior[s][l][i]),
                                    actual_locus_position(l));
          ownp(state_series);
        }
        ownp(state_graph);
      }
      PSGraph * const model_graph{new PSGraph{
          *phase_page, "", Bounds{0, 1, 0.0, 1.0 * max_pos()}}};
      PSXYSeries * model_series{new PSXYSeries{*model_graph, flip_marker}};
      PSGraph * const flip_graph{new PSGraph{
          *phase_page, "", Bounds{-200.0, 50.0, 0.0, 1.0 * max_pos()}}};
      flip_graph->ps("1 0 0 c 1 lw np " + std::to_string(flip_threshold) +
                     " 0 gc moveto " + std::to_string(flip_threshold) + " " +
                     std::to_string(max_pos()) + " gc lineto stroke");
      PSXYSeries * flip_series{new PSXYSeries{*flip_graph, flip_marker}};
      for (unsigned int h{0}; h != nHets; ++h) {
        const double p_mom_a{hets.flip(h) ? 1 - exp(p_mom_as[h]) :
              exp(p_mom_as[h])};
        model_series->add_point(p_mom_a, hets[h].position());
        flip_series->add_point(flip_diffs[h], hets[h].position());
      }
      ownp(flip_series);
      ownp(flip_graph);
      ownp(model_series);
      ownp(model_graph);
      ownp(phase_page);
      if (!show_all) break;
    }
  }

  unsigned int max_pos() const {
    return hets[nHets - 1].position();
  }

 private:
  // removed fragments leads to fragment misidentification for some hets!
  const std::string title;
  ThreadPool & pool;
  PSDoc & plots;
  const HetTable & hets;
  const CoverageHMM & coverage_hmm;

  const double init0{coverage_hmm.init0};
  const double init1{1 - init0};
  const double background{coverage_hmm.background};
  const double fragment_length{coverage_hmm.fragment_length};

  const FragmentTable & fragments{coverage_hmm.fragments};
  const uint64_t nFragments{coverage_hmm.nFragments};
  const uint64_t nSamples{hets.n_samples()};
  const uint64_t nHets{hets.n_hets()};
  const std::vector<LocusInfo> loci{get_loci()};
  const uint64_t nLoci{loci.size()};
  const std::vector<unsigned int> locus_distances{get_distances()};
  const std::vector<double> coverages{get_coverages()};

  Matrix<double, 0, nHaplos> initial_haplo{nSamples};
  std::vector<Trans> transitions{std::vector<Emit>(nLoci)};
  std::vector<Emit> mom_a_emissions{std::vector<Emit>(nHets)};
  std::vector<Emit> het_emissions{std::vector<Emit>(nHets)};
  Matrix<HaploProbs> het_sample_emissions{nSamples, nHets};

  Matrix<HaploProbs> forwards{nSamples, nLoci};
  Matrix<HaploProbs> backwards{nSamples, nLoci};
  Matrix<HaploProbs> posterior{nSamples, nLoci};
  Matrix<double> p_flip{nSamples, nHets};
  std::vector<double> flip_diffs{std::vector<double>(nHets)};
  std::vector<double> sample_probs{std::vector<double>(nSamples)};
  std::vector<double> p_mom_as{};

  static constexpr double p_mom_a_tol{0.001};
  double flip_threshold{log(1.001)};
  static constexpr uint64_t first_flip_iteration{10};
  static constexpr uint64_t min_iterations{first_flip_iteration + 10};
  static constexpr uint64_t max_iterations{500};
  uint64_t iteration{0};
  bool converged{false};

  uint64_t n_flips{0};
  uint64_t flip_index{0};
  double flip_diff{0};
  double p_mom_a_diffs{0};

  // Check this carefully!
  std::vector<LocusInfo> get_loci() const {
    std::vector<LocusInfo> result;
    // Get full list of sites of both locus types
    for (uint64_t f{0}, h{0};
         f != fragments.n_fragments() || h != hets.n_hets();) {
      if (f != fragments.n_fragments() && h != hets.n_hets() &&
          hets[h].position() == fragments[f].start()) {
        result.emplace_back(h++, true);
        ++f;
      } else if (f == fragments.n_fragments() ||
                 (h != hets.n_hets() &&
                  hets[h].position() < fragments[f].start())) {
        result.emplace_back(h++, true);
      } else {
        result.emplace_back(f++, false);
      }
    }
    return result;
  }

  unsigned int locus_position(const uint64_t index) const {
    return actual_locus_position(index);
    // return fragments[loci[index].fragment].start();  OLD
  }
  unsigned int actual_locus_position(const uint64_t index) const {
    const LocusInfo & locus{loci[index]};
    if (locus.is_het) {
      return hets[locus.index].position();
    } else {
      return fragments[locus.index].start();
    }
  }
  std::vector<unsigned int> get_distances() const {
    // Get locus distances
    std::vector<unsigned int> result;
    result.reserve(loci.size());
    unsigned int last_position{locus_position(0)};
    result.push_back(0);
    for (uint64_t index{1}; index != loci.size(); ++index) {
      const unsigned int position{locus_position(index)};
      result.push_back(position - last_position);
      last_position = position;
    }
    return result;
  }

  std::vector<double> get_coverages() const {
    std::vector<double> result(nHets);
    const LowessCoverageLookup coverage_lookup{fragments,
          coverage_hmm.coverages};
    const PreMappedVector<FragmentInfo> & all_fragments{
      fragments.all_fragments()};
    for (uint64_t h{0}; h != nHets; ++h) {
      auto found = std::upper_bound(all_fragments.begin(), all_fragments.end(),
                                    hets[h].position(),
                                    [](const unsigned int pos,
                                       const FragmentInfo & frag) {
                                      return pos < frag.start();
                                    });
      if (found == all_fragments.begin())
        throw Error("Unexpected begin() result in get_coverages()");
      result[h] = coverage_lookup((--found)->length());
    }
    return result;
  }

 public:
  std::vector<double> afterburner(const std::string & chr_name,
                                  const Reference & ref,
                                  const ChromosomeIndexLookup & lookup,
                                  const double tolerance) const {
    const unsigned int chr{lookup[chr_name]};
    // measure plots
    const Marker flip_marker{paa::circle(), 0.2, "0 0 0", 1, true, "0 0 0"};
    PSPage * const measure_page{new PSPage{plots,
            "Measure for " + chr_name, "1 3 (0.5 0.25)"}};
    ownp(measure_page);
    PSGraph * const measure_graph{new PSGraph{
        *measure_page, ";Position;Measure",
            Bounds{0.0, 1.0 * ref.size(chr), -0.5, 0.5}}};
    ownp(measure_graph);
    PSXYSeries * const measure_series{new PSXYSeries{
        *measure_graph, flip_marker}};
    ownp(measure_series);

    const Marker flip_marker2{paa::circle(), 0.2, "1 0 0", 1, true, "1 0 0"};
    PSGraph * const flip_measure_graph{new PSGraph{*measure_page,
            ";Position;Flip Difference",
            Bounds{0.0, 1.0 * ref.size(chr), 2 * tolerance, 10.0}}};
    std::ostringstream tol_ps;
    tol_ps << "0 0 1 c 2 lw np 0 xfc " << tolerance << " yc m "
           << " 1 xfc " << tolerance << " yc l sp ";
    flip_measure_graph->ps(tol_ps.str());
    ownp(flip_measure_graph);
    PSXYSeries * const flip_measure_series{new PSXYSeries{
        *flip_measure_graph, flip_marker2}};
    ownp(flip_measure_series);

    PSGraph * const p_mom_a_graph{new PSGraph{*measure_page,
            ";Position;Probability that Mom is Allele A",
            Bounds{0.0, 1.0 * ref.size(chr), -0.2, 1.2}}};
    ownp(p_mom_a_graph);
    PSXYSeries * const p_mom_a_series{new PSXYSeries{
        *p_mom_a_graph, flip_marker}};
    ownp(p_mom_a_series);

    std::random_device rd{};
    std::mt19937_64 mersenne{rd()};
    std::function<double()> gen{std::bind(
        std::uniform_real_distribution<double>(-0.1, 0.1), mersenne)};

    std::vector<std::vector<Call>> calls(1);

    for (unsigned int h{0}; h != hets.n_hets(); ++h) {
      const double prob{hets.flip(h) ? 1 - exp(p_mom_as[h]) : exp(p_mom_as[h])};
      if (prob < 0.5 || prob > 0.5) {
        const Call call{h, chr, hets[h].position(), prob};
        if (call.confidence_score() >= 7) {
          p_mom_a_series->add_point(call.pos, call.p_mom_a + gen());
        }
        if (flip_diffs[h] > 2 * tolerance)
          flip_measure_series->add_point(call.pos, flip_diffs[h]);
        if (flip_diffs[h] > tolerance) {
          if (calls.back().size()) calls.emplace_back();
          continue;
        }
        calls.back().push_back(call);
      }
    }

    std::vector<std::vector<Call>> test_calls;
    for (uint64_t b{0}; b != calls.size(); ++b) {
      const std::vector<Call> & s_calls{calls[b]};
      std::vector<Call> t_calls;
      for (const Call & call : s_calls)
        if (call.confidence_score() >= 7)
          t_calls.push_back(call);
      if (t_calls.size() > 200) {
        test_calls.resize(test_calls.size() + 1);
        test_calls.back().swap(t_calls);
      }
    }

    std::ostringstream calls_ps;

    uint64_t color_index{0};
    std::vector<unsigned int> flip_points;
    for (unsigned int b{0}; b + 1 < test_calls.size(); ++b) {
      int64_t opp{0};
      int64_t val{0};
      const std::vector<Call> & left_calls{test_calls[b]};
      const std::vector<Call> & right_calls{test_calls[b + 1]};
      for (uint64_t s{0}; s != hets.n_samples(); ++s) {
        for (const Call & left : left_calls) {
          const uint64_t id1{left.id};
          const uint64_t obs1{hets.unflipped(s, id1)};
          if (obs1 == 0 || obs1 == 3) continue;
          const bool mom_a1{left.p_mom_a > 0.5};
          for (const Call & right : right_calls) {
            const uint64_t id2{right.id};
            const uint64_t obs2{hets.unflipped(s, id2)};
            if (obs2 == 0 || obs2 == 3) continue;
            const bool mom_a2{right.p_mom_a > 0.5};
            ++opp;
            if ((obs1 == obs2 && mom_a1 == mom_a2) ||
                (obs1 != obs2 && mom_a1 != mom_a2)) {
              ++val;
            } else {
              --val;
            }
          }
        }
      }
      std::cout << chr_name << " " << opp << " " << val << " "
                << 1.0 * val / opp << std::endl;
      const double measure{1.0 * val / opp};
      const double height_scale{1.0 / 10000 / 2};
      const double lh{left_calls.size() * height_scale};
      const double rh{right_calls.size() * height_scale};
      calls_ps << colors[color_index++] << " c 2 lw np "
          // left line
               << left_calls.front().pos << " " << measure << " gc m "
               << left_calls.back().pos << " " << measure << " gc l "
          // right line
               << right_calls.front().pos << " " << measure << " gc m "
               << right_calls.back().pos << " " << measure << " gc l "
          // trapezoid
               << left_calls.back().pos << " " << measure - lh << " gc m "
               << right_calls.front().pos << " " << measure - rh << " gc l "
               << right_calls.front().pos << " " << measure + rh << " gc l "
               << left_calls.back().pos << " " << measure + lh << " gc l cp "
               << "sp\n";
      if (measure < 0) {
        flip_points.push_back(
            (left_calls.back().id + right_calls.front().id) / 2);
      }
    }
    measure_graph->ps(calls_ps.str());
    std::vector<double> new_phase{p_mom_as};
    for (const unsigned int f : flip_points)
      for (unsigned int h{f}; h != hets.n_hets(); ++h)
        new_phase[h] = log(1 - exp(new_phase[h]));
    return new_phase;
  }
};


// phase = afterburner(haha_hmm, tolerance);


}  // namespace paa

#endif  // PAA_HAHA_H_
