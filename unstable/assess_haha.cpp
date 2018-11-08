//
// assess_haha
//
// How good was a haha run?
//
// Copyright 2018 Peter Andrews @ CSHL
//

#include <algorithm>
#include <exception>
#include <iostream>
#include <fstream>
#include <future>
#include <memory>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

#include "error.h"
#include "genes.h"
#include "haha.h"
#include "mumdex.h"
#include "psplot.h"
#include "threads.h"

using std::cout;
using std::cerr;
using std::endl;
using std::exception;
using std::future;
using std::ifstream;
using std::make_unique;
using std::mutex;
using std::ostringstream;
using std::set;
using std::string;
using std::to_string;
using std::unique_lock;
using std::unique_ptr;
using std::vector;

using paa::mkdir;
using paa::readable;
using paa::Bounds;
using paa::ChromosomeIndexLookup;
using paa::Error;
using paa::Gap;
using paa::Gaps;
using paa::HetTable;
using paa::Marker;
using paa::PSDoc;
using paa::PSGraph;
using paa::PSPage;
using paa::PSXYSeries;
using paa::Reference;
using paa::ThreadPool;

  const vector<string> colors{
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

class DerivativeGraphs {
 public:
  template <class Container>
  DerivativeGraphs(PSDoc & doc,
                   const Container & original_values,
                   const unsigned int n_derivatives) {
    Container values{original_values.begin(), original_values.end()};
    Container derivative(values.size());
    const Marker marker{paa::circle(), 0.5, "0 0 0", 1, true, "0 0 0"};
    if (values.size() < n_derivatives + 1)
      throw Error("Too little data for derivatives");
    for (uint64_t d{0}; d != n_derivatives + 1; ++d) {
      PSPage * page{new PSPage{doc, "Derivative " + to_string(d) +
              " of Flip Improvements", "1 2"}};
      ownp(page);
      PSGraph * graph{new PSGraph{*page, ";Flip Number;Value"}};
      ownp(graph);
      PSXYSeries * series{new PSXYSeries{*graph, marker}};
      ownp(series);
      PSGraph * zgraph{new PSGraph{*page, ";Flip Number;Value"}};
      ownp(zgraph);
      PSXYSeries * zseries{new PSXYSeries{*zgraph, marker}};
      ownp(zseries);
      for (uint64_t f{d}; f != values.size(); ++f) {
        const typename Container::value_type val{values[f]};
        series->add_point(f, val);
        if (f < 100) zseries->add_point(f, val);
        derivative[f] = f ? val - values[f - 1] : 0;
      }
      values.swap(derivative);
    }
  }
};

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

struct CallPair {
  CallPair(const Call & call1, const Call & call2) :
      confidence{std::max(call1.confidence(), call2.confidence())},
    distance{call1.pos > call2.pos ? call1.pos - call2.pos :
        call2.pos - call1.pos},
    agree{call1.agree(call2)} { }
  double confidence;
  unsigned int distance;
  unsigned int agree;
  static constexpr unsigned int max_confidence_bin{10};
  static constexpr unsigned int distance_subdivision{2};
  static unsigned int confidence_bin(const double confidence) {
    return std::min(1.0 * max_confidence_bin, -log10(confidence));
  }
  unsigned int confidence_bin() const {
    return confidence_bin(confidence);
  }
  static unsigned int distance_bin(const double distance) {
    return log10(distance) * distance_subdivision;
  }
  unsigned int distance_bin() const {
    return distance_bin(distance);
  }
  static double bin_distance(const unsigned int dist_bin) {
    return pow(10.0, 1.0 * dist_bin / distance_subdivision);
  }
};

class Counts {
 public:
  Counts() : confidence_hist(CallPair::max_confidence_bin + 1),
             correct_counts(confidence_hist.size()),
             call_counts(confidence_hist.size()) { }
  Counts & operator+=(const Counts & other) {
    chr_boundaries.push_back(all_calls.size());
    all_calls.insert(all_calls.end(),
                     other.all_calls.begin(), other.all_calls.end());
    bin_counts.resize(std::max(bin_counts.size(),
                               other.bin_counts.size()));
    agree_counts.resize(std::max(agree_counts.size(),
                                 other.agree_counts.size()));
    bin_counts_dist.resize(std::max(bin_counts_dist.size(),
                                    other.bin_counts_dist.size()));
    agree_counts_dist.resize(std::max(agree_counts_dist.size(),
                                      other.agree_counts_dist.size()));
    for (unsigned int i{0};
         i != std::min(bin_counts.size(), other.bin_counts.size());
         ++i) {
      bin_counts[i] += other.bin_counts[i];
      agree_counts[i] += other.agree_counts[i];
      bin_counts_dist[i].resize(std::max(bin_counts_dist[i].size(),
                                         other.bin_counts_dist[i].size()));
      agree_counts_dist[i].resize(std::max(agree_counts_dist[i].size(),
                                           other.agree_counts_dist[i].size()));
      for (unsigned int j{0};
           j != std::min(bin_counts_dist[i].size(),
                         other.bin_counts_dist[i].size());
           ++j) {
        bin_counts_dist[i][j] += other.bin_counts_dist[i][j];
        agree_counts_dist[i][j] += other.agree_counts_dist[i][j];
      }
    }
    for (unsigned int b{0}; b != confidence_hist.size(); ++b)
      confidence_hist[b] += other.confidence_hist[b];
    distance_hist.resize(std::max(distance_hist.size(),
                                  other.distance_hist.size()));
    for (unsigned int b{0};
         b != std::min(distance_hist.size(), other.distance_hist.size());
         ++b)
      distance_hist[b] += other.distance_hist[b];
    segment_hist.resize(std::max(segment_hist.size(),
                                  other.segment_hist.size()));
    for (unsigned int b{0};
         b != std::min(segment_hist.size(), other.segment_hist.size());
         ++b)
      segment_hist[b] += other.segment_hist[b];

    for (unsigned int b{0}; b != call_counts.size(); ++b) {
      correct_counts[b] += other.correct_counts[b];
      call_counts[b] += other.call_counts[b];
    }
    return *this;
  }
  void flip_calls_if_needed() {
    const uint64_t tot_call{
      accumulate(call_counts.begin(), call_counts.end(), 0ul)};
    const uint64_t correct_call{
      accumulate(correct_counts.begin(), correct_counts.end(), 0ul)};
    if (correct_call * 2 < tot_call) {
      for (unsigned int b{0}; b != call_counts.size(); ++b) {
        correct_counts[b] = call_counts[b] - correct_counts[b];
      }
      for (unsigned int p{0}; p != all_calls.size(); ++p) {
        all_calls[p].p_mom_a = 1 - all_calls[p].p_mom_a;
      }
    }
  }
  void add_confidence(const double confidence) {
    ++confidence_hist[CallPair::confidence_bin(confidence)];
  }
  void add_correct(const double confidence, const bool correct) {
    const uint64_t confidence_bin{CallPair::confidence_bin(confidence)};
    ++call_counts[confidence_bin];
    if (correct) ++correct_counts[confidence_bin];
  }
  void increment_distance_bin(const double distance_bin) {
    if (distance_bin >= distance_hist.size())
      distance_hist.resize(distance_bin + 1);
    ++distance_hist[distance_bin];
  }
  void increment_segment_bin(const double distance_bin) {
    if (distance_bin >= segment_hist.size())
      segment_hist.resize(distance_bin + 1);
    ++segment_hist[distance_bin];
  }

  vector<uint64_t> confidence_hist{};
  vector<uint64_t> distance_hist{};
  vector<uint64_t> segment_hist{};
  vector<uint64_t> bin_counts{};
  vector<uint64_t> agree_counts{};
  vector<uint64_t> correct_counts{};
  vector<uint64_t> call_counts{};
  vector<Call> all_calls{};
  vector<uint64_t> chr_boundaries{};
  vector<vector<uint64_t>> bin_counts_dist{};
  vector<vector<uint64_t>> agree_counts_dist{};
};

const bool show_flip_graphs{false};

struct KnownSegment {
  explicit KnownSegment(const uint64_t start_, const bool split_ = false) :
      start{start_}, split{split_} { }
  KnownSegment(const uint64_t start_, const uint64_t stop_,
               const bool split_, const bool increasing_) :
      start{start_}, stop{stop_}, split{split_}, increasing{increasing_} { }
  uint64_t start{0};
  uint64_t stop{0};
  bool split{false};
  bool increasing{false};
};

int main(int argc, char * argv[]) try {
  // cout << min_confidence << '\n';
  // return 0;
  const std::string usage{"assess_haha ref data_dir flip_tol phase_files ..."};
  --argc;
  if (argc < 4) throw Error(usage);
  const Reference ref{argv[1]};
  const ChromosomeIndexLookup chr_lookup{ref};
  const Gaps gaps{ref.fasta_file() + ".bin/gaps.txt",
        chr_lookup, "centromere"};
  const string data_dir{argv[2]};
  const double flip_tol{strtod(argv[3], nullptr)};

  mutex plots_mutex;
  auto chr_fun = [&ref, &data_dir, &chr_lookup, &gaps, flip_tol, &plots_mutex]
      (const string & phase_name, PSDoc & plots) {
    // Open phase file
    ifstream phase_file{phase_name.c_str()};
    if (!phase_file) throw Error("Problem opening phase file") << phase_name;

    // Determine chromosome
    const size_t chr_idx{phase_name.find("chr")};
    const size_t period{phase_name.find(".phase")};
    const string chr_name{phase_name.substr(chr_idx + 3, period - chr_idx - 3)};
    const unsigned int chr{chr_lookup["chr" + chr_name]};

    unique_lock<mutex> plots_lock{plots_mutex};
    // Series to show flips
    PSGraph * const flip_graph{show_flip_graphs ?
          new PSGraph{plots,
            string("Flip diff for ") + phase_name + ";Position;Flip Diff"} :
      nullptr};
    ownp(flip_graph);
    const Marker flip_marker{paa::circle(), 0.2, "0 0 0", 1, true, "0 0 0"};
    PSXYSeries * const flip_series{show_flip_graphs ?
          new PSXYSeries{*flip_graph, flip_marker} : nullptr};
    ownp(flip_series);
    PSPage * const measure_page{new PSPage{plots,
            "Measure for chr" + chr_name, "1 3 (0.5 0.25)"}};
    ownp(measure_page);
    PSGraph * const measure_graph{new PSGraph{*measure_page,
            ";Position;Measure",
            Bounds{0.0, 1.0 * ref.size(chr), -0.5, 0.5}}};
    ownp(measure_graph);
    PSXYSeries * const measure_series{new PSXYSeries{
        *measure_graph, flip_marker}};
    ownp(measure_series);
    const Marker flip_marker2{paa::circle(), 0.2, "1 0 0", 1, true, "1 0 0"};
    PSGraph * const flip_measure_graph{new PSGraph{*measure_page,
            ";Position;Flip Difference",
            Bounds{0.0, 1.0 * ref.size(chr), 2 * flip_tol, 10.0}}};
    ostringstream tol_ps;
    tol_ps << "0 0 1 c 2 lw np 0 xfc " << flip_tol << " yc m "
    << " 1 xfc " << flip_tol << " yc l sp ";
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
    plots_lock.unlock();

    // Get centromere start and stop
    const Gap & gap{[&gaps, chr]() {
        for (const Gap & result : gaps) if (result.chr == chr) return result;
        throw Error("Gap not found");
      }()};
    const unsigned int gap_start{gap.start};
    const unsigned int gap_stop{gap.stop};

    // Load Het table for chromosome
    const string table_name{data_dir + "/hets_chr" + ref.name(chr) + ".txt"};
    const HetTable hets{table_name, chr_lookup};

    std::random_device rd{};
    std::mt19937_64 mersenne{rd()};
    std::function<double()> gen{std::bind(
        std::uniform_real_distribution<double>(-0.1, 0.1), mersenne)};

    unsigned int id;
    unsigned int pos;
    char phase1;
    char phase2;
    double p_mom_a;
    double confidence;
    double prob_diff;
    vector<vector<Call>> calls(1);
    Counts counts;
    while (phase_file >> id >> pos >> phase1 >> phase2
           >> p_mom_a >> confidence >> prob_diff) {
      try {
        const Call call{id, chr, pos, p_mom_a};
        if (false) cout << call.pos
                        << '\t' << call.p_mom_a
                        << '\t' << call.confidence()
                        << '\t' << call.confidence_score() << '\n';
        if (show_flip_graphs) flip_series->add_point(pos, prob_diff);
        counts.add_correct(call.confidence(), p_mom_a >= 0.5);
        if (call.confidence_score() >= 7) {
          counts.all_calls.push_back(call);
          p_mom_a_series->add_point(call.pos, call.p_mom_a + gen());
        }
        if (prob_diff > 2 * flip_tol)
          flip_measure_series->add_point(call.pos, prob_diff);
        if (prob_diff > flip_tol) {
          if (calls.back().size()) calls.emplace_back();
          continue;
        }
        if (false && call.pos >= gap_start && call.pos <= gap.stop) continue;
        calls.back().push_back(call);
        counts.add_confidence(call.confidence());
      } catch (...) {}
    }

    for (uint64_t c{0}; c != calls.size(); ++c) {
      if (calls[c].size() > 1) {
        const unsigned int distance_bin{CallPair::distance_bin(
            calls[c].back().pos - calls[c].front().pos)};
        counts.increment_segment_bin(distance_bin);
      }
    }
    if (false) cout << chr_name << '\t' << chr << '\t'
                    << gap_start << '\t' << gap_stop << '\n';
    for (uint64_t k{0}; k != calls.size(); ++k) {
      const vector<Call> & arm_calls{calls[k]};
      for (uint64_t i{0}; i != arm_calls.size(); ++i) {
        for (uint64_t j{i + 1}; j != arm_calls.size(); ++j) {
          const CallPair call_pair{arm_calls[i], arm_calls[j]};
          const unsigned int confidence_bin{call_pair.confidence_bin()};
          if (confidence_bin >= counts.bin_counts.size()) {
            counts.bin_counts.resize(confidence_bin + 1);
            counts.agree_counts.resize(confidence_bin + 1);
            counts.bin_counts_dist.resize(confidence_bin + 1);
            counts.agree_counts_dist.resize(confidence_bin + 1);
          }
          ++counts.bin_counts[confidence_bin];
          const unsigned int distance_bin{call_pair.distance_bin()};
          counts.increment_distance_bin(distance_bin);
          if (distance_bin >= counts.bin_counts_dist[confidence_bin].size()) {
            counts.bin_counts_dist[confidence_bin].resize(distance_bin + 1);
            counts.agree_counts_dist[confidence_bin].resize(distance_bin + 1);
          }
          ++counts.bin_counts_dist[confidence_bin][distance_bin];
          if (call_pair.agree) {
            ++counts.agree_counts[confidence_bin];
            ++counts.agree_counts_dist[confidence_bin][distance_bin];
          }
        }
      }
    }

    vector<vector<Call>> test_calls;
    for (uint64_t b{0}; b != calls.size(); ++b) {
      const vector<Call> & s_calls{calls[b]};
      vector<Call> t_calls;
      for (const Call & call : s_calls)
        if (call.confidence_score() >= 7)
          t_calls.push_back(call);
      if (t_calls.size() > 200) {
        test_calls.resize(test_calls.size() + 1);
        test_calls.back().swap(t_calls);
      }
    }

    ostringstream calls_ps;
#if 0
    const uint64_t call_block{250};
    vector<call> & ac{counts.all_calls};
    for (uint64_t le{1}; le != ac.size(); ++le) {
      const uint64_t lb{le > call_block ? le - call_block : 0};
      const uint64_t rb{le};
      const uint64_t re{std::min(rb + call_block, ac.size())};
    }
#endif

    uint64_t color_index{0};
    for (unsigned int b{0}; b + 1 < test_calls.size(); ++b) {
      int64_t opp{0};
      int64_t val{0};
      const vector<Call> & left_calls{test_calls[b]};
      const vector<Call> & right_calls{test_calls[b + 1]};
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
      cout << chr_name << " " << opp << " " << val << " " << 1.0 * val / opp
           << endl;
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
    }
    measure_graph->ps(calls_ps.str());
    cerr << chr_name << endl;
    counts.flip_calls_if_needed();
    return counts;
    };

  PSDoc plots{"assess"};
  plots.pdf(false);

  PSGraph assess{plots, "HaHa Assessment;Distance;% Correct"};
  assess.log_x(true);

  PSGraph flips_graph{plots, "Segments unlikely to flip;Length;N"};
  PSXYSeries flips{flips_graph};
  flips_graph.log_x(true);
  flips_graph.range().yl(0);

  PSXYSeries distance{plots, "Assessed HaHa Het Distances;log10(Distance);N"};
  distance.parents().front()->log_y(true);

  PSXYSeries confidence{
    plots, "HaHa Phase Call Confidence Scores;Confidence Score;N"};

  PSXYSeries correct{
    plots, "HaHa Phasing Evaluation;Confidence Score;% Calls Correct"};

  const Marker calls_marker{paa::circle(), 0.2, "0 0 0", 1, true, "0 0 0"};
  PSGraph calls_graph{plots,
        "HaHa Phasing Evaluation;Het #;Cumulative Correct - Incorrect"};
  calls_graph.do_ticks(false);
  PSXYSeries calls{calls_graph, calls_marker};

  PSGraph flip_improve_graph{plots,
        "Flip Improvements to Truth;N Flips;% Correct"};
  const Marker bad_marker{paa::circle(), 0.2, "1 0 0", 1, true, "1 0 0"};
  PSXYSeries bad_flip_improve{flip_improve_graph, bad_marker};
  const Marker best_marker{paa::circle(), 0.2, "0 1 0", 1, true, "0 1 0"};
  PSXYSeries best_flip_improve{flip_improve_graph, best_marker};

  // Load haha phase data
  ThreadPool pool{24};
  std::vector<std::future<Counts>> futures;
  argc -= 2; argv += 2;
  while (--argc && ++argv) {
    sleep(1);
    futures.push_back(pool.run(chr_fun, argv[1], std::ref(plots)));
  }

  // Retrieve call and count info from threads
  const Counts counts{[&futures]() {
      Counts result;
      for (auto & fut : futures) result += fut.get();
      return result;
    }()};

  cerr << "plots" << endl;
  for (unsigned int b{0}; b != counts.call_counts.size(); ++b)
    correct.add_point(b, 100.0 * counts.correct_counts[b] /
                      counts.call_counts[b]);
  for (unsigned int b{0}; b != counts.confidence_hist.size(); ++b)
    confidence.add_point(b, counts.confidence_hist[b]);
  for (unsigned int b{0}; b != counts.distance_hist.size(); ++b)
    if (counts.distance_hist[b])
      distance.add_point(b / 2.0, counts.distance_hist[b]);
  for (unsigned int b{0}; b != counts.segment_hist.size(); ++b)
    flips.add_point(CallPair::bin_distance(b), counts.segment_hist[b]);

  // Bad flipping method
  cerr << "bad flip" << endl;
  int64_t count{0};
  vector<int64_t> correct_counts(counts.all_calls.size());
  for (unsigned int p{0}; p != counts.all_calls.size(); ++p) {
    const Call & call{counts.all_calls[p]};
    if (call.p_mom_a >= 0.5) {
      ++count;
    } else {
      --count;
    }
    correct_counts[p] = count;
    if (p % 1 == 0) calls.add_point(p, count);
  }
  ostringstream calls_ps;
  // Put centromeres on graph - an ugly hack
  set<unsigned int> chromosomes;
  for (const Call & call : counts.all_calls) chromosomes.insert(call.chr);
  calls_ps << "0.7 0.7 0.7 c 0.5 lw\n";
  for (const unsigned int chr : chromosomes) {
    const Gap & gap{[&gaps, chr]() {
        for (const Gap & result : gaps)
          if (result.chr == chr)
            return result;
        throw Error("Gap not found");
      }()};
    for (unsigned int pos : {gap.start, gap.stop}) {
      uint64_t idx = lower_bound(
          counts.all_calls.begin(), counts.all_calls.end(), pos,
          [chr](const Call & call, const unsigned int p) {
            if (chr == call.chr) {
              return call.pos < p;
            } else {
              return call.chr < chr;
            }
          }) - counts.all_calls.begin();
      if (idx == counts.all_calls.size()) continue;
      calls_ps << "np " << idx << " xc " << 0 << " yfc m "
               << idx << " xc " << 1 << " yfc l sp\n";
    }
  }
  // Chromosome boundaries
  calls_ps << "1 lw 8 sf\n";
  for (unsigned int b{0}; b != counts.chr_boundaries.size(); ++b) {
    if (b) calls_ps << "np 1 0 0 c "
                    << counts.chr_boundaries[b] << " xc 0 yfc m "
                    << counts.chr_boundaries[b] << " xc 1 yfc l sp\n";
    const uint64_t end{b + 1 == counts.chr_boundaries.size() ?
          correct_counts.size() : counts.chr_boundaries[b + 1]};
    const uint64_t mid{(counts.chr_boundaries[b] + end) / 2};
    calls_ps << "0 0 0 c " << mid << " xc " << 0.01 << " yfc m "
             << "(" << b + 1 << ") jc s\n";
  }
  calls_graph.ps(calls_ps.str());
  vector<uint64_t> decreasing_intervals(1);
  int64_t last_count{0};
  for (uint64_t p{0}; p != correct_counts.size(); ++p) {
    if (last_count > correct_counts[p]) {
      ++decreasing_intervals.back();
    } else {
      if (decreasing_intervals.back() != 0) {
        decreasing_intervals.push_back(0);
      }
    }
    last_count = correct_counts[p];
  }
  sort(decreasing_intervals.begin(), decreasing_intervals.end(),
       [](const uint64_t lhs, const uint64_t rhs) {
         return rhs < lhs;
       });
  bad_flip_improve.add_point(0, 100.0 * (counts.all_calls.size() + count) /
                             (2 * counts.all_calls.size()));
  for (unsigned int f{0}; f != decreasing_intervals.size(); ++f) {
    count += 2 * decreasing_intervals[f];
    const double per_correct{100.0 * (counts.all_calls.size() + count) /
          (2 * counts.all_calls.size())};
    bad_flip_improve.add_point(f + 1, per_correct);
  }

  // Best flipping method
  cerr << "good flip" << endl;
  const uint64_t n_derivatives{4};
  vector<vector<double>> flip_derivatives(n_derivatives + 2);
  best_flip_improve.add_point(
      0, 100.0 * (correct_counts.size() + correct_counts.back()) /
      (2 * correct_counts.size()));
  uint64_t flip_number{0};
  ostringstream flips_ps;
  flips_ps << "10 sf 0 0 1 c\n";
  vector<KnownSegment> known_segments;
  vector<int64_t> orig_counts{correct_counts.begin(), correct_counts.end()};
  while (correct_counts.back() != static_cast<int64_t>(correct_counts.size())) {
    uint64_t low_index{correct_counts.size() - 1};
    uint64_t best_index{low_index};
    uint64_t best_low_index{low_index};
    int64_t best_diff{0};
    for (uint64_t fc{0}; fc != correct_counts.size(); ++fc) {
      const uint64_t rc{correct_counts.size() - fc - 1};
      if (correct_counts[rc] < correct_counts[low_index]) {
        low_index = rc;
      }
      const int64_t diff{correct_counts[rc] - correct_counts[low_index]};
      if (diff > best_diff) {
        best_diff = diff;
        best_index = rc;
        best_low_index = low_index;
      }
    }
    for (uint64_t i{best_index + 1}; i != correct_counts.size(); ++i) {
      if (i < best_low_index) {
        correct_counts[i] = 2 * correct_counts[best_index] -
            correct_counts[i];
      } else {
        correct_counts[i] += 2 * best_diff;
      }
    }
    const double val{100.0 * (correct_counts.size() + correct_counts.back()) /
          (2 * correct_counts.size())};
    best_flip_improve.add_point(++flip_number, val);
    flip_derivatives[0].push_back(val);

    if (flip_number <= 40) {
      const double yoff{orig_counts.back() / 200.0};
      flips_ps << best_index << " "
               << orig_counts[best_index] + yoff << " gc m "
               << "(" << flip_number << ") jc s\n";
      known_segments.emplace_back(best_index);
      known_segments.emplace_back(best_low_index);
    }
  }
  calls_graph.ps(flips_ps.str());

  for (uint64_t d{0}; d != 1; ++d) {
    const vector<double> & values{flip_derivatives[d]};
    for (uint64_t f{0}; f != values.size(); ++f) {
      if (f) {
        flip_derivatives[d + 1].push_back(values[f] - values[f - 1]);
      } else {
        flip_derivatives[d + 1].push_back(0);
      }
    }
  }
  if (false) {
    const DerivativeGraphs der_graphs{
      plots, flip_derivatives[0], n_derivatives};
  }
  cerr << "confidence" << endl;
  for (unsigned int bin{0}; bin != counts.bin_counts.size(); ++bin) {
    const double agree_frac{
      1.0 * counts.agree_counts[bin] / counts.bin_counts[bin]};
    cout << bin
         << '\t' << counts.bin_counts[bin]
         << '\t' << counts.agree_counts[bin]
         << '\t' << agree_frac
         << '\n';
    const Marker conf_marker{
      paa::circle(), 1, colors[bin], 1, true, colors[bin]};
    PSXYSeries * confidence_series{new PSXYSeries{assess, conf_marker}};
    ownp(confidence_series);
    for (unsigned int dist_bin{0};
         dist_bin != counts.agree_counts_dist[bin].size();
         ++dist_bin) {
      const double agreement{counts.bin_counts_dist[bin][dist_bin] ?
            1.0 * counts.agree_counts_dist[bin][dist_bin] /
            counts.bin_counts_dist[bin][dist_bin] : 0.0};
      cout << '\t' << dist_bin
           << '\t' << CallPair::bin_distance(dist_bin)
           << '\t' << counts.bin_counts_dist[bin][dist_bin]
           << '\t' << counts.agree_counts_dist[bin][dist_bin]
           << '\t' << agreement
           << '\n';
      confidence_series->add_point(CallPair::bin_distance(dist_bin),
                                   100 * agreement);
    }
  }

  // Compare adjacent segments
  cerr << "comparison" << endl;
  for (const uint64_t bound : counts.chr_boundaries)
    known_segments.emplace_back(bound);
  auto sort_start = [](const KnownSegment & lhs, const KnownSegment & rhs) {
    return lhs.start < rhs.start;
  };
  sort(known_segments.begin(), known_segments.end(), sort_start);
  auto start_eq = [](const KnownSegment & lhs, const KnownSegment & rhs) {
    return lhs.start == rhs.start;
  };
  auto end = unique(known_segments.begin(), known_segments.end(), start_eq);
  known_segments.erase(end, known_segments.end());
  // Split big segments
  const uint64_t target_seg_size{500};
  const uint64_t stop_seg{known_segments.size()};
  for (uint64_t i{0}; i + 1 != stop_seg; ++i) {
    KnownSegment & seg{known_segments[i]};
    KnownSegment & nseg{known_segments[i + 1]};
    const uint64_t seg_size{nseg.start - seg.start};
    if (seg_size > target_seg_size * 1.5) {
      seg.split = true;
      const uint64_t n_new_size{static_cast<uint64_t>(
          1.0 * seg_size / (seg_size / target_seg_size))};
      cerr << i << " "
           << seg.start << " " << nseg.start << " "
           << seg_size << " " << n_new_size;
      uint64_t ostart{seg.start};
      uint64_t start{seg.start};
      while ((start += n_new_size) + target_seg_size / 2 < nseg.start) {
        known_segments.emplace_back(start, true);
        cerr << " " << start - ostart;
        ostart = start;
      }
      cerr << " " << nseg.start - ostart;
      cerr << endl;
    }
  }
  sort(known_segments.begin(), known_segments.end(), sort_start);
  const uint64_t min_seg_size{250};
  vector<KnownSegment> ok_segments;
  for (uint64_t i{0}; i != known_segments.size(); ++i) {
    KnownSegment & seg{known_segments[i]};
    const uint64_t b{seg.start};
    const uint64_t e{i + 1 == known_segments.size() ?
          counts.all_calls.size() - 1 : known_segments[i + 1].start - 1};
    const Call & call{counts.all_calls[b]};
    const Call & ecall{counts.all_calls[e]};
    if (call.chr != ecall.chr) throw Error("Segment chromosome mismatch");
    if (e - b + 1 >= min_seg_size) {
      ok_segments.emplace_back(b, e, seg.split,
                               orig_counts[b] < orig_counts[e]);
      cout << "segs " << i << '\t' << b << '\t' << e << '\t'
           << e - b + 1 << '\t'
           << ref.name(call.chr) << '\t' << call.pos << '\t'
           << ref.name(ecall.chr) << '\t' << ecall.pos << '\t'
           << orig_counts[e] - orig_counts[b] << '\t'
           << 1.0 * (orig_counts[e] - orig_counts[b]) / (e - b) << '\t'
           << seg.split
           << endl;
    }
  }
  for (uint64_t i{0}; i != ok_segments.size(); ++i) {
    const KnownSegment & seg{ok_segments[i]};
    const Call & call{counts.all_calls[seg.start]};
    cout << i << " " << seg.start << " " << seg.stop << " "
         << seg.stop - seg.start << " "
         << seg.split << " " << seg.increasing << " "
         << ref.name(call.chr) << " " << call.pos << endl;
  }

  return 0;

  const vector<unique_ptr<HetTable>> het_tables{
    [&ref, &chr_lookup, &data_dir]() {
      vector<unique_ptr<HetTable>> result(ref.n_chromosomes());
      for (unsigned int c{0}; c != result.size(); ++c) {
        const string table_name{data_dir + "/hets_chr" + ref.name(c) + ".txt"};
        if (readable(table_name))
          result[c] = make_unique<HetTable>(table_name, chr_lookup);
      }
      return result;
    }()};


  PSGraph measure_graph{plots, "HaHa Assessment;Distance;Measure"};
  const Marker different_marker{paa::circle(), 0.4, "1 0 0", 1, true, "1 0 0"};
  PSXYSeries different_chrs{measure_graph, different_marker};
  const Marker same_marker{paa::circle(), 0.4, "0 1 0", 1, true, "0 1 0"};
  PSXYSeries same_chr{measure_graph, same_marker};
  const Marker flipped_marker{paa::circle(), 0.4, "0 0 1", 1, true, "0 0 1"};
  PSXYSeries flipped_chr{measure_graph, flipped_marker};
  const uint64_t n_adjacent_segments{20};
  for (uint64_t s1{0}; s1 != ok_segments.size(); ++s1) {
    const KnownSegment & seg1{ok_segments[s1]};
    const unsigned int chr1{counts.all_calls[seg1.start].chr};
    const HetTable & hets1{*het_tables[chr1]};
    for (uint64_t s2{s1 + 1};
         s2 != ok_segments.size() && s2 != s1 + n_adjacent_segments; ++s2) {
      const KnownSegment & seg2{ok_segments[s2]};
      const unsigned int chr2{counts.all_calls[seg2.start].chr};
      const HetTable & hets2{*het_tables[chr2]};
      int64_t opp{0};
      int64_t val{0};
      double dist{0};
      for (uint64_t s{0}; s != hets1.n_samples(); ++s) {
        for (uint64_t i{seg1.start}; i <= seg1.stop; ++i) {
          const uint64_t id1{counts.all_calls[i].id};
          const uint64_t obs1{hets1.unflipped(s, id1)};
          const unsigned int pos1{hets1[id1].position()};
          const bool mom_a1{counts.all_calls[i].p_mom_a > 0.5};
          if (obs1 == 0 || obs1 == 3) continue;
          for (uint64_t j{seg2.start}; j <= seg2.stop; ++j) {
            const uint64_t id2{counts.all_calls[j].id};
            const uint64_t obs2{hets2.unflipped(s, counts.all_calls[j].id)};
            const unsigned int pos2{hets2[id2].position()};
            const bool mom_a2{counts.all_calls[j].p_mom_a > 0.5};
            if (obs2 == 0 || obs2 == 3) continue;
            ++opp;
            dist += pos2 > pos1 ? pos2 - pos1 : pos1 - pos2;
            if ((obs1 == obs2 && mom_a1 == mom_a2) ||
                (obs1 != obs2 && mom_a1 != mom_a2)) {
              ++val;
            } else {
              --val;
            }
          }
        }
      }
      dist /= opp;
      const double measure{1.0 * val / opp};
      cerr << "Seg compare " << s1 << " " << s2
           << " " << chr1
           << " " << chr2
           << " " << (seg1.increasing == seg2.increasing)
           << " " << dist
           << " " << val << " " << opp << " " << measure << endl;
      if (chr1 == chr2) {
        if (seg1.increasing == seg2.increasing) {
          same_chr.add_point(dist, measure);
        } else {
          flipped_chr.add_point(dist, measure);
        }
      } else {
        different_chrs.add_point(dist, measure);
      }
    }
  }
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
