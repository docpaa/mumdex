//
// afterburner
//
// Fix up a haha run
//
// Copyright 2018 Peter Andrews @ CSHL
//

#include <algorithm>
#include <exception>
#include <iostream>
#include <fstream>
#include <memory>
#include <random>
#include <set>
#include <sstream>
#include <string>
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
using std::ifstream;
using std::make_unique;
using std::ostringstream;
using std::set;
using std::string;
using std::to_string;
using std::unique_ptr;
using std::vector;

using paa::mkdir;
using paa::readable;
using paa::Bounds;
using paa::ChromosomeIndexLookup;
using paa::Error;
using paa::HetTable;
using paa::Marker;
using paa::PSDoc;
using paa::PSGraph;
using paa::PSPage;
using paa::PSXYSeries;
using paa::Reference;

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

int main(int argc, char * argv[]) try {
  const std::string usage{"afterburner ref data_dir flip_tol chr"};
  if (--argc != 4) throw Error(usage);
  const Reference ref{argv[1]};
  const ChromosomeIndexLookup chr_lookup{ref};
  const string data_dir{argv[2]};
  const double flip_tol{strtod(argv[3], nullptr)};
  const string chr_name{argv[4]};
  const unsigned int chr{chr_lookup[chr_name]};

  ostringstream phase_name;
  phase_name << data_dir << "/" << chr_name << ".phase.txt";
  ifstream phase_file{phase_name.str().c_str()};
  if (!phase_file)
    throw Error("Problem opening phase file") << phase_name.str();

  PSDoc plots{"afterburner"};
  plots.pdf(true);

  const Marker flip_marker{paa::circle(), 0.2, "0 0 0", 1, true, "0 0 0"};

  PSPage * const measure_page{PSPage::create(plots, "Measure for " + chr_name,
                                             "1 3 (0.5 0.25)")};
  PSGraph * const measure_graph{PSGraph::create(
      *measure_page, ";Position;Measure",
      Bounds{0.0, 1.0 * ref.size(chr), -0.5, 0.5})};

  const Marker flip_marker2{paa::circle(), 0.2, "1 0 0", 1, true, "1 0 0"};
  PSGraph * const flip_measure_graph{PSGraph::create(
      *measure_page, ";Position;Flip Difference",
      Bounds{0.0, 1.0 * ref.size(chr), 2 * flip_tol, 10.0})};
  ostringstream tol_ps;
  tol_ps << "0 0 1 c 2 lw np 0 xfc " << flip_tol << " yc m "
         << " 1 xfc " << flip_tol << " yc l sp ";
  flip_measure_graph->ps(tol_ps.str());

  PSXYSeries * const flip_measure_series{PSXYSeries::create(
      *flip_measure_graph, flip_marker2)};

  PSGraph * const p_mom_a_graph{PSGraph::create(
      *measure_page, "HMM Output;Position;Probability that Mom is Allele A",
      Bounds{0.0, 1.0 * ref.size(chr), -0.2, 1.2})};
  PSXYSeries * const p_mom_a_series{PSXYSeries::create(
      *p_mom_a_graph, flip_marker)};

#if 0
  PSPage * const fixed_page{PSPage::create(
      plots, "Haha Afterburner Improvement for " + chr_name, "1 4")};
  p_mom_a_graph->add(fixed_page);
  PSXYSeries * const fixed_p_mom_a_series{PSXYSeries::create(
      *fixed_page,
      "Locus Afterburner;Position;Probability that Mom is Allele A",
      flip_marker)};
  PSXYSeries * const segment_p_mom_a_series{PSXYSeries::create(
      *fixed_page,
      "Segment Afterburner;Position;Probability that Mom is Allele A",
      flip_marker)};
  PSXYSeries * const window_p_mom_a_series{PSXYSeries::create(
      *fixed_page,
      "Growing Afterburner;Position;Probability that Mom is Allele A",
      flip_marker)};
#endif

  // Load Het table for chromosome
  const string table_name{data_dir + "/hets_" + ref.name(chr) + ".txt"};
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
  const unsigned int calls_per_bin{100};
  vector<vector<Call>> confident_calls(1);
  while (phase_file >> id >> pos >> phase1 >> phase2
         >> p_mom_a >> confidence >> prob_diff) {
    try {
      const Call call{id, chr, pos, p_mom_a};
      if (call.confidence_score() >= 7) {
        p_mom_a_series->add_point(call.pos, call.p_mom_a + gen());
        if (true || pos < 5000000) {
          confident_calls.back().push_back(call);
          if (confident_calls.back().size() == calls_per_bin)
            confident_calls.resize(confident_calls.size() + 1);
        }
      }
      if (prob_diff > 2 * flip_tol)
        flip_measure_series->add_point(call.pos, prob_diff);
      if (prob_diff > flip_tol) {
        if (calls.back().size()) calls.emplace_back();
        continue;
      }
      calls.back().push_back(call);
    } catch (...) {}
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
    const double measure{1.0 * val / std::max(opp, static_cast<int64_t>(1))};
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

  // afterburner heatmap...
  if (confident_calls.back().size() < calls_per_bin)
    confident_calls.pop_back();
  PSGraph * const afterburner_graph{PSGraph::create(
      plots, "Segment Concordance;Segment 1; Segment 2",
      Bounds{0.0, 1.0 * confident_calls.size(),
            0.0, 1.0 * confident_calls.size()})};

  PSXYSeries * const self_graph{PSXYSeries::create(
      plots, "Self-agreement;Fraction correct;Flip Measure", flip_marker)};

  uint64_t n_correct{0};
  uint64_t n_calls{0};
  for (unsigned int s1{0}; s1 != confident_calls.size(); ++s1) {
     const vector<Call> & left_calls{confident_calls[s1]};
     for (const Call & left : left_calls) {
       ++n_calls;
       const bool mom_a1{left.p_mom_a > 0.5};
       n_correct += mom_a1;
     }
  }
  cout << "initial fraction correct " << 1.0 * n_correct / n_calls << endl;

  unsigned int best_segment{0};
  double best_self_measure{0};
  ostringstream after_ps;
  vector<vector<double>> measures(confident_calls.size(),
                                  vector<double>(confident_calls.size()));
  for (unsigned int s1{0}; s1 != confident_calls.size(); ++s1) {
    const vector<Call> & left_calls{confident_calls[s1]};
    for (unsigned int s2{0}; s2 != confident_calls.size(); ++s2) {
      const vector<Call> & right_calls{confident_calls[s2]};
      int64_t opp{0};
      int64_t val{0};
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
      const double measure{1.0 * val / std::max(opp, static_cast<int64_t>(1))};
      measures[s1][s2] = measure;
      const double frac{(measure + 1) / 2};
      after_ps << "newpath " << s1 << " " << s2 << " gc moveto "
               << s1 + 1 << " " << s2 << " gc lineto "
               << s1 + 1 << " " << s2 + 1 << " gc lineto "
               << s1 << " " << s2 + 1 << " gc lineto closepath "
               << frac << " " << " 0 " << 1 - frac << " setrgbcolor "
               << "fill\n";
      if (s1 == s2) {
        uint64_t n_calls_correct{0};
        for (const Call & left : left_calls) {
          if (left.p_mom_a > 0.5) ++n_calls_correct;
        }
        self_graph->add_point(1.0 * n_calls_correct / left_calls.size(),
                              measure);
        if (measure > best_self_measure) {
          best_self_measure = measure;
          best_segment = s1;
        }
      }
    }
  }
  afterburner_graph->ps(after_ps.str());
  if (0) cout << best_segment;
#if 0
  // Is best segment itself flipped?
  uint64_t best_n_correct{0};
  const vector<Call> & left_calls{confident_calls[best_segment]};
  for (const Call & left : left_calls) {
    const bool mom_a1{left.p_mom_a > 0.5};
    best_n_correct += mom_a1;
  }
  const bool best_correct{best_n_correct * 2 > left_calls.size()};

  uint64_t n_fixed{0};
  uint64_t n_broken{0};
  uint64_t locus_correct{0};
  uint64_t segment_correct{0};
  uint64_t window_correct{0};
  // correct phase of all loci to have good measure with best segment
  // this includes best_segment loci
  for (unsigned int s2{0}; s2 != confident_calls.size(); ++s2) {
    const vector<Call> & right_calls{confident_calls[s2]};
    for (const Call & right : right_calls) {
      const uint64_t id2{right.id};
      const bool mom_a2{right.p_mom_a > 0.5};
      int64_t opp{0};
      int64_t val{0};
      for (uint64_t s{0}; s != hets.n_samples(); ++s) {
        for (const Call & left : left_calls) {
          const uint64_t id1{left.id};
          const uint64_t obs1{hets.unflipped(s, id1)};
          if (obs1 == 0 || obs1 == 3) continue;
          const bool mom_a1{left.p_mom_a > 0.5};
          const uint64_t obs2{hets.unflipped(s, id2)};
          if (obs2 == 0 || obs2 == 3) continue;
          ++opp;
          if ((obs1 == obs2 && mom_a1 == mom_a2) ||
              (obs1 != obs2 && mom_a1 != mom_a2)) {
            ++val;
          } else {
            --val;
          }
        }
      }
      const double measure{1.0 * val / std::max(opp, 1l)};
      if (measure < 0) {
        if (mom_a2 < 0.5) {
          ++n_fixed;
          ++locus_correct;
        } else {
          ++n_broken;
        }
      } else {
        if (mom_a2 >= 0.5) {
          ++locus_correct;
        }
      }
      fixed_p_mom_a_series->add_point(
          right.pos,
          (measure < 0 ? 1 - right.p_mom_a : right.p_mom_a) + gen());
      const bool segme_p_mom_a{
        measures[best_segment][s2] < 0 ? measure < 0 ? 1 - mom_a1 : mom_a1};
      window_p_mom_a_series->add_point(
          to_change.pos, window_p_mom_a + gen());
      window_correct += segment_p_mom_a > 0.5;
    }
  }

  if (!best_correct) {
    std::swap(n_fixed, n_broken);
    locus_correct = n_calls - locus_correct;
  }
  cout << "Fixed " << n_fixed << " broken " << n_broken << endl;
  cout << "locus fraction correct " << 1.0 * locus_correct / n_calls << endl;
  cout << "segment fraction correct " << 1.0 * segment_correct / n_calls
       << endl;
  cout << "window fraction correct " << 1.0 * window_correct / n_calls << endl;
#endif

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
