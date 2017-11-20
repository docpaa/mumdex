//
// cancer_pop
//
// Bridge counts for cancer and other population
// to look for differences, variable regions, etc
//
// Copyright 2017 Peter Andrews CSHL
//

#include <algorithm>
#include <deque>
#include <exception>
#include <iostream>
#include <queue>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "bridges.h"
#include "error.h"
#include "files.h"
#include "mumdex.h"
#include "population.h"
#include "repeats.h"

using std::cerr;
using std::cout;
using std::deque;
using std::endl;
using std::exception;
using std::lock_guard;
using std::lower_bound;
using std::move;
using std::mutex;
using std::ostringstream;
using std::priority_queue;
using std::string;
using std::vector;

using paa::BridgeInfo;
using paa::ChromosomeIndexLookup;
using paa::Error;
using paa::ExactRepeat;
using paa::ExactRepeats;
using paa::MaskerRepeat;
using paa::Mappability;
using paa::MergeHelper;
using paa::MergeHelperCompare;
using paa::Population;
using paa::Reference;
using paa::RepeatMasker;
using paa::Sample;

#if 1
using MUMdex = paa::FileMUMdex;
using Bridges = paa::FileBridges;
#else
using MUMdex = paa::MUMdex;
using Bridges = paa::Bridges;
#endif

using Anchors = MUMdex::Anchors;
using Bridge = Bridges::Bridge;

const bool verbose{false};

bool require_two_or_less{true};
bool require_one_family{true};
bool require_just_kids{true};

unsigned int rare_denovo_family_limit{10000};
const unsigned int min_support_length{20};
const unsigned int min_mate_support{5};
const unsigned int min_excess_mappability{2};

uint64_t n_no_cut{0};
uint64_t n_not_snp{0};
uint64_t n_two_or_less{0};
uint64_t n_one_family{0};
uint64_t n_children{0};
uint64_t n_children_one_family{0};
uint64_t n_no_rare_inherited{0};
uint64_t n_one_parent{0};
uint64_t n_big_count{0};
uint64_t n_good_support{0};
uint64_t n_good_mate_support{0};
uint64_t n_good_excess{0};
uint64_t n_good_pop{0};
uint64_t n_adjacent_cut{0};
uint64_t n_parent_coverage_cut{0};

int main(int argc, char* argv[])  try {
  if (--argc != 10) {
    throw Error("usage: cancer_pop ref family_file bridges_dir "
                "samples_dir chromosome start stop min_count masker exact");
  }

  // Command line arguments
  const string reference_file{argv[1]};
  const Reference ref{reference_file};
  const ChromosomeIndexLookup chr_lookup{ref};
  const Mappability map{reference_file, true};
  const string family_file{argv[2]};
  const Population pop{family_file};
  const string bridges_dir{argv[3]};
  const string samples_dir{argv[4]};
  const string chromosome_name{argv[5]};
  const unsigned int chromosome{chr_lookup[chromosome_name]};
  const unsigned int chromosome_offset{ref.offset(chromosome)};
  const unsigned int start{static_cast<unsigned int>(atoi(argv[6]))};
  const unsigned int stop{static_cast<unsigned int>(atoi(argv[7]))};
  const unsigned int min_count{static_cast<unsigned int>(atoi(argv[8]))};
  const RepeatMasker masker{argv[9], chr_lookup};
  const ExactRepeats exact{argv[10], chr_lookup};

  // Set up queue of bridge files
  deque<MergeHelper> helpers;
  using PQueue = priority_queue<MergeHelper*, vector<MergeHelper *>,
      MergeHelperCompare>;
  PQueue queue{MergeHelperCompare()};
  unsigned int n_samples_skipped{0};
  for (const auto s : pop.samples()) {
    ostringstream bridges_name;
    bridges_name << bridges_dir << "/"
                 << pop.sample(s) << "/chrbridges."
                 << chromosome << ".bin";
    MergeHelper helper{bridges_name.str(), s, start, stop};
    if (helper.valid_start()) {
      helpers.push_back(move(helper));
      queue.push(&helpers.back());
    } else {
      ++n_samples_skipped;
    }
    cerr << pop.family(pop.family(s)) << " "
         << pop.sample(s) << " "
         << pop.member(s) << " "
         << pop.is_cancer(s) << " "
         << pop.is_matched(s) << " "
         << pop.is_normal(s) << " "
         << pop.is_parent(s) << " "
         << bridges_name.str() << " "
         << queue.size() << endl;
  }

  if (n_samples_skipped)
    throw Error("Samples skipped") << n_samples_skipped;

  cout << "chrA\tposA\thighA\tchrB\tposB\thighB\tinvariant\toffset\tmapA\tmapB"
      "\tsA\tsB\tmsA\tmsB"
      "\tco1\tcoM\tcoH\tmo1\tmoM\tmoH\tm1\tmM\tmH\tp1\tpM\tpH"
      "\tmasker\tmasker_length\tmotif\tmotif_length\tmotif_copies";
  for (const auto s : pop.samples()) {
    if (pop.is_parent(s)) break;
    cout << "\t" << pop.sample(s);
  }
  cout << "\n";

  // Loop, adding one bridge at a time to all_bridges if the same event
  vector<BridgeInfo> all_bridges;
  vector<Sample> all_samples;
  bool reset{true};
  while (true) {
    // New bridge expected - start fresh
    if (reset) {
      all_bridges.clear();
      all_samples.clear();
      reset = false;
    }

    // Add a bridge to the list if the same event
    bool bridge_done{true};
    if (queue.size()) {
      MergeHelper * top{queue.top()};
      const BridgeInfo current{top->current()};
      if (all_bridges.empty() || !(all_bridges.back() < current)) {
        ++n_no_cut;
        bridge_done = false;
        all_bridges.push_back(current);
        const Sample sample{top->sample()};
        all_samples.emplace_back(sample);
        queue.pop();
        if (top->advance()) {
          queue.push(top);
        }
      }
      if (queue.empty()) bridge_done = true;
      if (bridge_done) reset = true;
    } else {
      if (all_samples.empty()) {
        if (verbose) cerr << "bridge queue empty" << endl;
        break;
      }
    }

    // No new bridges added, so bridge list is ready to process
    if (bridge_done) {
      if (false && !is_sorted(all_samples.begin(), all_samples.end()))
        throw Error("Samples unsorted");

      // The first bridge, the exemplar since all are the same
      const BridgeInfo & fbridge{all_bridges[0]};

      // SNP cut
      const bool is_snp{fbridge.chr1() == fbridge.chr2() &&
            fbridge.invariant() == 0 && fbridge.high1() != fbridge.high2()};
      if (is_snp) {
        continue;
      }

      // Inter-chromosomal cut
      if (fbridge.chr1() != fbridge.chr2()) continue;

      // Do bridge samples meet eligibility criteria for reporting?
      bool counts_good{false};
      for (unsigned int s{0}; s != all_samples.size(); ++s) {
        const Sample sample{all_samples[s]};
        // Ignore matched here
        if (pop.is_matched(sample)) continue;

        // See if matched for a cancer sample has counts and ignore if so
        if (pop.is_cancer(sample) && s + 1 != all_samples.size() &&
            pop.family(sample) == pop.family(all_samples[s + 1])) continue;

        const BridgeInfo & bridge{all_bridges[s]};
        if (bridge.bridge_count() < min_count) continue;
        counts_good = true;
        break;
      }

      if (!counts_good) continue;

      // Gather max support lengths
      unsigned int max_anchor1_support{0};
      unsigned int max_anchor2_support{0};
      unsigned int max_mate_anchor1_support{0};
      unsigned int max_mate_anchor2_support{0};
      for (const BridgeInfo & sample_bridge : all_bridges) {
        if (max_anchor1_support < sample_bridge.anchor1_length()) {
          max_anchor1_support = sample_bridge.anchor1_length();
        }
        if (max_anchor2_support < sample_bridge.anchor2_length()) {
          max_anchor2_support = sample_bridge.anchor2_length();
        }
        if (max_mate_anchor1_support < sample_bridge.mate_anchor1_length()) {
          max_mate_anchor1_support = sample_bridge.mate_anchor1_length();
        }
        if (max_mate_anchor2_support < sample_bridge.mate_anchor2_length()) {
          max_mate_anchor2_support = sample_bridge.mate_anchor2_length();
        }
      }

      // Max support cuts
      // Short anchor support cut
      if (max_anchor1_support < min_support_length ||
          max_anchor2_support < min_support_length) continue;
      ++n_good_support;

      // Short mate support cut
      if (max_mate_anchor1_support < min_mate_support ||
          max_mate_anchor2_support < min_mate_support) continue;
      ++n_good_mate_support;

      // Get mappability information
      const unsigned int mum1_abspos{chromosome_offset + fbridge.pos1()};
      const unsigned int mum2_abspos{ref.offset(
          fbridge.chr2()) + fbridge.pos2()};
      const unsigned int mum1_map{map.low_high(fbridge.high1(), mum1_abspos)};
      const unsigned int mum2_map{map.low_high(fbridge.high2(), mum2_abspos)};

      // Short excess mappability cut
      if (max_anchor1_support < min_excess_mappability + mum1_map) {
        continue;
      } else {
        if (max_anchor2_support < min_excess_mappability + mum2_map) {
          continue;
        }
      }
      ++n_good_excess;

      // Gather counts for bridge
      // cancer only / matched only / in matched / pop x 1 / med / high
      static vector<vector<unsigned int>> counts(4, vector<unsigned int>(3));
      for (vector<unsigned int> & vec : counts) {
        for (unsigned int & val : vec) {
          val = 0;
        }
      }
      for (unsigned int s{0}; s != all_samples.size(); ++s) {
        const Sample sample{all_samples[s]};
        // Ignore cancer + matched here
        if (pop.is_cancer(sample) &&
            (s + 1 != all_samples.size() &&
             pop.family(sample) == pop.family(all_samples[s + 1]))) continue;

        const BridgeInfo & bridge{all_bridges[s]};
        const unsigned int count_index{bridge.bridge_count() == 1 ?
              0U : (bridge.bridge_count() >= min_count ? 2U : 1U)};

        const unsigned int type_index{pop.is_parent(sample) ? 3U :
              (pop.is_cancer(sample) ? 0U :
               (pop.is_matched(sample) ? 2U : 4U))};
        if (type_index == 4) throw Error("Bad type index");

        ++counts[type_index][count_index];
        // Seen in matched but not cancer
        if (pop.is_matched(sample) &&
            (s == 0 || pop.family(sample) != pop.family(all_samples[s - 1])))
          ++counts[1][count_index];
      }

      // Output count information
      cout << ref.name(fbridge.chr1()) << "\t"
           << fbridge.pos1() << "\t"
           << fbridge.high1() << "\t"
           << ref.name(fbridge.chr2()) << "\t"
           << fbridge.pos2() << "\t"
           << fbridge.high2() << "\t"
           << fbridge.invariant() << "\t"
           << fbridge.offset() << "\t"
           << mum1_map << "\t"
           << mum2_map << "\t"
           << max_anchor1_support << "\t"
           << max_anchor2_support << "\t"
           << max_mate_anchor1_support << "\t"
           << max_mate_anchor2_support;

      // Summary counts
      for (const vector<unsigned int> & vec : counts) {
        for (const unsigned int val : vec) {
          cout << "\t" << val;
        }
      }

      // Repeat masker annotation
      cout << "\t";
      const vector<const MaskerRepeat *> masker_result{masker(
          fbridge.chr1(), fbridge.pos1(), fbridge.chr2(), fbridge.pos2())};
      if (masker_result.size()) {
        const MaskerRepeat & repeat{*masker_result.front()};
        cout << repeat.repeat << "-" << repeat.family << "\t"
             << repeat.stop - repeat.start;
      } else {
        cout << "-\t0";
      }

      // Exact repeats annotation
      cout << "\t";
      const vector<const ExactRepeat *> exact_result{exact(
          fbridge.chr1(), fbridge.pos1(), fbridge.chr2(), fbridge.pos2())};
      if (exact_result.size()) {
        cout << exact_result.front()->motif << "\t"
             << exact_result.front()->motif.size() << "\t"
             << exact_result.front()->n_copies;
      } else {
        cout << "-\t0\t0";
      }

      // Individual counts
      unsigned int csi{0};
      for (const auto s : pop.samples()) {
        if (pop.is_parent(s)) break;
        if (csi != all_samples.size() && all_samples[csi] == s) {
          cout << "\t" << all_bridges[csi].bridge_count();
          ++csi;
        } else {
          cout << "\t" << 0;
        }
      }
      cout << "\n";
    }
  }

  cerr << "finished bridge loop" << endl;

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

