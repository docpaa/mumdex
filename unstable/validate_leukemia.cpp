//
// validate_mumdex.cpp
//
// Copyright 2017 Peter Andrews @ CSHL
//

#include <algorithm>
#include <iostream>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "bridges.h"
#include "error.h"
#include "mumdex.h"
#include "strings.h"
#include "utility.h"

using std::cout;
using std::cerr;
using std::endl;
using std::ifstream;
using std::istringstream;
using std::map;
using std::max;
using std::string;
using std::vector;

using paa::replace;
using paa::sout;
using paa::Bridge;
using paa::BridgeCounts;
using paa::Bridges;
using paa::ChromosomeIndexLookup;
using paa::Error;
using paa::Mappability;
using paa::MUMdex;
using paa::Reference;

void parse_counts(const string & in_string, unsigned int vals[4]) {
  istringstream in{in_string.c_str()};
  for (const unsigned int i : {0, 1, 2, 3}) {
    in >> vals[i];
    in.get();
  }
}

class CandInfo {
 public:
  class ParseError {};
  explicit CandInfo(const string & line_) {
    line = line_;
    istringstream input{line.c_str()};
    input >> event_class
          >> chrA >> posA >> highA >> chrB >> posB >> highB
          >> invariant >> offset >> ori
          >> supA >> supB >> mCA >> mCB >> mSA >> mSB
          >> count_141 >> count_149 >> count_267
          >> p_value >> log_p_value >> refbase >> altbase;
    replace(event_class, '0', 'L');
    if (!input) {
      throw ParseError{};
    }
  }
  string line{};
  string event_class{};

  string chrA{};
  unsigned int posA{};
  bool highA{};
  string chrB{};
  unsigned int posB{};
  bool highB{};
  int invariant{};
  int offset{};
  char ori{};

  unsigned int supA{};
  unsigned int supB{};
  unsigned int mCA{};
  unsigned int mSA{};
  unsigned int mCB{};
  unsigned int mSB{};

  unsigned int count_141{};
  unsigned int count_149{};
  unsigned int count_267{};

  double p_value{};
  double log_p_value{};
  string refbase{};
  string altbase{};
};

int main(int argc, char ** argv)  try {
  // Check arguments
  if (--argc < 2)
    throw Error("usage: validate_mumdex cand mumdexes...");

  // Input file
  const string cand_file_name{argv[1]};
  ifstream cand_file{cand_file_name.c_str()};
  if (!cand_file) throw Error("Could not open file") << cand_file_name;

  // Read Input file into candidate objects
  string line;
  vector<CandInfo> info;
  while (cand_file) {
    getline(cand_file, line);
    if (line.size()) {
      try {
        CandInfo i{line};
        info.push_back(i);
      } catch (CandInfo::ParseError & e) {
        if (line[0] != 'r' || line[1] != 'a') {
          cerr << "parse error for line: " << line << endl;
        }
      }
    }
  }

  argc -= 1;
  argv += 1;

  // Load mumdexes for family members
  vector<MUMdex> mumdexes;
  mumdexes.reserve(argc);
  while (argc) {
    mumdexes.emplace_back(argv[1]);
    --argc;
    ++argv;
  }
  const Reference & ref{mumdexes.front().reference()};
  const ChromosomeIndexLookup lookup{ref};
  const Mappability mappability{ref};

  const unsigned int early{152};

  // Header line
  cout << "chrA\tposA\thighA\tchrB\tposB\thighB\tinv\toff";
  cout << "\t141\t149\t267";
  const vector<string> counts{"cc", "ac", "bc", "br"};
  const string anchor{"AB"};
  const vector<string> members{"141", "149", "267", "268"};
  if (0) {
    for (const bool ab : {false, true}) {
      for (const unsigned int c : {0, 1, 2, 3}) {
        for (unsigned int m{0}; m != mumdexes.size(); ++m) {
          if (ab || c < 2) {
            cout << "\t" << members[m] << counts[c];
            if (c < 2) cout << anchor[ab];
          }
        }
      }
    }
  }
  cout << "\tref\talt";
  for (unsigned int m{0}; m != mumdexes.size(); ++m) {
    for (const string & base : {"A", "C", "G", "T"}) {
      cout << "\t" << members[m] << base;
    }
  }
  for (unsigned int m{0}; m != mumdexes.size(); ++m) {
    cout << "\t" << members[m] << "ar";
  }
  cout << "\tclass\toclass\tsuccess\tratio" << endl;

  unsigned int n_success{0};

  for (const CandInfo & i : info) {
    unsigned int max_base_seen{0};
    const string chrs[2]{i.chrA, i.chrB};
    const unsigned int chr[2]{lookup[i.chrA], lookup[i.chrB]};
    const unsigned int pos[2]{i.posA, i.posB};
    const unsigned int high[2]{i.highA, i.highB};
    const int npos[2]{static_cast<int>(pos[0]) + (high[0] ? 1 : -1),
          static_cast<int>(pos[1]) + (high[1] ? 1 : -1)};
    const unsigned int abspos[2]{ref.abspos(chr[0], pos[0]),
          ref.abspos(chr[1], pos[1])};
    const unsigned int map[2]{mappability.low_high(high[0], abspos[0]),
          mappability.low_high(high[1], abspos[1])};

    const vector<BridgeCounts> bridge_counts2{
      [&mumdexes, &ref, &mappability, &i, chr, pos, high]() {
        vector<BridgeCounts> result;
        for (unsigned int m{0}; m != mumdexes.size(); ++m) {
          const MUMdex & mumdex{mumdexes[m]};
          result.emplace_back(mumdex, ref, mappability,
                              chr, pos, high, i.offset, true);
        }
        return result;
      }()};

    for (const bool b : {false, true}) {
      cout << chrs[b] << "\t" << pos[b] << "\t" << high[b] << "\t";
    }
    cout << i.invariant << "\t" << i.offset;
    cout << "\t" << i.count_141 << "\t" << i.count_149 << "\t" << i.count_267;

    vector<unsigned int> bridge_counts(mumdexes.size());
    vector<unsigned int> ref_counts(mumdexes.size());

    // How many mums contain the anchor base
    for (const bool b : {false, true}) {
      // Find counts
      vector<unsigned int> base_counts(mumdexes.size());
      vector<unsigned int> anchor_counts(mumdexes.size());
      for (unsigned int m{0}; m != mumdexes.size(); ++m) {
        const MUMdex & mumdex{mumdexes[m]};

        // Count bridges
        if (b) {
          for (const Bridge & bridge : Bridges{
              chr[0], pos[0], static_cast<bool>(high[0]),
                  chr[1], pos[1], static_cast<bool>(high[1]),
                  i.invariant, mumdex}) {
            if (0) std::cerr << bridge.pair_index();
            ++bridge_counts[m];
          }
        }

        const auto begin_index = mumdex.lower_bound(
            chr[b], pos[b] > early ? pos[b] - early : 0);
        const auto end_index = mumdex.lower_bound(
            chr[b], pos[b] + 1);
        for (auto index = begin_index; index != end_index; ++index) {
          // const auto pair = mumdex.pair(*index);
          const auto mum = mumdex.mum(*index);
          const auto pair = mumdex.pair(*index);
          const int read_length(pair.length(mum.read_2()));

          if (pos[b] >= mum.position0() &&
              pos[b] < mum.position0() + mum.length()) {
            ++base_counts[m];
            if ((pos[b] == mum.position0() && !high[b] &&
                 (mum.flipped() ? !mum.touches_end() : mum.offset())) ||
                (pos[b] == mum.position0() + mum.length() - 1 && high[b] &&
                 (mum.flipped() ? mum.offset() : !mum.touches_end()))) {
              ++anchor_counts[m];
            }
            // Bridge limits in mum space
            const int bp[2]{static_cast<int>(pos[b]) + (high[b] ? -1 : 1) *
                  static_cast<int>(map[b] - 1),
                  static_cast<int>(pos[b]) + (high[b] ? 1 : -1) *
                  static_cast<int>(i.offset + map[1 - b] - 1)};
            bool good{true};
            for (const bool e : {false, true}) {
              const int p{mum.position0_to_read(bp[e])};
              if (p < 0 || p >= read_length) {
                good = false;
              }
            }
            if (good) {
              // Is reference continuous at anchor position?
              if (npos[b] >= static_cast<int>(mum.position0()) &&
                  npos[b] < static_cast<int>(mum.position0() + mum.length())) {
                ++ref_counts[m];
              }
            }
          }
        }
      }

      if (0) {
        // Output counts
        for (unsigned int m{0}; m != mumdexes.size(); ++m) {
          max_base_seen = max(base_counts[m], max_base_seen);
          cout << "\t" << base_counts[m];
        }
        for (unsigned int m{0}; m != mumdexes.size(); ++m) {
          cout << "\t" << anchor_counts[m];
        }
      }
    }
    if (0) {
      for (unsigned int m{0}; m != mumdexes.size(); ++m) {
        cout << "\t" << bridge_counts[m];
      }
      for (unsigned int m{0}; m != mumdexes.size(); ++m) {
        cout << "\t"
             << bridge_counts[m] / (bridge_counts[m] + ref_counts[m] / 2.0);
      }
    }
    cout << "\t" << i.refbase << "\t" << i.altbase;
    for (unsigned int m{0}; m != mumdexes.size(); ++m) {
      for (const string & base : {"A", "C", "G", "T"}) {
        const unsigned int count{base[0] == ref[chr[0]][i.posB - 1] ?
              bridge_counts2[m].bridge.ref_count :
              bridge_counts2[m].bridge.get_allele_count(base)};
        cout << "\t" << count;
      }
    }
    for (unsigned int m{0}; m != mumdexes.size(); ++m) {
      const unsigned int remission_count{
        bridge_counts2[m].bridge.get_allele_count(i.altbase)};
      cout << "\t" << 1.0 * remission_count /
          (remission_count + bridge_counts2[m].bridge.ref_count);
    }

    // Success or failure?
    const unsigned int min_count{max(1U, *min_element(bridge_counts.begin(),
                                                      bridge_counts.end()))};
    const unsigned int max_count{*max_element(bridge_counts.begin(),
                                              bridge_counts.end())};
    const double ratio{1.0 * max_count / min_count};
    string event_class;
    if (max_count > 50) {
      for (unsigned int s{0}; s != 3; ++s) {
        if (10 * bridge_counts[s] > max_count) {
          event_class += "H";
        } else if (50 * bridge_counts[s] < max_count) {
          event_class += "L";
        } else {
          event_class += "M";
        }
      }
    } else {
      event_class = "XXX";
    }

    const bool success{event_class == i.event_class};
    if (success) ++n_success;
    cout << "\t" << event_class << "\t" << i.event_class
         << "\t" << success << "\t" << ratio;
    cout << endl;
  }

  cerr << n_success << " of " << info.size() << " succeeded, or "
       << 100.0 * n_success / info.size() << " percent" << endl;

  return 0;
} catch (Error & e) {
  cerr << "paa::Error:" << endl;
  cerr << e.what() << endl;
  return 1;
} catch(std::exception & e) {
  cerr << e.what() << '\n';
  return 1;
} catch(...) {
  cerr << "Some exception was caught.\n";
  return 1;
}
