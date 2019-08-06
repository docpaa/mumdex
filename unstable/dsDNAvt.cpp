//
// dsDNAvt.cpp
//
// Check properties of a special sequencing run
// smash library with two varietal tags per pair
// look for template jumping
//
// Copyright 2017 Peter Andrews @ CSHL
//

#include <algorithm>
#include <exception>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "error.h"
#include "mumdex.h"

using std::cout;
using std::cerr;
using std::endl;
using std::exception;
using std::ifstream;
using std::istringstream;
using std::make_pair;
using std::map;
using std::min;
using std::pair;
using std::set;
using std::string;
using std::vector;

using paa::reverse_complement;
using paa::ChromosomeIndexLookup;
using paa::Error;
using paa::Mappability;
using paa::Reference;

using VTS = pair<string, string>;

class MapInfo {
 public:
  string chr{};
  unsigned int pos{};
  unsigned int abspos{};
  bool read_2{};
  unsigned int offset{};
  unsigned int length{};
  bool flipped{};
};

class PairInfo {
 public:
  string read_name{};
  string tag[2]{"", ""};
  string read[2]{"", ""};
  vector<MapInfo> maps{};
  bool pair_flipped{false};
};

unsigned int hamming_dist(const string & lhs, const string & rhs) {
  // if (lhs.size() != rhs.size()) throw Error("Hamming size mismatch");
  const unsigned int msize{static_cast<unsigned int>(
      std::min(lhs.size(), rhs.size()))};
  const unsigned int sdiff{static_cast<unsigned int>(
      lhs.size() > rhs.size() ? lhs.size() - rhs.size() :
      rhs.size() - lhs.size())};
  unsigned int result{sdiff};
  for (unsigned int b{0}; b != msize; ++b) {
    if (lhs[b] != rhs[b]) ++result;
  }
  return result;
}

int main(int argc, char * argv[]) try {
  --argc;
  if (argc != 2) throw Error("usage: dsDNAvt ref maps_file");

  const Reference ref{argv[1]};
  const ChromosomeIndexLookup lookup{ref};
  const Mappability mappability{ref};

  ifstream maps_file{argv[2]};
  if (!maps_file) throw Error("Problem opening maps file") << argv[2];

  string read_name;
  string vt1;
  string vt2;
  unsigned int vtm1;
  unsigned int vtm2;
  string read_name2{""};
  string read;
  string read1;
  string read2;
  string read_name3{""};
  unsigned int read_2;
  string chromosome;
  unsigned int pos;
  unsigned int offset;
  unsigned int length;
  char orient;


  vector<PairInfo> pairs;
  map<string, PairInfo> mpairs;
  while (maps_file >> read_name >> read_2 >> vtm1 >> vtm2
         >> vt1 >> vt2 >> read1 >> read2 >> chromosome >> pos >> offset
         >> length >> orient) {
    if (!maps_file) continue;
    if (vtm1 || vtm2) continue;
    if (read1.size() < 50 || read2.size() < 50) continue;
    if (length < 20) continue;
    const unsigned int abspos{ref.abspos(lookup[chromosome], pos)};
    const unsigned int map_len{mappability.low(abspos)};
    if (length < map_len + 4) continue;

    PairInfo info;
    info.read_name = read_name;
    info.tag[0] = vt1;
    info.tag[1] = vt2;
    if (vt2 < vt1) {
      using std::swap;
      swap(vt1, vt2);
      const string rcr2{reverse_complement(read2)};
      read2 = reverse_complement(read1);
      read1 = rcr2;
      info.pair_flipped = true;
    }
    info.read[0] = read1;
    info.read[1] = read2;

    auto pi = mpairs.insert(make_pair(read_name, info));

    MapInfo minfo;
    minfo.chr = chromosome;
    minfo.pos = pos;
    minfo.abspos = abspos;
    minfo.read_2 = read_2;
    minfo.offset = offset;
    minfo.length = length;
    minfo.flipped = orient == '-';
    pi.first->second.maps.push_back(minfo);
  }

  for (auto & pair : mpairs) {
    sort(pair.second.maps.begin(), pair.second.maps.end(),
         [](const MapInfo & lhs, const MapInfo & rhs) {
           if (lhs.read_2 == rhs.read_2) {
             if (lhs.read_2) {
               return lhs.offset > rhs.offset;
             } else {
               return lhs.offset < rhs.offset;
             }
           } else {
             return lhs.read_2 < rhs.read_2;
           }
         });
    pairs.push_back(pair.second);
  }

  cerr << "Loaded " << pairs.size() << " pairs of complete info" << endl;

  // Sort by vts, look for different sequences
  sort(pairs.begin(), pairs.end(),
       [](const PairInfo & lhs, const PairInfo & rhs) {
         if (lhs.tag[0] == rhs.tag[0]) {
           return lhs.tag[1] < rhs.tag[1];
         } else {
           return lhs.tag[0] < rhs.tag[0];
         }
    });

  unsigned int n_identical_vts{1};
  set<string> sequences;
  vector<unsigned int> hamming_counts(150);
  map<pair<unsigned int, unsigned int>, unsigned int> vt_counts;
  for (unsigned int p{0}; p != pairs.size(); ++p) {
    const PairInfo & info{pairs[p]};
    if (p) {
      const PairInfo & last{pairs[p - 1]};
      if (info.tag[0] == last.tag[0] &&
          info.tag[1] == last.tag[1]) {
        ++n_identical_vts;
      } else {
        pair<unsigned int, unsigned int> key{n_identical_vts, sequences.size()};
        if (n_identical_vts == 2 && sequences.size() == 2) {
          const string & seq1{*sequences.begin()};
          const string & seq2{*++sequences.begin()};
          unsigned int hamming{0};
          for (unsigned int s{0}; s != seq1.size() && s != seq2.size(); ++s) {
            if (seq1[s] != seq2[s]) {
              ++hamming;
            }
          }
          ++hamming_counts[hamming];
          if (hamming == 39 || hamming == 40) {
            cout << "Read mismatch sample " << info.pair_flipped << endl;
            bool seen_read_2{false};
            for (const auto minfo : pairs[p - 1].maps) {
              const unsigned int abspos{minfo.abspos};
              if (!seen_read_2 && minfo.read_2) {
                seen_read_2 = true;
                cout << " |";
              }
              cout << " " << abspos;
              bool found{false};
              for (const auto & minfo2 : pairs[p - 2].maps) {
                const unsigned int abspos2{minfo2.abspos};
                const unsigned int apdiff{abspos2 > abspos ? abspos2 - abspos :
                      abspos - abspos2};
                if (apdiff < 10) found = true;
              }
              if (found) {
                cout << "*(" << minfo.offset - 10 << ")";
              }
            }
            cout << endl;
            seen_read_2 = false;
            for (const auto minfo : pairs[p - 2].maps) {
              if (!seen_read_2 && minfo.read_2) {
                seen_read_2 = true;
                cout << " |";
              }
              cout << " " << minfo.abspos;
            }
            cout << endl;
            cout << pairs[p - 1].read[0] << endl;
            cout << pairs[p - 2].read[0] << " "
                 << hamming_dist(pairs[p - 1].read[0], pairs[p - 2].read[0])
                 << endl;
            cout << pairs[p - 1].read[1] << endl;
            cout << pairs[p - 2].read[1] << " "
                 << hamming_dist(pairs[p - 1].read[1], pairs[p - 2].read[1])
                 << endl;
            cout << endl;
          }
        }
        ++vt_counts[key];
        sequences.clear();
        n_identical_vts = 1;
      }
    }
    const string sequence{info.read[0] + info.read[1]};
    sequences.insert(sequence);
  }

  cout << "Varietal tag counts" << endl;
  for (const auto elem : vt_counts) {
    cout << elem.first.first << " "
         << elem.first.second << " "
         << elem.second << endl;
  }

  cout << "Hamming counts for 2 2" << endl;
  for (unsigned int h{0}; h != hamming_counts.size(); ++h) {
    cout << h << " " << hamming_counts[h] << endl;
  }

  // Sort by sequences, look for different vts
  sort(pairs.begin(), pairs.end(),
       [](const PairInfo & lhs, const PairInfo & rhs) {
         if (lhs.read[0] == rhs.read[0]) {
           return lhs.read[1] < rhs.read[1];
         } else {
           return lhs.read[0] < rhs.read[0];
         }
    });

  unsigned int n_identical_reads{1};
  set<string> vt_sequences;
  map<pair<unsigned int, unsigned int>, unsigned int> seq_counts;
  for (unsigned int p{0}; p != pairs.size(); ++p) {
    const PairInfo & info{pairs[p]};
    if (p) {
      const PairInfo & last{pairs[p - 1]};
      if (info.read[0] == last.read[0] &&
          info.read[1] == last.read[1]) {
        ++n_identical_reads;
      } else {
        pair<unsigned int, unsigned int> key{
          n_identical_reads, vt_sequences.size()};
        ++seq_counts[key];
        vt_sequences.clear();
        n_identical_reads = 1;
      }
    }
    const string sequence{info.read[0] + info.read[1]};
    vt_sequences.insert(sequence);
  }

  cout << "Sequence counts" << endl;
  for (const auto elem : seq_counts) {
    cout << elem.first.first << " "
         << elem.first.second << " "
         << elem.second << endl;
  }

  return 0;
} catch (Error & e) {
  cerr << "paa::Error:" << endl;
  cerr << e.what() << endl;
  return 1;
}
catch (exception & e) {
  cerr << "std::exception" << endl;
  cerr << e.what() << endl;
  return 1;
} catch (...) {
  cerr << "unknown exception was caught" << endl;
  return 1;
}
