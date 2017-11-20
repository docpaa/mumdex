//
// compare_yoonha.cpp
//
// Copyright 2016 Peter Andrews @ CSHL
//

#include <algorithm>
#include <array>
#include <iostream>
#include <fstream>
#include <functional>
#include <map>
#include <memory>
#include <random>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "error.h"
#include "mumdex.h"
#include "psplot.h"
#include "utility.h"

using std::accumulate;
using std::array;
using std::cout;
using std::cerr;
using std::endl;
using std::ifstream;
using std::istringstream;
using std::function;
using std::map;
using std::max;
using std::min;
using std::mt19937_64;
using std::pair;
using std::random_device;
using std::sort;
using std::string;
using std::uniform_real_distribution;
using std::vector;

using paa::serr;
using paa::sout;
using paa::nunset;
using paa::unset;
using paa::Bounds;
using paa::ChromosomeIndexLookup;
using paa::Error;
using paa::Marker;
using paa::PosInfo;
using paa::PSDoc;
using paa::PSPage;
using paa::PSGraph;
using paa::PSSeries;
using paa::PSXYSeries;
using paa::PSXYMSeries;
using paa::Reference;

void parse_counts(const string & in_string, unsigned int vals[4]) {
  istringstream in{in_string.c_str()};
  for (const unsigned int i : {0, 1, 2, 3}) {
    in >> vals[i];
    if (!in) throw Error("Parse error in parse_counts");
    in.get();
  }
}

const int min_inv{-47};
const int max_inv{32};

class CandInfo {
 public:
  class ParseError {};
  explicit CandInfo(const string & line) {
    istringstream input{line.c_str()};
    string bridges_;
    string anchorsA_;
    string anchorsB_;
    string coveragesA_;
    string coveragesB_;
    input >> family >> parent >> kids >> nFam >> nSam >> nInFam
          >> chr >> chrA >> posA >> highA >> chrB >> posB >> highB
          >> type >> ori >> invariant >> offset >> bridge
          >> supA >> supB >> mCA >> mSA >> mCB >> mSB >> mapA >> mapB
          >> minParCov >> bridges_ >> anchorsA_ >> anchorsB_
          >> coveragesA_ >> coveragesB_
          >> emA >> b0A >> b1A >> b2A >> b3A
          >> emB >> b0B >> b1B >> b2B >> b3B
          >> rlen >> mn >> mlen >> motif;

    if (!input) {
      throw ParseError{};
    }
    parse_counts(bridges_, bridges);
    parse_counts(anchorsA_, anchorsA);
    parse_counts(anchorsB_, anchorsB);
    parse_counts(coveragesA_, coveragesA);
    parse_counts(coveragesB_, coveragesB);

    if (bridge >= 5 &&
        supA >= 25 && supB >= 25 &&
        mSA >= 20 && mSB >= 20 &&
        minParCov >= 10) {
      strong = true;
      ++n_strong;
      if (chrA == chrB && highA != highB &&
          invariant >= -47 && invariant <= 32) {
        strong_indel = true;
        ++n_strong_indel;
      }
    }
    if (chrA == chrB && highA != highB &&
        invariant >= min_inv && invariant <= max_inv) {
      indel = true;
      ++n_indel;
    }
  }
  string family{};
  string parent{};
  string kids{};
  unsigned int nFam{};
  unsigned int nSam{};
  unsigned int nInFam{};
  unsigned int chr{};
  string chrA{};
  unsigned int posA{};
  bool highA{};
  string chrB{};
  unsigned int posB{};
  bool highB{};
  string type{};
  char ori{};
  int invariant{};
  int offset{};
  unsigned int bridge{};
  unsigned int supA{};
  unsigned int supB{};
  unsigned int mCA{};
  unsigned int mSA{};
  unsigned int mCB{};
  unsigned int mSB{};
  unsigned int mapA{};
  unsigned int mapB{};
  unsigned int minParCov{};
  unsigned int bridges[4]{0, 0, 0, 0};
  unsigned int anchorsA[4]{0, 0, 0, 0};
  unsigned int anchorsB[4]{0, 0, 0, 0};
  unsigned int coveragesA[4]{0, 0, 0, 0};
  unsigned int coveragesB[4]{0, 0, 0, 0};
  unsigned int emA{};
  unsigned int b0A{};
  unsigned int b1A{};
  unsigned int b2A{};
  unsigned int b3A{};
  unsigned int emB{};
  unsigned int b0B{};
  unsigned int b1B{};
  unsigned int b2B{};
  unsigned int b3B{};
  unsigned int rlen{};
  unsigned int mn{};
  unsigned int mlen{};
  string motif{};

  bool strong{false};
  bool strong_indel{false};
  bool indel{false};
  static unsigned int n_strong;
  static unsigned int n_strong_indel;
  static unsigned int n_indel;
  bool operator<(const CandInfo & rhs) const {
    if (chr == rhs.chr) {
      if (posA == rhs.posA) {
        return invariant < rhs.invariant;
      } else {
        return posA < rhs.posA;
      }
    } else {
      return chr < rhs.chr;
    }
  }
  bool operator<(const PosInfo & rhs) const {
    if (chr == rhs.chr) {
      return posA < rhs.pos;
    } else {
      return chr < rhs.chr;
    }
  }
};

unsigned int CandInfo::n_strong{0};
unsigned int CandInfo::n_strong_indel{0};
unsigned int CandInfo::n_indel{0};

class YoonhaInfo {
 public:
  class ParseError {};
  explicit YoonhaInfo(const string & line,
                      const ChromosomeIndexLookup & lookup) {
    istringstream input{line.c_str()};
    string dummy;
    input >> family >> location >> variant >> property;
    input.get();
    getline(input, Count, '\t');
    input >> bestStateScore >> denovoScore >> chi2Score;
    input.get();
    getline(input, chi2expectedCount, '\t');
    getline(input, bestState, '\t');
    getline(input, alpha, '\t');
    input >> dummy;
    input >> dummy;
    input.get();
    getline(input, fromParent, '\t');
    input >> effectType >> effectGene >> effectDetails;
    input >> numFamIndel >> numParentsIndel;
    input >> totalParentsCount >> totalFamCounts;
    input >> yxCall;
    input.get();
    getline(input, yxCount, '\t');
    input >> yuCall;
    input.get();
    getline(input, yuCount, '\t');
    input >> dummy;
    input >> dummy;
    input >> dummy;
    input >> dummy;
    input >> dummy;
    input >> strength;
    input.get();

    if (!input) {
      throw ParseError{};

      serr << "f" << family << '\n';
      serr << "l" << location << '\n';
      serr << "v" << variant << '\n';
      serr << "p" << property << '\n';
      serr << "c" << Count << '\n';
      serr << "b" << bestStateScore << '\n';
      serr << "d" << denovoScore << '\n';
      serr << "c" << chi2Score << '\n';
      serr << "c" << chi2expectedCount << '\n';
      serr << "b" << bestState << '\n';
      serr << "a" << alpha << '\n';
      serr << "f" << fromParent << '\n';
      serr << "e" << effectType << '\n';
      serr << "e" << effectGene << '\n';
      serr << "e" << effectDetails << '\n';
      serr << "n" << numFamIndel << '\n';
      serr << "n" << numParentsIndel << '\n';
      serr << "t" << totalParentsCount << '\n';
      serr << "t" << totalFamCounts << '\n';
      serr << "y" << yxCall << '\n';
      serr << "y" << yxCount << '\n';
      serr << "y" << yuCall << '\n';
      serr << "y" << yuCount << '\n';
      serr << "s" << strength << '\n';
      serr << "s" << SSC_freq << '\n';
      serr << "e" << EVS_freq << '\n';
      serr << "e" << E65_freq << '\n';
      // if (family != "familyId") exit(1);
    }
    getline(input, dummy, '\t');
    getline(input, SSC_freq, '\t');
    getline(input, EVS_freq, '\t');
    getline(input, E65_freq);

    istringstream pos_stream{location.c_str()};
    getline(pos_stream, dummy, ':');
    chr = lookup[dummy];
    pos_stream >> pos;

    istringstream inv_stream{variant.c_str()};
    string indel;
    getline(inv_stream, indel, '(');
    if (indel == "ins") {
      getline(inv_stream, indel, ')');
      invariant = static_cast<int>(indel.size());
    } else if (indel == "del") {
      inv_stream >> invariant;
      invariant = -invariant;
    } else {
      throw Error("variant is not ins or del");
    }

    if (strength == "S") {
      strong = true;
      ++n_strong;
    }
  }

  unsigned int chr{};
  unsigned int pos{};
  int invariant{};

  string family{};
  string location{};
  string variant{};
  string property{};
  string Count{};
  double bestStateScore{};
  double denovoScore{};
  double chi2Score{};
  string chi2expectedCount{};
  string bestState{};
  string alpha{};
  string fromParent{};
  string effectType{};
  string effectGene{};
  string effectDetails{};
  unsigned int numFamIndel{};
  unsigned int numParentsIndel{};
  unsigned int totalParentsCount{};
  unsigned int totalFamCounts{};
  double yxCall{};
  string yxCount{};
  string yuCall{};
  string yuCount{};
  string strength{};
  string SSC_freq{};
  string EVS_freq{};
  string E65_freq{};
  bool strong{false};
  static unsigned int n_strong;

  bool operator<(const YoonhaInfo & rhs) const {
    if (chr == rhs.chr) {
      if (pos == rhs.pos) {
        return invariant < rhs.invariant;
      } else {
        return pos < rhs.pos;
      }
    } else {
      return chr < rhs.chr;
    }
  }
  bool operator<(const PosInfo & rhs) const {
    if (chr == rhs.chr) {
      return pos < rhs.pos;
    } else {
      return chr < rhs.chr;
    }
  }
};

unsigned int YoonhaInfo::n_strong{0};

const unsigned int default_distance{300000000};
unsigned int distance(const CandInfo & mumdex, const YoonhaInfo & yoonha) {
  if (mumdex.family != yoonha.family) throw Error("Distance family");
  if (mumdex.posA > mumdex.posB) throw Error("anchors unordered");
  if (mumdex.chr != yoonha.chr) return default_distance;
  if (yoonha.pos >= mumdex.posA && yoonha.pos <= mumdex.posB) return 0;
  if (yoonha.pos < mumdex.posA) return mumdex.posA - yoonha.pos;
  return yoonha.pos - mumdex.posB;
}

struct CandMatch {
  CandMatch(const CandInfo * m, const YoonhaInfo * y, const unsigned int d) :
      mumdex{m}, yoonha{y}, distance{d} {}
  const CandInfo * mumdex;
  const YoonhaInfo * yoonha;
  unsigned int distance;
  bool operator<(const CandMatch & rhs) const {
    if (distance == rhs.distance) {
      return *mumdex < *rhs.mumdex;
    } else {
      return distance < rhs.distance;
    }
  }
};

struct YoonhaMatch {
  YoonhaMatch(const YoonhaInfo * y, const CandInfo * m, const unsigned int d) :
      yoonha{y}, mumdex{m}, distance{d} {}
  const YoonhaInfo * yoonha;
  const CandInfo * mumdex;
  unsigned int distance;
  bool operator<(const YoonhaMatch & rhs) const {
    if (distance == rhs.distance) {
      return *yoonha < *rhs.yoonha;
    } else {
      return distance < rhs.distance;
    }
  }
};

int main(int argc, char ** argv)  try {
  // Check arguments
  if (--argc != 3)
    throw Error("usage: compare_yoonha ref mumdex_cand yoonha_cand");

  // Reference
  const Reference ref{argv[1]};
  const ChromosomeIndexLookup lookup{ref};

  // Input file
  const string cand_file_name{argv[2]};
  ifstream cand_file{cand_file_name.c_str()};
  if (!cand_file) throw Error("Could not open file") << cand_file_name;

  // Read Input file into candidate objects
  string line;
  vector<CandInfo> mumdex_info;
  unsigned int nline{0};
  while (cand_file) {
    ++nline;
    getline(cand_file, line);
    if (line.size()) {
      try {
        CandInfo i{line};
        mumdex_info.push_back(i);
      } catch (CandInfo::ParseError & e) {
        if (line[0] != 'f' || line[1] != 'a') {
          cerr << "parse error for line: " << nline << " " << line << endl;
        }
      }
    }
  }

  serr << "Read" << mumdex_info.size() << "mumdex candidates,"
       << CandInfo::n_strong << "of which were strong,"
       << CandInfo::n_indel << "of which were indels,"
       << CandInfo::n_strong_indel << "of which were strong indels"
       << endl;

  // Yoonha candidates
  const string yoonha_file_name{argv[3]};
  ifstream yoonha_file{yoonha_file_name.c_str()};
  if (!yoonha_file) throw Error("Could not open file") << yoonha_file_name;

  // Read Input file into candidate objects
  vector<YoonhaInfo> yoonha_info;
  nline = 0;
  while (yoonha_file) {
    ++nline;
    getline(yoonha_file, line);
    if (line.size()) {
      try {
        YoonhaInfo i{line, lookup};
        if (i.numParentsIndel == 0)
          yoonha_info.push_back(i);
      } catch (YoonhaInfo::ParseError & e) {
        if (line[0] != 'f' || line[1] != 'a') {
          cerr << "parse error for line: " << nline << " " << line << endl;
        }
      }
    }
  }

  serr << "Read" << yoonha_info.size() << "yoonha candidates,"
       << YoonhaInfo::n_strong << "of which were strong" << endl;

  // Group events by family
  map<string, vector<CandInfo>> mumdex_by_family;
  map<string, vector<CandInfo>> mumdex_by_family_strong;
  map<string, vector<CandInfo>> mumdex_by_family_weak;
  for (const CandInfo & info : mumdex_info) {
    if (info.indel)
      mumdex_by_family[info.family].push_back(info);
    if (info.strong_indel)
      mumdex_by_family_strong[info.family].push_back(info);
    if (info.indel && !info.strong_indel)
      mumdex_by_family_weak[info.family].push_back(info);
  }
  for (pair<const string, vector<CandInfo>> & info : mumdex_by_family) {
    sort(info.second.begin(), info.second.end());
  }
  for (pair<const string, vector<CandInfo>> & info : mumdex_by_family_strong) {
    sort(info.second.begin(), info.second.end());
  }
  for (pair<const string, vector<CandInfo>> & info : mumdex_by_family_weak) {
    sort(info.second.begin(), info.second.end());
  }
  map<string, vector<YoonhaInfo>> yoonha_by_family;
  map<string, vector<YoonhaInfo>> yoonha_by_family_strong;
  for (const YoonhaInfo & info : yoonha_info) {
    yoonha_by_family[info.family].push_back(info);
    if (info.strong)
      yoonha_by_family_strong[info.family].push_back(info);
  }
  for (pair<const string, vector<YoonhaInfo>> & info : yoonha_by_family) {
    sort(info.second.begin(), info.second.end());
  }
  for (pair<const string, vector<YoonhaInfo>> & info :
           yoonha_by_family_strong) {
    sort(info.second.begin(), info.second.end());
  }

  if (0) {
    // Quick look
    for (const CandInfo & info : mumdex_by_family_strong.begin()->second) {
      sout << info.chr << info.posA << info.invariant << endl;
    }
    sout << endl;
    for (const YoonhaInfo & info : yoonha_by_family_strong.begin()->second) {
      sout << info.chr << info.pos << info.invariant << info.variant << endl;
    }
  }

  // Randomization
  random_device rd;
  mt19937_64 mersenne{rd()};
  function<double()> unitGen{
    bind(uniform_real_distribution<double>(-0.5, 0.5), mersenne)};

  // Match candidate lists
  const unsigned int close_limit{50};
  array<map<string, vector<CandInfo>> *, 3> mumdex_lists{{
      &mumdex_by_family_strong, &mumdex_by_family, &mumdex_by_family_weak}};
  array<string, 3> mumdex_descr{{"strong short indel", "short indel",
        "weak short indel"}};
  array<map<string, vector<YoonhaInfo>> *, 2> yoonha_lists{{
      &yoonha_by_family, &yoonha_by_family_strong}};
  array<string, 2> yoonha_descr{{"all", "strong"}};
  for (unsigned int m{0}; m != mumdex_lists.size(); ++m) {
    map<string, vector<CandInfo>> & mumdex_list{*mumdex_lists[m]};
    for (unsigned int y{0}; y != yoonha_lists.size(); ++y) {
      // if (m != 0 || y != 0) continue;
      map<string, vector<YoonhaInfo>> & yoonha_list{*yoonha_lists[y]};
      vector<CandMatch> my_matches;
      unsigned int n_matches{0};
      unsigned int n_close_matches{0};
      unsigned int n_mumdex_tried{0};
      for (const pair<const string, vector<CandInfo>> & info : mumdex_list) {
        const string & family{info.first};
        const vector<CandInfo> & mumdex_candidates{info.second};
        for (const CandInfo & mumdex_cand : mumdex_candidates) {
          ++n_mumdex_tried;
          auto lower(lower_bound(yoonha_list[family].begin(),
                                 yoonha_list[family].end(),
                                 PosInfo(mumdex_cand.chr, 0)));
          unsigned int best_dist{default_distance};
          int best_invariant_dist{10000000};
          const YoonhaInfo * best_match{nullptr};
          while (lower != yoonha_list[family].end()) {
            const YoonhaInfo & match{*lower++};
            if (match.chr != mumdex_cand.chr) break;
            const unsigned int match_dist{distance(mumdex_cand, match)};
            const int invariant_dist{
              abs(match.invariant - mumdex_cand.invariant)};
            if (match_dist < best_dist ||
                (match_dist == best_dist &&
                 best_invariant_dist > invariant_dist)) {
              best_match = &match;
              best_dist = match_dist;
              best_invariant_dist = invariant_dist;
            }
          }
          if (best_match) {
            ++n_matches;
            if (best_dist < close_limit) {
              ++n_close_matches;
            }
          }
          my_matches.emplace_back(&mumdex_cand, best_match, best_dist);
        }
      }
      serr << "Matched" << n_matches << "of" << n_mumdex_tried
           << mumdex_descr[m] << "mumdex candidates to" << yoonha_descr[y]
           << "Yoon_ha candidates and" << n_close_matches << "or"
           << 100.0 * n_close_matches / n_mumdex_tried;
      cerr << "% were within";
      serr << close_limit << "bases" << endl;


      vector<YoonhaMatch> ym_matches;
      n_matches = 0;
      n_close_matches = 0;
      unsigned int n_yoonha_tried{0};
      for (const pair<const string, vector<YoonhaInfo>> & info : yoonha_list) {
        const string & family{info.first};
        const vector<YoonhaInfo> & yoonha_candidates{info.second};
        for (const YoonhaInfo & yoonha_cand : yoonha_candidates) {
          ++n_yoonha_tried;
          auto lower(lower_bound(mumdex_list[family].begin(),
                                 mumdex_list[family].end(),
                                 PosInfo(yoonha_cand.chr, 0)));
          unsigned int best_dist{default_distance};
          int best_invariant_dist{10000000};
          const CandInfo * best_match{nullptr};
          while (lower != mumdex_list[family].end()) {
            const CandInfo & match{*lower++};
            if (match.chr != yoonha_cand.chr) break;
            const unsigned int match_dist{distance(match, yoonha_cand)};
            const int invariant_dist{
              abs(match.invariant - yoonha_cand.invariant)};
            if (match_dist < best_dist ||
                (match_dist == best_dist &&
                 best_invariant_dist > invariant_dist)) {
              best_match = &match;
              best_dist = match_dist;
              best_invariant_dist = invariant_dist;
            }
          }
          if (best_match) {
            ++n_matches;
            if (best_dist < close_limit) {
              ++n_close_matches;
            }
          }
          ym_matches.emplace_back(&yoonha_cand, best_match, best_dist);
        }
      }
      serr << "Matched" << n_matches << "of" << n_yoonha_tried
           << yoonha_descr[y] << "yoonha candidates to" << mumdex_descr[m]
           << "MUMdex candidates and" << n_close_matches << "or"
           << 100.0 * n_close_matches / n_yoonha_tried;
      cerr << "% were within";
      serr << close_limit << "bases" << endl;

      if (m == 0 && y == 0) {
        sort(my_matches.begin(), my_matches.end());
        Marker marker{paa::circle(), 0.2, "0 0 0", 1, true};
        PSDoc ps{"mumdex_yoonha", "MUMdex - Yoon-Ha Comparison"};
        paa::PSHSeries<int, uint64_t> dist_hist{ps,
              ";Event Distance (<20);N",
              Bounds{0.0, 20.0}, 20};
        PSXYSeries distance_vs_rank{ps, ";Distance Rank;Distance",
              Bounds(0, my_matches.size() + 1, 0, 25), marker};
        PSXYSeries bridge_vs_rank{ps, ";Distance Rank;Bridge Count",
              Bounds(0, my_matches.size() + 1, 0, 40), marker};
        PSXYSeries nfam_vs_rank{ps, ";Distance Rank;Number of Families",
              Bounds(0, my_matches.size() + 1, 0, 6), marker};
        PSXYSeries nsam_vs_rank{ps, ";Distance Rank;Number of Samples",
              Bounds(0, my_matches.size() + 1, 0, 15), marker};
        PSXYSeries ninfam_vs_rank{ps, ";Distance Rank;Number in Family",
              Bounds(0, my_matches.size() + 1, 0, 3), marker};
        PSXYSeries inv_vs_rank{ps, ";Distance Rank;Invariant",
              Bounds(0, my_matches.size() + 1, min_inv - 1, max_inv + 1),
              marker};
        PSXYSeries off_vs_rank{ps, ";Distance Rank;Anchor Offset",
              Bounds(0, my_matches.size() + 1, -75, 55),
              marker};
        PSXYSeries off_vs_inv{ps, ";Invariant;Anchor Offset",
              Bounds(min_inv - 1, max_inv + 1, -75, 55),
              marker};
        PSXYSeries sup_vs_rank{ps, ";Distance Rank;Minimum MUM support",
              Bounds(0, my_matches.size() + 1, 0, 152), marker};
        PSXYSeries msup_vs_rank{ps, ";Distance Rank;Minimum Mate MUM support",
              Bounds(0, my_matches.size() + 1, 0, 152), marker};
        PSXYSeries msupc_vs_rank{ps,
              ";Distance Rank;Minimum Mate support count",
              Bounds(0, my_matches.size() + 1, 0, 20), marker};
        PSXYSeries parcov_vs_rank{ps, ";Distance Rank;Minimum Parent Coverage",
              Bounds(0, my_matches.size() + 1, 0, 200), marker};
        PSXYSeries map_vs_rank{ps, ";Distance Rank;Minimum Anchor Mappability",
              Bounds(0, my_matches.size() + 1, 10, 100), marker};
        PSXYSeries max_map_vs_rank{ps,
              ";Distance Rank;Maximum Anchor Mappability",
              Bounds(0, my_matches.size() + 1, 10, 100), marker};
        PSXYSeries emap_vs_rank{ps,
              ";Distance Rank;Minimum Anchor Excess Mappability",
              Bounds(0, my_matches.size() + 1, 0, 130), marker};
        PSXYSeries max_emap_vs_rank{ps,
              ";Distance Rank;Maximum Anchor Excess Mappability",
              Bounds(0, my_matches.size() + 1, 0, 130), marker};
        PSXYSeries anc_vs_rank{ps,
              ";Distance Rank;Total Parent Anchor Count",
              Bounds(0, my_matches.size() + 1, 0, 300), marker};
        PSXYSeries inv_vs_inv{ps, ";Yoon-Ha Invariant;MUMdex Invariant",
              Bounds(-75, 55, -75, 55), marker};
        PSXYSeries b1A_vs_rank{ps,
              ";Distance Rank;Number of Alignments Allowing 1 Mismatch",
              Bounds(0, my_matches.size() + 1, unset(), nunset()), marker};
        b1A_vs_rank.parents().front()->log_y(true);
        PSXYSeries b2A_vs_rank{ps,
              ";Distance Rank;Number of Alignments Allowing 2 Mismatches",
              Bounds(0, my_matches.size() + 1, unset(), nunset()), marker};
        b2A_vs_rank.parents().front()->log_y(true);
        PSXYSeries b3A_vs_rank{ps,
              ";Distance Rank;Number of Alignments Allowing 3 Mismatches",
              Bounds(0, my_matches.size() + 1, unset(), nunset()), marker};
        b3A_vs_rank.parents().front()->log_y(true);
        PSXYSeries rlen_vs_rank{ps,
              ";Distance Rank;Total Repeat Length",
              Bounds(0, my_matches.size() + 1, -1, 100), marker};
        PSXYSeries mlen_vs_rank{ps,
              ";Distance Rank;Repeat Motif Length",
              Bounds(0, my_matches.size() + 1, -1, 50), marker};
        PSXYSeries mn_vs_rank{ps,
              ";Distance Rank;Repeat Motif Copies",
              Bounds(0, my_matches.size() + 1, -1, 60), marker};

        for (unsigned int r{0}; r != my_matches.size(); ++r) {
          const CandMatch match{my_matches[r]};
          dist_hist.add_point(match.distance);
          distance_vs_rank.add_point(r + 1, match.distance + unitGen());
          bridge_vs_rank.add_point(r + 1, match.mumdex->bridge + unitGen());
          nfam_vs_rank.add_point(r + 1, match.mumdex->nFam + unitGen());
          nsam_vs_rank.add_point(r + 1, match.mumdex->nSam + unitGen());
          ninfam_vs_rank.add_point(r + 1, match.mumdex->nInFam + unitGen());
          inv_vs_rank.add_point(r + 1, match.mumdex->invariant + unitGen());
          off_vs_rank.add_point(r + 1, match.mumdex->offset + unitGen());
          off_vs_inv.add_point(match.mumdex->invariant + unitGen(),
                               match.mumdex->offset + unitGen());
          sup_vs_rank.add_point(
              r + 1, min(match.mumdex->supA, match.mumdex->supB) + unitGen());
          msup_vs_rank.add_point(
              r + 1, min(match.mumdex->mSA, match.mumdex->mSB) + unitGen());
          msupc_vs_rank.add_point(
              r + 1, min(match.mumdex->mCA, match.mumdex->mCB) + unitGen());
          parcov_vs_rank.add_point(r + 1, match.mumdex->minParCov + unitGen());
          map_vs_rank.add_point(
              r + 1, min(match.mumdex->mapA, match.mumdex->mapB) + unitGen());
          max_map_vs_rank.add_point(
              r + 1, max(match.mumdex->mapA, match.mumdex->mapB) + unitGen());
          emap_vs_rank.add_point(
              r + 1, min(match.mumdex->supA - match.mumdex->mapA,
                         match.mumdex->supB - match.mumdex->mapB) + unitGen());
          max_emap_vs_rank.add_point(
              r + 1, max(match.mumdex->supA - match.mumdex->mapA,
                         match.mumdex->supB - match.mumdex->mapB) + unitGen());
          anc_vs_rank.add_point(
              r + 1, accumulate(match.mumdex->anchorsA,
                                match.mumdex->anchorsA + 2, 0) +
              accumulate(match.mumdex->anchorsB,
                         match.mumdex->anchorsB + 2, 0) + unitGen());
          b1A_vs_rank.add_point(
              r + 1, max(match.mumdex->b1A + match.mumdex->b1B, 1U) +
              unitGen());
          b2A_vs_rank.add_point(
              r + 1, max(match.mumdex->b2A + match.mumdex->b2B, 1U) +
              unitGen());
          b3A_vs_rank.add_point(
              r + 1, max(match.mumdex->b3A + match.mumdex->b3B, 1U) +
              unitGen());
          if (match.yoonha && match.distance < 20)
            inv_vs_inv.add_point(match.mumdex->invariant + unitGen(),
                                 match.yoonha->invariant + unitGen());
          rlen_vs_rank.add_point(r + 1, match.mumdex->rlen + unitGen());
          mlen_vs_rank.add_point(r + 1, match.mumdex->mlen + unitGen());
          mn_vs_rank.add_point(r + 1, match.mumdex->mn + unitGen());
        }
      }

      if (m == 1 && y == 1) {
        sort(ym_matches.begin(), ym_matches.end());
        Marker marker{paa::circle(), 0.2, "0 0 0", 1, true};
        PSDoc ps{"yoonha_mumdex", "Yoon-Ha - MUMdex Comparison"};
        paa::PSHSeries<int, uint64_t> dist_hist{ps,
              ";Event Distance (<20);N",
              Bounds{0.0, 20.0}, 20};
        PSXYSeries distance_vs_rank{ps, ";Distance Rank;Distance",
              Bounds(0, ym_matches.size() + 1, 0, 25), marker};
        PSXYSeries count_vs_rank{ps, ";Distance Rank;Count",
              Bounds(0, ym_matches.size() + 1, 0, 100), marker};
        PSXYSeries state_score_vs_rank{ps, ";Distance Rank;State Score",
              Bounds(0, ym_matches.size() + 1, 0, 400), marker};
        PSXYSeries denovo_score_vs_rank{ps, ";Distance Rank;Denovo Score",
              Bounds(0, ym_matches.size() + 1, 50, 200), marker};
        PSXYSeries chi2_score_vs_rank{ps, ";Distance Rank;Chi2 Score",
              Bounds(0, ym_matches.size() + 1, 0, 50), marker};
        PSXYSeries nFam_vs_rank{ps, ";Distance Rank;Number of Families",
              Bounds(0, ym_matches.size() + 1, 0, 100), marker};
        PSXYSeries nPar_vs_rank{ps, ";Distance Rank;Number of Parents",
              Bounds(0, ym_matches.size() + 1, 0, 100), marker};
        PSXYSeries tFam_vs_rank{ps, ";Distance Rank;Total Family Count",
              Bounds(0, ym_matches.size() + 1, 0, 500), marker};
        PSXYSeries tPar_vs_rank{ps, ";Distance Rank;Total Parent Count",
              Bounds(0, ym_matches.size() + 1, 0, 300), marker};
        for (unsigned int r{0}; r != ym_matches.size(); ++r) {
          const YoonhaMatch match{ym_matches[r]};
          dist_hist.add_point(match.distance);
          distance_vs_rank.add_point(r + 1, match.distance + unitGen());
          // count_vs_rank.add_point(r + 1, match.yoonha->Count + unitGen());
          state_score_vs_rank.add_point(
              r + 1, match.yoonha->bestStateScore + unitGen());
          denovo_score_vs_rank.add_point(
              r + 1, match.yoonha->denovoScore + unitGen());
          chi2_score_vs_rank.add_point(
              r + 1, match.yoonha->chi2Score + unitGen());
          nFam_vs_rank.add_point(
              r + 1, match.yoonha->numFamIndel + unitGen());
          nPar_vs_rank.add_point(
              r + 1, match.yoonha->numParentsIndel + unitGen());
          tFam_vs_rank.add_point(
              r + 1, match.yoonha->totalFamCounts + unitGen());
          tPar_vs_rank.add_point(
              r + 1, match.yoonha->totalParentsCount + unitGen());
        }
      }
    }
  }



  std::cerr << "done" << endl;

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
