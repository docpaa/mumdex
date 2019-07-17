//
// bad_haha_loci.cpp
//
// find properties of bad haha loci
//
// Copyright 2019 Peter Andrews @ CSHL
//

#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "error.h"
#include "mumdex.h"
#include "psplot.h"
#include "utility.h"

using std::cerr;
using std::cin;
using std::cout;
using std::endl;
using std::exception;
using std::ifstream;
using std::map;
using std::max;
using std::pair;
using std::set;
using std::string;
using std::vector;

using paa::Bounds;
using paa::ChromosomeIndexLookup;
using paa::Error;
using paa::Mappability;
using paa::Marker;
using paa::PSDoc;
using paa::PSGraph;
using paa::PSHSeries;
using paa::PSPage;
using paa::PSXYSeries;
using paa::Reference;

template <class Value>
class CumVs {
 public:
  using Hist = PSHSeries<Value, uint64_t>;
  CumVs(PSDoc & plots,
        const std::string & value_name,
        const std::string & name1,
        const std::string & name2,
        const Bounds & bounds = Bounds{},
        const unsigned int n_bins = 0) :
      page{plots, "", "1 2"},
      hist1{page, name1 + ";" + value_name + ";N", bounds,
            (n_bins ? n_bins :
             static_cast<unsigned int>(bounds.xh() - bounds.xl())),
            "1 0 0", true},
      hist2{page, name2 + ";" + value_name + ";N", bounds,
            (n_bins ? n_bins :
             static_cast<unsigned int>(bounds.xh() - bounds.xl())),
            "1 0 0", true} {}
  void add1(const Value value) { hist1.add_point(value); }
  void add2(const Value value) { hist2.add_point(value); }

 private:
  PSPage page;
  Hist hist1;
  Hist hist2;
};
using iCumVs = CumVs<unsigned int>;
using dCumVs = CumVs<double>;

int main(int argc, char ** argv) try {
  paa::exit_on_pipe_close();
  std::ios_base::sync_with_stdio(false);
  cin.tie(nullptr);

  const string usage{
    "usage: bad_haha_loci ref loci_info good_loci bad_loci min_bad"};
  if (--argc != 5) throw Error(usage);

  const Reference ref{argv[1]};
  const Mappability mappability{ref};
  const ChromosomeIndexLookup chr_lookup{ref};

  // Read input info on how hets were chosen
  const string loci_info_name{argv[2]};
  ifstream loci_info_file{loci_info_name.c_str()};
  if (!loci_info_file)
    throw Error("Problem opening loci info file") << loci_info_name;
  loci_info_file.ignore(100000, '\n');
  string dummy;
  using ChrPos = pair<string, unsigned int>;
  ChrPos chr_pos;
  using LocusInfo = pair<unsigned int, double>;
  LocusInfo locus_info;
  map<ChrPos, LocusInfo> loci_info;
  while (loci_info_file >> chr_pos.first >> chr_pos.second
         >> dummy >> dummy >> dummy >> dummy
         >> dummy >> dummy >> dummy >> dummy
         >> locus_info.first >> locus_info.second) {
    --chr_pos.second;
    if (false)
      cout << chr_pos.first << " " << chr_pos.second
           << " " << locus_info.first << " " << locus_info.second << '\n';
    loci_info[chr_pos] = locus_info;
  }

  // considered loci
  set<ChrPos> considered_loci;

  // Read info on good loci
  const string good_loci_name{argv[3]};
  ifstream good_loci_file{good_loci_name.c_str()};
  if (!good_loci_file)
    throw Error("Problem opening good loci file") << good_loci_name;
  good_loci_file.ignore(100000, '\n');
  unsigned int bin_size;
  unsigned int cover;
  map<ChrPos, unsigned int> good_loci;
  while (good_loci_file >> chr_pos.first >> chr_pos.second >> bin_size
         >> cover) {
    if (false)
      cout << chr_pos.first << " " << chr_pos.second << " " << bin_size << '\n';
    good_loci[chr_pos] = bin_size;
    considered_loci.insert(chr_pos);
  }

  // Read info on bad loci
  const string bad_loci_name{argv[4]};
  const unsigned int min_bad{static_cast<unsigned int>(atoi(argv[5]))};
  ifstream bad_loci_file{bad_loci_name.c_str()};
  if (!bad_loci_file)
    throw Error("Problem opening bad loci file") << bad_loci_name;
  using BadInfo = pair<unsigned int, unsigned int>;
  map<ChrPos, BadInfo> bad_loci;
  BadInfo bad_info;
  unsigned int extra;
  while (bad_loci_file >> bad_info.first
         >> chr_pos.first >> chr_pos.second
         >> bad_info.second >> extra) {
    if (false) cout << chr_pos.first << " " << chr_pos.second
                    << " " << bad_info.first << " " << bad_info.second << '\n';
    bad_loci[chr_pos] = bad_info;
    auto inserted = considered_loci.insert(chr_pos);
    if (inserted.second == false || bad_info.first < min_bad)
      considered_loci.erase(inserted.first);
  }

  // Generate output, mostly plots
  cerr << "Considering " << considered_loci.size()
       << " with " << bad_loci.size() << " bad loci "
       << "and " << good_loci.size() << " good loci "
       << "of " << loci_info.size() << endl;
  PSDoc plots{"bad_loci"};
  // Mappability plots
  using uHist = PSHSeries<unsigned int, unsigned int>;
  const unsigned int max_map{255};
  PSPage map_page{plots, "", "1 2"};
  uHist good_map{map_page, "Good loci;Mappability;N", Bounds{0, max_map},
        max_map, "1 0 0", true};
  uHist bad_map{map_page, "Bad loci;Mappability;N", Bounds{0, max_map},
        max_map, "1 0 0", true};
  iCumVs map_vs{plots, "Mappability", "Good loci", "Bad loci",
        Bounds{0, max_map}};

  // Size plots
  const unsigned int max_size{10000000};
  PSPage size_page{plots, "", "1 2"};
  uHist good_size{size_page, "Good loci;Bin Size;N", Bounds{0, max_size},
        100, "1 0 0", true};
  uHist bad_size{size_page, "Bad loci;Bin Size;N", Bounds{0, max_size},
        100, "1 0 0", true};
  // Total count plots
  PSPage total_page{plots, "", "1 2"};
  uHist good_total{total_page, "Good loci;Pop Count;N",
        Bounds{1875, 2025}, 75};
  uHist bad_total{total_page, "Bad loci;Pop Count;N",
        Bounds{1875, 2025}, 75};
  // p-value plots
  using dHist = PSHSeries<double, unsigned int>;
  PSPage pval_page{plots, "", "1 2"};
  dHist good_pval{pval_page, "Good loci;P-Value;N", Bounds{-0.1, 1.1},
        100, "1 0 0", true};
  dHist bad_pval{pval_page, "Bad loci;P-Value;N", Bounds{-0.1, 1.1},
        100, "1 0 0", true};
  PSPage score_page{plots, "", "1 2"};
  const double max_score{5};
  dHist good_score{score_page, "Good loci;Score;N", Bounds{0, max_score},
        100, "1 0 0", true};
  dHist bad_score{score_page, "Bad loci;Score;N", Bounds{0, max_score},
        100, "1 0 0", true};
  for (const ChrPos locus : considered_loci) {
    const unsigned int chr{chr_lookup[locus.first]};
    const unsigned int abspos{ref.abspos(chr, locus.second)};
    const unsigned int map_val{
      max(mappability.low(abspos), mappability.high(abspos))};
    locus_info = loci_info.at(locus);
    auto bad_result = bad_loci.find(locus);
    if (bad_result != bad_loci.end()) {
      bad_info = bad_result->second;
      // const unsigned int bad_count{bad_info.first};
      map_vs.add2(map_val);
      bad_map.add_point(map_val);
      bin_size = bad_info.second;
      bad_size.add_point(std::min(bin_size, max_size - 1));
      bad_total.add_point(locus_info.first);
      bad_pval.add_point(locus_info.second);
      bad_score.add_point(-log10(locus_info.second));
    } else {
      map_vs.add1(map_val);
      good_map.add_point(map_val);
      bin_size = good_loci.at(locus);
      good_size.add_point(std::min(bin_size, max_size - 1));
      good_total.add_point(locus_info.first);
      good_pval.add_point(locus_info.second);
      good_score.add_point(-log10(locus_info.second));
    }
  }

  return 0;
} catch (Error & e) {
  cerr << "paa::Error:" << endl;
  cerr << e.what() << endl;
  return 1;
} catch(exception & e) {
  cerr << e.what() << endl;
  return 1;
} catch(...) {
  cerr << "Some exception was caught." << endl;
  return 1;
}
