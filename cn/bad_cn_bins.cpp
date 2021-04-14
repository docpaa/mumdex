//
// bad_cn_bins.cpp
//
// Set bad cn bins based on population variability, etc
//
// Copyright 2016 Peter Andrews @ CSHL
//

#include <algorithm>
#include <exception>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "cn.h"
#include "error.h"
#include "mumdex.h"
#include "psplot.h"
#include "stats.h"
#include "paastrings.h"

using std::cin;
using std::cout;
using std::cerr;
using std::endl;
using std::exception;
using std::ifstream;
using std::map;
using std::ofstream;
using std::ostringstream;
using std::pair;
using std::string;
using std::vector;

using paa::remove_substring;
using paa::Bin;
using paa::Bounds;
using paa::ChromosomeIndexLookup;
using paa::CN_abspos;
using paa::CN_Bin;
using paa::CN_Bins;
using paa::Error;
using paa::MAD;
using paa::Marker;
using paa::PSDoc;
using paa::PSGraph;
using paa::PSHSeries;
using paa::PSXYSeries;
using paa::Reference;

int main(int argc, char* argv[]) try {
  if (--argc < 3)
    throw Error("usage: bad_cn_bins ref bin_file results_file ...");

  // Process arguments
  const Reference ref{argv[1]};
  const ChromosomeIndexLookup lookup{ref};
  const string bins_name{argv[2]};
  const vector<Bin> all_bins{load_bins(bins_name, ref)};
  // List of previously determined good bins
  const vector<Bin> bins{[&all_bins]() {
      vector<Bin> result;
      for (const Bin & bin : all_bins) {
        if (!bin.bad()) {
          result.push_back(bin);
        }
      }
      return result;
    }()};

  argc -= 2;
  argv += 3;
  const vector<CN_Bins> results{[argv, argc, &bins] () {
      vector<CN_Bins> result;
      result.reserve(argc);
      for (int a{0}; a != argc; ++a) {
        result.emplace_back(argv[a]);
        if (result.back().size() != bins.size()) {
          throw Error("Bin size mismatch");
        }
      }
      return result;
    }()};

  const string out_base{remove_substring(bins_name, ".txt") + ".cut"};

  // Graphs
  PSDoc ps{out_base};
  ps.pdf(false);
  const Marker small_red_marker{paa::circle(), 0.1, "1 0 0", 1, true};

  // Ratio MAD histogram
  PSGraph mad_hist_graph{ps, ";Ratio MAD;N", Bounds(0, 0.2)};
  PSHSeries<double, unsigned int> mad_hist{mad_hist_graph, 200, "1 0 0", false};

  // Ratio MAD vs abspos
  PSGraph mad_vs_abspos_graph{ps, ";Absolute Position (GB); Ratio MAD"};
  mad_vs_abspos_graph.log_y(true);
  // mad_vs_abspos_graph.range().yh(0.5);
  PSXYSeries mad_vs_abspos{mad_vs_abspos_graph, small_red_marker};

  // Get result statistics
  vector<MAD> ratio_mads;
  ratio_mads.reserve(bins.size());
  vector<double> ratios;
  ratios.reserve(results.size());
  const unsigned int x_chr{ref.find_x_chromosome()};
  const unsigned int y_chr{ref.find_y_chromosome()};
  vector<vector<double>> mads(3);
  const double billion{1000000000};
  vector<double> position_mads(bins.size());
  vector<unsigned int> position_regions(bins.size());
  for (unsigned int b{0}; b != bins.size(); ++b) {
    ratios.clear();
    for (unsigned int r{0}; r != results.size(); ++r) {
      const CN_Bin & result_bin{results[r][b]};
      ratios.push_back(result_bin.ratio());
    }
    sort(ratios.begin(), ratios.end());
    ratio_mads.emplace_back(ratios);
    const MAD & mad{ratio_mads.back()};
    const Bin & bin{bins[b]};
    // Fill graphs
    mad_hist.add_point(mad.mad());
    mad_vs_abspos.add_point(bin.abspos() / billion, mad.mad());
    const unsigned int region{bin.chromosome() == x_chr ? 1U :
          (bin.chromosome() == y_chr ? 2U : 0U)};
    mads[region].push_back(mad.mad());
    position_mads[b] = mad.mad();
    position_regions[b] = region;
  }

  // Mad of mads...
  vector<MAD> mads_mads;
  const CN_abspos cn_abspos{ref};
  const unsigned int region_chromosome_start[3]{
    cn_abspos.chromosomes().front(), x_chr, y_chr};
  vector<double> cutoffs;
  for (unsigned int c{0}; c != mads.size(); ++c) {
    vector<double> & vals{mads[c]};
    sort(vals.begin(), vals.end());
    mads_mads.emplace_back(vals);
    const unsigned int chr{region_chromosome_start[c]};
    const MAD & mads_mad{mads_mads.back()};
    const unsigned int start{cn_abspos.ref_offset(chr)};
    const unsigned int stop{chr == cn_abspos.chromosomes().front() ?
          cn_abspos.ref_offset(x_chr) : cn_abspos.ref_offset(chr + 1)};
    const double cutoff{mads_mad.median() + 6 * mads_mad.mad()};
    auto lower = lower_bound(vals.begin(), vals.end(), cutoff);
    cout << ref.name(chr) << " " << start << " " << stop << " "
         << mads_mad.median() << " " << cutoff << " "
         << 1.0 * (lower - vals.begin()) / vals.size() << endl;

    ostringstream cutoff_line;
    cutoff_line << "0 0 0 c 2 lw np " << start / billion << " "
                << cutoff << " gc m "
                << stop / billion << " " << cutoff << " gc l sp";
    mad_vs_abspos_graph.ps(cutoff_line.str());
    cutoffs.push_back(cutoff);
  }

  // Output bins
  const string out_name{out_base + ".txt"};
  ofstream out{out_name.c_str()};
  if (!out) throw Error("Problem opening output bins file") << out_name;
  unsigned int gb{0};
  for (unsigned int b{0}; b != all_bins.size(); ++b) {
    const Bin & bin{all_bins[b]};
    if (bin.bad()) {
      bin.output(out, ref);
    } else {
      bin.output(out, ref, position_mads[gb] > cutoffs[position_regions[gb]]);
      ++gb;
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
