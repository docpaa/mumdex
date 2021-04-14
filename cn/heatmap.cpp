//
// heatmap.cpp
//
// Make a copy number heatmap
//
// Copyright 2020 Peter Andrews @ CSHL
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
using std::istringstream;
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
using paa::PSPage;
using paa::PSHSeries;
using paa::PSXYSeries;
using paa::Reference;

int main(int argc, char* argv[]) try {
  if (--argc < 4)
    throw Error("usage: heatmap ref out title results_file[:name] ...");

  // Process arguments
  const Reference ref{argv[1]};
  const ChromosomeIndexLookup lookup{ref};
  const CN_abspos cn_abspos{ref};
  const string out_name{argv[2]};
  const string title{argv[3]};

  argc -= 3;
  argv += 4;

  using Spec = pair<string, string>;
  const vector<Spec> file_specs{[argv, argc]() {
      vector<Spec> result;
      for (int a{0}; a != argc; ++a) {
        istringstream spec{argv[a]};
        string file;
        string name;
        getline(spec, file, ':');
        getline(spec, name);
        if (name.empty()) name = file;
        result.emplace_back(name, file);
      }
      return result;
    }()};

  // Postscript document
  PSDoc ps{out_name};
  ps.pdf(true);
  // PSXYSeries mapping{ps, "Color mapping;CN;value"};

  // Heatmap code
  ostringstream output;
  // Layout dimensions
  const unsigned int page_width{ps.width()};
  const unsigned int page_height{ps.height()};
  const double title_font_size{25};
  const double chr_font_size{12};
  const double name_font_size{12};
  const double trim{10};
  const double name_width{50};
  const double title_height{page_height - 1.5 * title_font_size};
  const double left{trim + name_width};
  const double right{page_width - trim};
  const double bottom{trim + chr_font_size + trim};
  const double top{title_height - trim};
  const double box_width{right - left};
  const double box_height{top - bottom};
  const uint64_t n_samples{file_specs.size()};
  const double sample_height{box_height / n_samples};

  // Title line
  output << page_width / 2 << " " << title_height
         << " m " << title_font_size << " sf (" << title << ") jc s\n";
  // Heatmap data
  for (uint64_t s{0}; s != n_samples; ++s) {
    const string & seg_name{file_specs[s].second};
    const string & sample{file_specs[s].first};
    ifstream seg_file{seg_name};
    if (!seg_file) throw Error("Problem opening seg file") << seg_name;
    seg_file.ignore(10000, '\n');
    string chr_name;
    unsigned int start;
    unsigned int stop;
    unsigned int n_bins;
    unsigned int bin_start;
    unsigned int bin_stop;
    double norm_count;
    double cn;
    double norm_cn;
    uint64_t total_n_bins{0};
    double score{0};
    while (seg_file >> chr_name >> start >> stop >> n_bins
           >> bin_start >> bin_stop >> norm_count >> cn >> norm_cn) {
      seg_file.ignore(10000, '\n');
      const unsigned int chr{lookup[chr_name]};
      const double deviation{cn - norm_cn};
      total_n_bins += n_bins;
      score += n_bins * fabs(deviation);
      const unsigned int abspos_start{cn_abspos(chr, start)};
      const unsigned int abspos_stop{cn_abspos(chr, stop)};
      const double xl{
        left + box_width * abspos_start / cn_abspos.n_positions()};
      const double xr{left + box_width * abspos_stop / cn_abspos.n_positions()};
      const double b{bottom + s * box_height / n_samples};
      const double t{bottom + (s + 1) * box_height / n_samples};
      if (fabs(deviation) < 0.5) continue;
      output << " np " << xl << " " << b << " m " << xr << " " << b << " l "
             << xr << " " << t << " l " << xl << " " << t << " l cp ";
      if (deviation >= 0) {
        const double value{1 - atan2(deviation, 1) / atan2(1, 0)};
        // mapping.add_point(deviation, value);
        output << value << " " << value << " 1 c fp\n";
      } else {
        const double value{1 - (norm_cn - cn) / norm_cn};
        // mapping.add_point(deviation, value);
        output << " 1 " << value << " " << value << " c fp\n";
      }
    }
    score /= total_n_bins;
    cout << sample << " " << score << endl;
  }

  // Heatmap border
  output << " 0 0 0 c 2 lw np "
         << left << " " << bottom << " m "
         << right << " " << bottom << " l "
         << right << " " << top << " l "
         << left << " " << top << " l cp sp\n";
  // Sample names
  output << name_font_size << " sf\n";
  for (uint64_t s{0}; s != file_specs.size(); ++s) {
    const Spec & spec{file_specs[s]};
    const string & name{spec.first};
    output << left - 6 << " " << bottom + (s + 0.3) * sample_height
           << " m (" << name << ") jr s\n";
  }
  // Sample dividers
  output << "0.5 lw\n";
  for (uint64_t s{1}; s != file_specs.size(); ++s) {
    output << "np " << left << " " << bottom + s * sample_height
           << " m " << right << " " << bottom + s * sample_height
           << " l sp\n";
  }
  // Chromosome dividers and names
  output << "0.5 lw " << chr_font_size << " sf\n";
  for (const unsigned int c : cn_abspos.chromosomes()) {
    const double frac{1.0 * cn_abspos.ref_offset(c) / cn_abspos.n_positions()};
    const double x{left + frac * box_width};
    output << "np " << x << " " << bottom << " m "
           << x << " " << top << " l sp\n";
    const double half{
      1.0 * ref.size(c) * box_width / cn_abspos.n_positions() / 2};
    output << x + half << " " << trim << " m "
           << "(" << remove_substring(ref.name(c), "chr") << ") jc s\n";
  }
  PSPage heatmap_page{ps, output.str()};

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
