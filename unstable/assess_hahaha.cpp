//
// assess_hahaha
//
// assess the hahaha method
//
// Copyright 2019 Peter Andrews @ CSHL
//

#include <algorithm>
#include <fstream>
#include <iostream>
#include <memory>
#include <random>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#include <climits>

#include "cn.h"
#include "error.h"
#include "mumdex.h"
#include "psplot.h"
#include "strings.h"
#include "utility.h"

using std::cerr;
using std::cin;
using std::cout;
using std::endl;
using std::exception;
using std::ifstream;
using std::istream;
using std::make_unique;
using std::min;
using std::ostringstream;
using std::pair;
using std::string;
using std::to_string;
using std::unique_ptr;
using std::vector;

using paa::replace_substring;
using paa::Bounds;
using paa::CN_abspos;
using paa::Error;
using paa::Marker;
using paa::PSDoc;
using paa::PSGraph;
using paa::PSHSeries;
using paa::PSPage;
using paa::PSXYMSeries;
using paa::PSXYSeries;
using paa::Reference;

int main(int argc, char ** argv) try {
  paa::exit_on_pipe_close();
  std::ios_base::sync_with_stdio(false);
  cin.tie(nullptr);

  const string usage{
    "usage: assess_hahaha ref [input_files ...]"};

  // process optional arguments
  --argc;
  const double scatter_amount{0.5};
  const bool scatter{true};
  // const bool do_abs{true};

  if (argc < 1) throw Error("Not enough arguments:") << usage;
  const Reference ref{argv[1]};
  const CN_abspos cn_abspos{ref};
  const string out_name{"panels"};
  const string title{";;"};
  const uint64_t n_panels{10};
  argc -= 1;
  argv += 1;

  vector<unique_ptr<ifstream>> input_files;
  vector<istream *> inputs;
  if (argc) {
    while (argc--) {
      const string file_name{(argv++)[1]};
      input_files.push_back(make_unique<ifstream>(file_name.c_str()));
      if (!*input_files.back())
        throw Error("Could not open file for input") << file_name;
      inputs.push_back(input_files.back().get());
    }
  } else {
    inputs.push_back(&cin);
  }

  vector<pair<double, double>> xy;
  {
    double x;
    double y;
    for (istream * in : inputs)
      while ((*in) >> x >> y) xy.emplace_back(x, y);
  }

  if (xy.size() < 2) throw Error("Not enough points");

  sort(xy.begin(), xy.end());
  const double min_x{xy.front().first};
  const double max_x{xy.back().first};
  const double width{max_x - min_x};
  if (min_x >= max_x) throw Error("min equals max");

  PSDoc plots{out_name};
  PSDoc plots2{out_name + ".all"};
  plots.pdf(false);
  PSPage plot_page{plots, "", "1 " + to_string(n_panels)};
  vector<unique_ptr<PSGraph>> graphs;
  vector<unique_ptr<PSXYMSeries>> series;
  const Marker good_marker{paa::circle(), 0.1, "0 0.8 0", 0.1, true, "0 0.8 0"};
  const Marker bad_marker{paa::circle(), 0.1, "0.8 0 0", 0.1, true, "1 0.8 0"};
  for (uint64_t p{0}; p != n_panels; ++p) {
    graphs.push_back(make_unique<PSGraph>(plot_page, title));
    graphs.back()->do_ticks(false);
    series.push_back(make_unique<PSXYMSeries>(*graphs.back()));
    ostringstream ps;
    ps << "1 lw 0 0 0 c 12 sf\n";
    const double min_pos{p * width / n_panels};
    const double max_pos{(p + 1) * width / n_panels};

    for (unsigned int c{0}; c != cn_abspos.chromosomes().size(); ++c) {
      const unsigned int chr{cn_abspos.chromosomes()[c]};
      if (ref.name(chr).find("X") != string::npos) break;
      const unsigned int offset{cn_abspos.offsets()[c]};
      if (offset < min_pos) continue;
      if (offset > max_pos) break;
      ps << "np " << offset << " xc 0 yfc m " << offset << " xc 1 yfc l sp "
         << offset << " xc 1 yfc m 5 -12 rm ("
         << replace_substring(ref.name(chr), "chr", "") << ") s\n";
    }
    graphs.back()->ps(ps.str());
  }
  std::random_device rd{};
  std::mt19937_64 mersenne{rd()};
  std::function<double()> gen{
    std::bind(std::uniform_real_distribution<double>(
        -scatter_amount, scatter_amount), mersenne)};

  PSHSeries<double, uint64_t> scores{plots2, ";Score;N", 100, Bounds{0, 11}};
  for (uint64_t p{0}; p != xy.size(); ++p) {
    const double x{xy[p].first};
    const double y{xy[p].second};
    const uint64_t panel{static_cast<uint64_t>(
        min(1.0 * n_panels - 1, n_panels * (x - min_x) / width))};
    const double scatter_y{(scatter ? gen() : 0.0) + y};
    // (do_abs ? fabs(y - 0.5) : y)};
    const Marker & marker{y > 0.5 ? good_marker : bad_marker};
    if (-log10(0.5 - fabs(y - 0.5)) >= 2)
      series[panel]->add_point(x, scatter_y, marker);
    scores.add_point(std::max(10.0, -log10(0.5 - fabs(y - 0.5))));
  }

  PSGraph stats{plots2, "HAHA run statistics;Score;Fraction"};
  PSGraph success{plots2,
        "HAHA run statistics;Fraction Called;Fraction Incorrect"};
  success.log_y(true);
  const Marker called_marker{paa::circle(), 0.5, "0 0 0", 0.5, true, "0 0 0"};
  const Marker good_call_marker{
    paa::circle(), 0.5, "0.8 0 0", 0.5, true, "0.8 0 0"};

  PSXYSeries frac_called_series{stats, called_marker};
  PSXYSeries good_call_series{stats, good_call_marker};
  PSXYSeries success_series{success, good_call_marker};
  ostringstream thresh_ps;
  thresh_ps << "12 sf\n";
  for (double threshold{0}; threshold < 10; threshold += 0.1) {
    uint64_t n{0};
    uint64_t n_good{0};
    for (uint64_t p{0}; p != xy.size(); ++p) {
      const double y{xy[p].second};
      if (-log10(0.5 - fabs(y - 0.5)) >= threshold) {
        ++n;
        if (y > 0.5) ++n_good;
      }
    }
    frac_called_series.add_point(threshold, 1.0 * n / xy.size());
    good_call_series.add_point(threshold, 1.0 * n_good / n);
    success_series.add_point(1.0 * n / xy.size(), 1 - 1.0 * n_good / n);
    double integral;
    double fractional;
    fractional = modf(threshold + 0.001, &integral);
    if (fractional < 0.01)
      thresh_ps << 1.0 * n / xy.size() << " " << (1 - 1.0 * n_good / n)
                << " gc m 0 -15 rm (" << threshold << ") s\n";
  }
  success.ps(thresh_ps.str());

  return 0;
} catch (Error & e) {
  cerr << "paa::Error:" << endl;
  cerr << e.what() << endl;
  return 1;
} catch(exception & e) {
  cerr << e.what() << endl;
  return 1;
}
catch(...) {
  cerr << "Some exception was caught." << endl;
  return 1;
}
