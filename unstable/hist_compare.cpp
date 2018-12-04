//
// hist_compare
//
// read in table, plot hists split by first column
//
// Copyright 2018 Peter Andrews @ CSHL
//

#include <algorithm>
#include <exception>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "error.h"
#include "psplot.h"

using std::cerr;
using std::cin;
using std::cout;
using std::endl;
using std::exception;
using std::istringstream;
using std::string;
using std::vector;

using paa::Bounds;
using paa::Error;
using paa::PSDoc;
using paa::PSHSeries;
using paa::PSPage;

int main(int argc, char* argv[])  try {
  --argc;
  if (argc != 1) throw Error("data_file | hist_compare out_name");

  const string out_name{argv[1]};

  string line;
  getline(cin, line);
  istringstream header_stream{line.c_str()};
  vector<string> header;
  while (header_stream >> line) header.push_back(line);
  vector<vector<double>> data(header.size());

  while (getline(cin, line)) {
    istringstream line_stream{line.c_str()};
    for (uint64_t c{0}; c != data.size(); ++c) {
      vector<double> & vals{data[c]};
      double input;
      line_stream >> input;
      if (!line_stream) throw Error("Bad input from line") << line;
      vals.push_back(input);
    }
    if (!line_stream) throw Error("Parse Error");
  }

  PSDoc plots{out_name};
  plots.pdf(true);

  using Hist = PSHSeries<int, unsigned int>;
  for (uint64_t c{1}; c != data.size(); ++c) {
    const bool suppress_zeros{header[c] == "num_echos"};
    vector<double> & values{data[c]};
    auto minmax = minmax_element(values.begin(), values.end());
    const double min{*minmax.first};
    const double max{*minmax.second};
    PSPage * page = PSPage::create(
        plots, header[c] +
        (suppress_zeros ? " (zeros suppressed)" : ""), "1 2");
    const unsigned int n_bins{static_cast<unsigned int>(max - min + 1.5)};
    Hist * hist0 = Hist::create(
        *page, header[0] + " = false;" + header[c] + ";N",
        Bounds{min, max + 1}, n_bins);
    Hist * hist1 = Hist::create(
        *page, header[0] + " = true;" + header[c] + ";N",
        Bounds{min, max + 1}, n_bins);
    vector<Hist *> hists{hist0, hist1};
    for (unsigned int v{0}; v != values.size(); ++v)
      if (fabs(values[v]) > 0.001 || !suppress_zeros)
        hists[data[0][v]]->add_point(values[v]);
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



