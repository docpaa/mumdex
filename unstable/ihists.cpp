//
// ihists
//
// make histograms for integer input columns
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
  if (argc != 1) throw Error("data_file | ihists out_name");

  const string out_name{argv[1]};

  string line;
  getline(cin, line);
  istringstream header_stream{line.c_str()};
  vector<string> header;
  while (header_stream >> line) header.push_back(line);
  vector<vector<int64_t>> data(header.size());

  while (getline(cin, line)) {
    istringstream line_stream{line.c_str()};
    for (uint64_t c{0}; c != data.size(); ++c) {
      vector<int64_t> & vals{data[c]};
      int64_t input;
      line_stream >> input;
      if (!line_stream) throw Error("Bad input from line") << line;
      vals.push_back(input);
    }
    if (!line_stream) throw Error("Parse Error");
    line_stream.get();
    if (line_stream) throw Error("unexpected character");
  }

  PSDoc plots{out_name};
  plots.pdf(true);

  using Hist = PSHSeries<int64_t, uint64_t>;
  for (uint64_t c{0}; c != data.size(); ++c) {
    vector<int64_t> & values{data[c]};
    auto minmax = minmax_element(values.begin(), values.end());
    const int64_t min{*minmax.first};
    const int64_t max{*minmax.second};
    PSPage * page{PSPage::create(plots)};
    const unsigned int n_bins{static_cast<unsigned int>(max - min + 1)};
    Hist * hist{Hist::create(*page, header[c] + ";" + header[c] + ";N",
                             Bounds{static_cast<double>(min),
                                   static_cast<double>(max) + 1}, n_bins)};
    for (unsigned int v{0}; v != values.size(); ++v)
      hist->add_point(values[v]);
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



