//
// plot.cpp
//
// plot series from a text file with format:
//
// x_name y1_name y2_name
// data...
//
// Copyright 2016 Peter Andrews @ CSHL
//

#include <iostream>
#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "error.h"
#include "psplot.h"
#include "utility.h"

using std::cout;
using std::cerr;
using std::endl;
using std::ifstream;
using std::istringstream;
using std::ostringstream;
using std::string;
using std::vector;

using paa::sout;
using paa::Bounds;
using paa::Error;
using paa::Marker;
using paa::PSDoc;
using paa::PSPage;
using paa::PSGraph;
using paa::PSXYSeries;

int main(int argc, char ** argv)  try {
  // Check arguments
  --argc;
  const string usage{
    "usage: plot -[log_y|legend|lines]... data_file title;x;y"
        " [xlow xhigh ylow yhigh]"};
  if (argc < 2) throw Error(usage);
  bool log_y{false};
  bool legend{false};
  bool acted {true};
  bool lines{false};
  while (acted) {
    acted = false;
    if (string(argv[1]) == "-log_y") log_y = acted = true;
    if (string(argv[1]) == "-legend") legend = acted = true;
    if (string(argv[1]) == "-lines") lines = acted = true;
    if (acted) {
      --argc;
      ++argv;
    }
  }
  if (argc != 2 && argc != 6) throw Error(usage);

  // Interpret arguments
  const string input_name{argv[1]};
  ifstream input{input_name};
  if (!input) throw Error("Problem opening data file") << input_name;
  const string title{argv[2]};
  const double x_low{argc == 6 ? strtod(argv[3], nullptr) : paa::unset()};
  const double x_high{argc == 6 ? strtod(argv[4], nullptr) : paa::nunset()};
  const double y_low{argc == 6 ? strtod(argv[5], nullptr) : paa::unset()};
  const double y_high{argc == 6 ? strtod(argv[6], nullptr) : paa::nunset()};

  // Graph
  PSDoc ps{input_name, input_name};
  ps.pdf(true);
  PSGraph graph{ps, title, Bounds{x_low, x_high, y_low, y_high}};
  if (log_y) graph.log_y(true);
  //  graph.ps(  // "np 0 0 gfc m 1 0 gfc l 1 1 gfc l 1 0 gfc l cp clip "
  // "0 0 0 c 0.5 lw np 0 0 gc m 4 4 gc l sp");
  // graph.log_x(true).log_y(true);
  // graph.log_x(true);
  vector<PSXYSeries> series;
  series.reserve(100);

  // Default markers
  const vector<string> colors{
    "0 0 0",  // black
        "0 0.717647 0",  // green
        "0.898039 0 0",  // red
        "0.145098 0.145098 1",  // blue
        "0.909804 0.690196 0",  // orange
        "0.439216 0.564706 0.564706",  // grey
        "0.972549 0.282353 0.878431",  // pink
        "0 0.909804 0.972549",  // light blue
        "0.38 0.2 0.07",  // brown
        "0.5 0.5 1",  //
        };
  const vector<string (*)()> shapes{
    paa::circle, paa::square, paa::triangle,
        paa::diamond, paa::circle, paa::plus,
        paa::xshape, paa::circle, paa::star};
  vector<Marker> markers;

  // Series names
  vector<string> series_names;
  string header;
  getline(input, header);
  istringstream header_stream{header.c_str()};
  header_stream >> header;
  unsigned int color{0};
  unsigned int shape{0};
  while (header_stream >> header) {
    series_names.push_back(header);
    series.emplace_back(graph, Marker{shapes[shape++ % shapes.size()](),
            1.0, colors[color++ % colors.size()]});  // , header);
    if (lines) series.back().do_lines(true);
  }

  // Read data
  double x;
  double y;
  while (input) {
    input >> x;
    if (!input) break;
    for (unsigned int s{0}; s != series.size(); ++s) {
      input >> y;
      // cerr << "add " << s << " " << x << " " << y << endl;
      series[s].add_point(x, y);
    }

    if (!input) throw Error("Problem parsing input");
  }

  if (legend) {
    ostringstream legend_stream;
    legend_stream << "70 15 m 15 sf (Legend:) s";
    for (uint64_t c{0}; c != series_names.size(); ++c) {
      if (c) legend_stream << " 0 0 0 c (,) s";
      legend_stream << " " << colors[c % colors.size()] << " c ( "
                    << series_names[c] << ") s gs 10 5 rm cxy tr "
                    << series[c].marker(1.0).commands() << " sp gr"
                    << " 17 0 rm ";
    }
    graph.ps(legend_stream.str());
  }

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
