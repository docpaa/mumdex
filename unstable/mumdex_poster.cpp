//
// mumdex_poster
//
// Generate a poster describing mumdex
//
// Copyright 2018 Peter Andrews @ CSHL
//

#include <exception>
#include <iostream>
#include <fstream>
#include <sstream>

#include "error.h"
#include "layout.h"

using std::cout;
using std::cerr;
using std::endl;
using std::exception;
using std::ofstream;
using std::ostringstream;

using paa::Error;
using paa::Geometry;
using paa::Point;

constexpr int inch{72};
constexpr int foot{12};
constexpr int points_per_foot{inch * foot};
constexpr int page_width{4 * points_per_foot};
constexpr int page_height{3 * points_per_foot};
constexpr Geometry page_geometry{page_width, page_height};
constexpr int header_height{4 * inch};
constexpr int header_pos{page_height - header_height};
constexpr int title_font_size{header_height / 2};
constexpr int ncol{3};
constexpr int colw{page_width / 3};

int main(int argc, char**) try {
  if (--argc != 0) throw Error("usage: mumdex_poster");

  ostringstream poster;
  poster
      << ""
      << "%!PS-Adobe-3.0\n"
      << "%%BoundingBox: 0 0 " << page_width << " " << page_height << "\n"
      << "%%Pages: 1\n"
      << "%%Title: MUMdex presentation\n"
      << "%%Creator: Peter Andrews\n"
      << "%%CreationDate: Today\n"
      << "%%EndComments\n"
      << "%%BeginProlog\n"
      << "/pgw " << page_width << " def\n"
      << "/pgh " << page_height << " def\n"
      << "/inch " << inch << " def\n"
      << "/colw " << colw << " def\n"
      << "%%ENDProlog\n"
      << "%%Page: 1 1\n"
      << "%%\n";

  poster << "1 0 0 setrgbcolor "  << inch / 2 << " setlinewidth\n"
         << "newpath "
         << "0 " << header_pos << " moveto\n" << page_width << " 0 rlineto\n";
  for (int c{1}; c != ncol; ++c) {
    // poster << "gsave \n";
    poster << c * colw << " 0 moveto 0 " << header_pos << " rlineto\n";
    // poster << "grestore\n";
    // poster << "stroke\n";
  }
  poster << "stroke\n";
  poster << "showpage\n"
         << "%%EndPage: 1\n"
         << "%%EOF\n";

  std::cout << poster.str() << std::endl;
  return 0;
} catch (Error & e) {
  cerr << "Error:" << endl;
  cerr << e.what() << endl;
  return 1;
} catch (exception & e) {
  cerr << "exception" << endl;
  cerr << e.what() << endl;
  return 1;
} catch (...) {
  cerr << "unknown exception was caught" << endl;
  return 1;
}
