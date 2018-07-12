//
// fastplot
//
// Try to make a faster plot! A work in progress.....
//
// Copyright 2018 Peter Andrews @ CSHL
//

#include <exception>
#include <iostream>
#include <sstream>
#include <string>

#include "error.h"
#include "psplot.h"

namespace paa {
using PSFGraph = PSGraph;
}  // namespace paa

using std::cout;
using std::cerr;
using std::endl;
using std::exception;
using std::ostringstream;
using std::string;

using paa::Error;
using paa::Marker;
using paa::PSDoc;
using paa::PSGraph;
using paa::PSFGraph;
using paa::PSXYSeries;

int main(int argc, char** argv) try {
  if (--argc != 0) throw Error("usage: video");
  if (0) cout << argv[1];

  PSDoc doc{"faster"};
  doc.pdf(true);
  const Marker marker{paa::circle(), 0.2, "0 0 0", 1, true, "0 0 0"};
  PSGraph old_plot{doc, "Old vs New Graph Method;X;Y"};
  PSXYSeries old_series{old_plot, marker};
  const uint64_t n_points{100000};
  for (uint64_t i{0}; i != n_points; ++i) {
    old_series.add_point(i, sqrt(i));
  }
  ostringstream test;
  test << R"xxx(
1 0 0 c
/circle {
  ucache
  -1 -1 1 1 setbbox
  0 0 1 0 360 arc
  closepath
} cvlit def

/circle2 {
  {-1 -1 1 1
  0 0 1 0 360}
  <0B 00 07 0A>
} cvlit def

/op {
translate circle2 ufill
} def

/circle3
<<
/FormType 1
/PaintProc {
  pop
  0 0 1 0 360 arc
  closepath fill
} bind

/BBox [-1 -1 1 1]
/Matrix [1 0 0 1 0 0]
>> def

0 0 gc translate

)xxx";
  //   {11 0 7 10}
  const double scale{0.0002};
  double last{0.0};
  const uint64_t n_block{1000};
  const uint64_t n_points_2{1000000};
  for (uint64_t i{0}; i != n_points_2; ++i) {
    const double current{i * scale};
    const double diff{current - last};
    last = current;
    test << diff << " " << diff << "\n";
    const uint64_t mod_i{(i + 1) % n_block};
    if (i && mod_i == 0) {
      test << n_block << " {translate circle2 ufill} bind repeat\n";
    } else if (i + 1 == n_points_2) {
      test << mod_i << " {translate circle2 ufill} bind repeat\n";
    }
  }

  test << R"xxx(
0 0 1 c
)xxx";

  last = 0;
  for (uint64_t i{0}; i != n_points_2; ++i) {
    const double current{i * scale};
    const double diff{current - last};
    last = current;
    test << diff << " " << diff << "\n";
    const uint64_t mod_i{(i + 1) % n_block};
    if (i && mod_i == 0) {
      test << n_block << " {translate circle3 execform} bind repeat\n";
    } else if (i + 1 == n_points_2) {
      test << mod_i << " {translate circle3 execform} bind repeat\n";
    }
  }

  old_plot.ps(test.str());

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
