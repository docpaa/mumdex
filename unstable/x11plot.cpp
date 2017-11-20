//
// x11plot
//
// data exploration
//
// Copyright 2016 Peter Andrews @ CSHL
//

#include <exception>
#include <iostream>
#include <fstream>

#include "error.h"
#include "tsv.h"
#include "x11plot.h"

using std::cout;
using std::cerr;
using std::endl;
using std::exception;
using std::ofstream;

using paa::Error;
using paa::TSV;
using paa::X11Plotter;

int main(int argc, char** argv) try {
  if (--argc != 1) throw Error("usage: x11plot input_file");

  TSV tsv{argv[1]};

  cerr << "create plotter" << endl;
  // Create an X11 plotter from TSV
  X11Plotter plotter{tsv};

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
