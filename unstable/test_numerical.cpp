//
// test_numerical
//
// test the numerical.h header
//
// Copyright Peter Andrews 2017 CSHL
//

#include <cmath>
#include <exception>
#include <iostream>

#include "error.h"
#include "numerical.h"

using std::cout;
using std::cerr;
using std::endl;
using std::exception;
using std::fabs;
using std::sqrt;

using paa::Error;
using paa::GoldenMinimizer;

double test_fun(const double x) {
  return (x - 5.000000001) * (x - 5.000000001) + 10;
}

int main(int argc, char * argv[]) try {
  if (--argc != 0) throw Error("usage: test_numerical");
  if (0) cout << argv[0] << endl;
  const GoldenMinimizer minimizer(test_fun, -1000, 1000);
  cout << "Minimum value is at " << minimizer.min() << endl;

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
