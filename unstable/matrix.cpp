//
// matrix
//
// test matrix class
//
// Copyright 2016 Peter Andrews @ CSHL
//

#include "matrix.h"

#include <exception>
#include <iostream>
#include <random>
#include <string>

#include "error.h"
#include "utility.h"

using std::cout;
using std::cerr;
using std::endl;
using std::exception;
using std::string;

using paa::DiagMatrix;
using paa::Error;
using paa::Matrix;
using paa::Progress;

int main(int argc, char* argv[], char * []) try {
  const string usage("usage: matrix gen nx ny nx2 dist [n_threads]"
                     "   or: matrix read l d r dist [n_threads]"
                     "   or: matrix row l r out dummy [n_threads]");
  --argc;
  if (argc != 5 && argc != 6) throw Error(usage);

  const string out_dir{"/home/paa/analysis/matrix"};

  const unsigned int n_threads{argc == 6 ?
        static_cast<unsigned int>(atoi(argv[6])) :
        std::thread::hardware_concurrency()};

  const string type{argv[1]};
  if (type == "gen") {
    const uint64_t nx{static_cast<uint64_t>(atol(argv[2]))};
    const uint64_t ny{static_cast<uint64_t>(atol(argv[3]))};
    const uint64_t nx2{static_cast<uint64_t>(atol(argv[4]))};
    const bool dist{static_cast<bool>(atoi(argv[5]))};

    auto mersenne = std::mt19937_64();
    mersenne.seed(time(nullptr));
    std::uniform_real_distribution<double> real{0, 1};
    auto gen = bind(real, mersenne);

    Matrix<double> l{nx, ny};
    DiagMatrix<double> d{nx};
    Matrix<double> r{nx2, nx, false};
    Progress progress1{nx * ny, 0.1, "Matrix fill 1"};
    for (uint64_t y{0}; y != ny; ++y) {
      for (uint64_t x{0}; x != nx; ++x) {
        l(x, y) = gen();
        progress1();
      }
    }
    for (uint64_t x{0}; x != nx; ++x) {
      d(x) = gen();
    }
    Progress progress2{nx2 * nx, 0.1, "Matrix fill 2"};
    for (uint64_t x{0}; x != nx2; ++x) {
      for (uint64_t y{0}; y != nx; ++y) {
        r(x, y) = gen();
        progress2();
      }
    }
#if 0
    l.write("l");
    d.write("d");
    r.write("r");
#endif
    const time_t start_time{time(nullptr)};
    const Matrix<double> m{dist ?
          l.diag_mul(d, n_threads).dist_mul(r, out_dir) :
          l.diag_mul(d, n_threads).thread_mul(r, n_threads)};
    const time_t stop_time{time(nullptr)};
    cerr << "Elapsed all calculation " << (stop_time - start_time) / 60.0
         << " minutes" << endl;
    m.save("gen_result");
    if (nx <= 10 && ny <= 10 && nx2 <= 10) {
      cout << "Left" << endl;
      cout << l << endl;
      cout << "Diagonal" << endl;
      cout << d << endl;
      cout << "Right" << endl;
      cout << r << endl;
      cout << "Result" << endl;
      cout << m << endl;
    }
  } else if (type == "read") {
    const Matrix<double> l{argv[2]};
    const DiagMatrix<double> d{argv[3]};
    const Matrix<double> r{argv[4], false};
    const bool dist{static_cast<bool>(atoi(argv[5]))};
    const time_t start_time{time(nullptr)};
    const Matrix<double> m{dist ?
          l.diag_mul(d, n_threads).dist_mul(r, out_dir) :
          l.diag_mul(d, n_threads).thread_mul(r, n_threads)};
    const time_t stop_time{time(nullptr)};
    cerr << "Elapsed all calculation " << (stop_time - start_time) / 60.0
         << " minutes" << endl;
    m.save("read_result");
  } else if (type == "row") {
    const Matrix<double> l{argv[2]};
    const Matrix<double> r{argv[3], false};
    const string out{argv[4]};
    l.thread_mul(r, n_threads).write(out);
  } else {
    throw Error("Unknown type") << type << "\n" << usage;
  }

  cerr << "all done" << endl;

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
