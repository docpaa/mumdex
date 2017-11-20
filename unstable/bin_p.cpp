//
// bin_p
//
// approximate binomial p value
//
// Copyright 2016 Peter Andrews @ CSHL
//

#include <exception>
#include <iostream>
#include <random>

#include "error.h"

using std::cout;
using std::cerr;
using std::endl;
using std::exception;

using paa::Error;

int main(int argc, char* argv[], char * []) try {
  if (--argc != 2) throw Error("usage: bin_p trials successes");

  const uint64_t n_trials{static_cast<uint64_t>(atol(argv[1]))};
  const uint64_t n_successes_read{static_cast<uint64_t>(atol(argv[2]))};
  const uint64_t n_successes{n_successes_read * 2 > n_trials ?
      n_successes_read : n_trials - n_successes_read};
  const uint64_t n_repeats{10000};

  auto mersenne = std::mt19937_64();
  mersenne.seed(time(nullptr));
  std::binomial_distribution<uint64_t> dist{n_trials, 0.5};

  uint64_t n_sim_trials{0};
  uint64_t n_sim_successes{0};
  while (n_sim_successes < n_repeats) {
    ++n_sim_trials;
    if (dist(mersenne) >= n_successes) {
      ++n_sim_successes;
    }
  }

  cout << "P value = " << 1.0 * n_sim_successes / n_sim_trials << endl;

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
