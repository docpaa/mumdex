//
// model_transmission.cpp
//
// Explain observed spmf patterns
//
// Copyright 2017 Peter Andrews @ CSHL
//

#include <exception>
#include <iostream>
#include <map>
#include <random>
#include <string>
#include <functional>
#include <vector>

#include "error.h"

using std::bind;
using std::cout;
using std::cerr;
using std::endl;
using std::exception;
using std::map;
using std::mt19937_64;
using std::random_device;
using std::string;
using std::vector;

using paa::Error;

int main(int argc, char* [], char * []) try {
  if (--argc != 0) throw Error("usage: model_transmission");

  const unsigned int n_runs{6200000};
  const double average_coverage{25};
  const double homozygous_rate{0.005};
  const double heterozygous_rate{1 - homozygous_rate};
  const double parent_share_rate{0.0003};
  const double father_rate{0.515};
  const double spoof_rate{0};
  const double denovo_rate{0};
  const unsigned int min_count{10};

  random_device rd;
  auto mersenne = mt19937_64(rd());
  auto realGen = bind(std::uniform_real_distribution<double>(0.0, 1.0),
                      mersenne);
  auto homoGen = bind(std::poisson_distribution<uint64_t>(average_coverage),
                      mersenne);
  auto heteroGen = bind(
      std::poisson_distribution<uint64_t>(average_coverage / 2),
      mersenne);
  const string members{"mfps"};

  map<string, unsigned int> counts;

  for (unsigned int r{0}; r != n_runs; ++r) {
    // const bool heterozygous{realGen() < heterozygous_rate};
    const bool both_parents{realGen() < parent_share_rate};
    const bool in_father{both_parents || realGen() < father_rate};
    const bool in_mother{both_parents || !in_father};
    vector<bool> hetero_status{realGen() < heterozygous_rate,
          realGen() < heterozygous_rate, false, false};
    vector<bool> in_status{in_mother, in_father, false, false};
    for (const unsigned int m : {2, 3}) {
      const bool from_mother{in_status[0] &&
            (!hetero_status[0] || realGen() < 0.5)};
      const bool from_father{in_status[1] &&
            (!hetero_status[1] || realGen() < 0.5)};
      in_status[m] = from_mother || from_father;
      hetero_status[m] = !(from_mother && from_father);
    }
    string out;
    unsigned int max_count{0};
    for (const unsigned int m : {3, 2, 1, 0}) {
      unsigned int count{static_cast<unsigned int>(
          in_status[m] ? (hetero_status[m] ? heteroGen() : homoGen()) : 0)};
      const bool denovo{realGen() < denovo_rate};
      if (count || realGen() < spoof_rate || denovo) {
        out += members[m];
        if (denovo) count = min_count;
        if (count > max_count) max_count = count;
      }
    }
    if (max_count >= min_count) {
      ++counts[out];
    }
  }

  for (const auto & elem : counts) {
    cout <<elem.second << " " << elem.first << endl;
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
