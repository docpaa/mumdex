//
// qsd_enzymes.cpp
//
// Design enzyme cutting experiment for QSD
//   Not completed
//
// Copyright 2019 Peter Andrews
//

#include <fstream>
#include <future>
#include <functional>
#include <iostream>
#include <memory>
#include <random>
#include <string>
#include <utility>
#include <vector>
#include <climits>

#include "error.h"
#include "threads.h"
#include "utility.h"

using std::bind;
using std::cerr;
using std::cin;
using std::cout;
using std::endl;
using std::exception;
using std::function;
using std::future;
using std::ifstream;
using std::istream;
using std::make_unique;
using std::mt19937_64;
using std::random_device;
using std::ref;
using std::string;
using std::uniform_int_distribution;
using std::unique_ptr;
using std::vector;

using paa::Error;
using paa::ThreadPool;

int main(int argc, char ** argv) try {
  paa::exit_on_pipe_close();
  std::ios_base::sync_with_stdio(false);
  cin.tie(nullptr);

  // read command line arguments
  const string usage{"usage: qsd_enzymes enzyme_list flanking_length"};
  if (--argc != 2) throw Error(usage);
  const string enzymes_file_name{argv[1]};
  const unsigned int flanking_length{static_cast<unsigned int>(atoi(argv[2]))};

  // read restriction enzyme sequences
  ifstream enzymes_file{enzymes_file_name.c_str()};
  if (!enzymes_file)
    throw Error("Cound not open enzymes file") << enzymes_file_name;
  vector<string> enzyme_names;
  vector<string> enzyme_sequences;
  string enzyme_name;
  string enzyme_sequence;
  while (enzymes_file >> enzyme_name >> enzyme_sequence) {
    enzyme_names.push_back(enzyme_name);
    enzyme_sequences.push_back(enzyme_sequence);
  }
  cerr << "Read " << enzyme_names.size() << " enzymes" << endl;

  const string bases{"ACGT"};
  const string ns{"NNNN"};

  // one trial result and function
  struct TrialResult {
    explicit TrialResult(const vector<string> & enzyme_sequences_) :
        exact_matches(enzyme_sequences_.size()) {}
    string sequence{};
    vector<unsigned int> exact_matches;
    void get_matches(const vector<string> &) {
    }
  };
  auto trial = [flanking_length, &enzyme_sequences, &ns, &bases]() {
    static thread_local mt19937_64 mersenne{random_device()()};
    static thread_local uniform_int_distribution<unsigned char> dist{
      0, static_cast<unsigned char>(bases.size() - 1)};
    static thread_local function<unsigned char()> gen{bind(dist, mersenne)};
    TrialResult result{enzyme_sequences};
    for (unsigned int b{0}; b != flanking_length; ++b)
      result.sequence.push_back(bases[gen()]);
    result.sequence += ns;
    for (unsigned int b{0}; b != flanking_length; ++b)
      result.sequence.push_back(bases[gen()]);

    return result;
  };

  // launch threads
  const unsigned int n_threads{192};
  ThreadPool pool{n_threads};
  const unsigned int n_trials{100};
  vector<future<TrialResult>> futures;
  for (uint64_t t{0}; t != n_trials; ++t) futures.push_back(pool.run(trial));

  // collect results
  for (future<TrialResult> & future : futures) {
    const TrialResult trial_result{future.get()};
    const string & sequence{trial_result.sequence};
    cerr << sequence << endl;
  }

  return 0;
} catch (Error & e) {
  cerr << "paa::Error:" << endl;
  cerr << e.what() << endl;
  return 1;
} catch(exception & e) {
  cerr << e.what() << endl;
  return 1;
} catch(...) {
  cerr << "Some exception was caught." << endl;
  return 1;
}
