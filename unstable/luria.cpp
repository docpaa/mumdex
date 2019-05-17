//
// luria
//
// play with luria-delbruck mutational process
//
// Copyright 2018 Peter Andrews @ CSHL
//

#include <exception>
#include <functional>
#include <iostream>
#include <fstream>
#include <limits>
#include <random>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "error.h"
#include "files.h"
#include "psplot.h"
#include "pstream.h"
#include "threads.h"
#include "utility.h"

using std::bernoulli_distribution;
using std::binomial_distribution;
using std::cout;
using std::cerr;
using std::endl;
using std::exception;
using std::function;
using std::ifstream;
using std::mt19937_64;
using std::numeric_limits;
using std::ofstream;
using std::ostringstream;
using std::pair;
using std::random_device;
using std::ref;
using std::string;
using std::to_string;
using std::uniform_int_distribution;
using std::vector;

using paa::commas;
using paa::get_next_file;
using paa::readable;
using paa::Bounds;
using paa::Error;
using paa::Marker;
using paa::PSDoc;
using paa::PSGraph;
using paa::PSXYSeries;
using paa::ThreadPool;

using redi::ipstream;

int main(int argc, char * argv[]) try {
  if (--argc < 3)
    throw Error("usage: luria [-(m|s)] n_generations n_trials n_threads p ...");

  const bool master{string(argv[1]) == "-m"};
  const bool slave{string(argv[1]) == "-s"};
  if (master || slave) {
    cerr << "Running in " << (master ? "master" : "slave") << " mode" << endl;
    --argc;
    ++argv;
  }

  const uint64_t n_generations{strtoul(argv[1], nullptr, 10)};
  const uint64_t n_trials{strtoul(argv[2], nullptr, 10)};
  const unsigned int n_threads{static_cast<unsigned int>(
      strtoul(argv[3], nullptr, 10))};
  argc -= 3;
  argv += 3;

  constexpr uint64_t n_bins{10001};
  const uint64_t trials_per_job{(n_bins - 1) * 10};

  PSDoc plots{"luria"};
  plots.pdf(true);
  ThreadPool pool{n_threads};

  while (argc--) {
    const double p{strtod((argv++)[1], nullptr)};
    auto bin_bin = [] (const double value) {
      if (value > 100) throw Error("Bad percentage > 100");
      if (value < 0) throw Error("Bad percentage < 0");
      return (n_bins - 1) * value / 100;
    };
    auto bin_value = [] (const uint64_t bin) {
      return 100 * (bin + 0.5) / (n_bins - 1);
    };

    using Result = vector<uint64_t>;  // pair<uint64_t, uint64_t>;
    auto do_trial = [bin_bin, n_generations, p]
        (const unsigned int seed, const uint64_t trials_this_job) {
      mt19937_64 mersenne{seed};
      vector<uint64_t> bin_data(n_bins);
      for (uint64_t t{0}; t != trials_this_job; ++t) {
        uint64_t n_mutations{0};
        uint64_t generation{0};
        uint64_t n_cells{1};
        while (++generation <= n_generations) {
          n_mutations +=
              binomial_distribution<uint64_t>{n_cells - n_mutations, p}(
                  mersenne);
          n_mutations -=
              binomial_distribution<uint64_t>{n_mutations, p}(mersenne);
          if (generation != n_generations) {
            n_cells *= 2;
            n_mutations *= 2;
          }
        }
        const double percent_mutation{100.0 * n_mutations / n_cells};
        ++bin_data[bin_bin(percent_mutation)];
      }
      return bin_data;
    };

    // output dir
    random_device rd;
    mt19937_64 mersenne{rd()};
    uniform_int_distribution<uint64_t> uniform{
      0, numeric_limits<uint64_t>::max()};
    ostringstream hist_dir;
    hist_dir << "hists." << n_generations << "." << p << "." << n_trials
             << "." << n_bins;
    if (system(("mkdir -p " + hist_dir.str()).c_str()) !=0)
      throw Error("Could not create directory") << hist_dir.str();

    Result bin_data(n_bins);

    // load old results
    uint64_t total_trials{n_trials};
    if (master) {
      ipstream file_stream{"find " + hist_dir.str() + " -name '*.counts.txt'"};
      string fn;
      while (file_stream >> fn) {
        // if (fn == hist_name) continue;
        if (readable(fn)) {
          ifstream old_data{fn.c_str()};
          if (!old_data) throw Error("problem reading old hist") << fn;
          uint64_t count{0};
          for (uint64_t bin{0}; bin != bin_data.size(); ++bin) {
            old_data >> count;
            bin_data[bin] += count;
          }
          if (!old_data) throw Error("parse error for") << fn;
        total_trials += n_trials;
        if (total_trials == 1000000000000) break;
        } else {
          break;
        }
      }
    } else {
      // run jobs
      ThreadPool::Results<Result> results;
      for (uint64_t trial{0}; trial < n_trials; trial += trials_per_job)
        pool.run(results, do_trial, rd(), trial + trials_per_job >= n_trials ?
                 n_trials - trial : trials_per_job);

      // combine results
      uint64_t trial{0};
      while (results.size()) {
        const Result result{results.get()};
        for (uint64_t b{0}; b != bin_data.size(); ++b)
          bin_data[b] += result[b];
        trial += trials_per_job;
        cerr << results.size() << " " << trial << endl;
      }

      // counts output
      const string hist_name{hist_dir.str() + "/hist." +
            to_string(uniform(mersenne)) + ".counts.txt"};
      if (system(("touch " + hist_name).c_str()) != 0)
        throw Error("Could not touch hist") << hist_name;
      ofstream hist_file{hist_name.c_str()};
      if (!hist_file) throw Error("Problem opening hist file") << hist_name;
      for (uint64_t bin{0}; bin != bin_data.size(); ++bin)
        hist_file << bin_data[bin] << "\n";
    }

    // plot
    if (!slave) {
      ostringstream title;
      title << "Luria-Delbruck distribution simulation for p = " << p
            << ", Ngen = " << n_generations
            << ", Ntrials = " << commas(total_trials);
      PSGraph * graph{PSGraph::create(
          plots, title.str() + ";Percent mutation;N")};
      graph->log_y(true);
      const Marker marker{paa::circle(), 0.1, "0 0 0", 0.1, true, "0 0 0"};
      PSXYSeries * hist{PSXYSeries::create(*graph, marker)};
      for (uint64_t bin{0}; bin != bin_data.size(); ++bin)
        hist->add_point(bin_value(bin), bin_data[bin]);
    }
  }

  cerr << "done" << endl;

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
