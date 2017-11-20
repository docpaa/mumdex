//
// mumdex_cn.cpp
//
// Copy number from a mumdex file
//
// Copyright 2016 Peter Andrews @ CSHL
//

#include <algorithm>
#include <exception>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "cn.h"
#include "error.h"
#include "files.h"
#include "mumdex.h"
#include "threads.h"
#include "utility.h"

using std::cout;
using std::cerr;
using std::endl;
using std::exception;
using std::istringstream;
using std::max;
using std::min;
using std::string;
using std::to_string;
using std::vector;

using paa::Bin;
using paa::CN_abspos;
using paa::Error;
using paa::Mappability;
using paa::MappedVector;
using paa::MUMdex;
using paa::PosInfo;
using paa::Progress;
using paa::Reference;
using paa::ThreadPool;

int main(int argc, char * argv[]) try {
  paa::exit_on_pipe_close();

  --argc;
  if (argc != 4 && argc != 5) throw Error(
          "usage: mumdex_cn mumdex bins bad sample [title]");

  paa::set_cn_parameters();

  const bool create_pdf{false};
  cerr << "Loading mumdex" << endl;
  const MUMdex mumdex{argv[1]};
  const string bins_string{argv[2]};
  const MappedVector<unsigned char> bad{argv[3]};
  const unsigned int sample{static_cast<unsigned int>(atoi(argv[4]))};
  if (!sample) throw Error("Sample must be a nonzero integer");
  const string title{argc == 5 ? argv[5] : ""};

  cerr << "Loading reference and mappability" << endl;
  const Reference & ref{mumdex.reference()};
  const CN_abspos cn_abspos{ref};
  const Mappability mappability{ref};

  cerr << "Loading bins" << endl;
  using AllBins = vector<vector<Bin> >;
  const AllBins bins{[&bins_string, &ref]() {
      vector<string> bins_names;
      istringstream bins_stream{bins_string.c_str()};
      string bins_name;
      while (getline(bins_stream, bins_name, ',')) {
        bins_names.push_back(bins_name);
      }
      AllBins result(bins_names.size());
      for (unsigned int n{0}; n != bins_names.size(); ++n) {
        result[n] = load_bins(bins_names[n], ref);
      }
      return result;
    }()};

  cerr << "Filtering mappings from " << mumdex.n_pairs()
       << " pairs and assigning to bins" << endl;

  const unsigned int n_threads{12};
  ThreadPool pool{n_threads};
  ThreadPool::Results<vector<vector<unsigned int>>> results;
  const uint64_t block_size{max(static_cast<uint64_t>(1000),
                                mumdex.n_pairs() / n_threads / 10)};

  for (uint64_t b_{0}; b_ < mumdex.n_pairs(); b_ += block_size) {
    auto block_fun = [&mumdex, &mappability, &cn_abspos, &bad, &bins,
                      block_size, sample] (const uint64_t b) {
      try {
        vector<vector<unsigned int>> counts(bins.size());
        for (unsigned int i{0}; i != bins.size(); ++i) {
          counts[i].resize(bins[i].size());
        }
        for (uint64_t p{b}; p < min(b + block_size, mumdex.n_pairs());
             p += sample) {
          const unsigned int abspos{pair_cn_abspos(mumdex, mappability,
                                                   cn_abspos, p)};
          if (abspos >= cn_abspos.bad_abspos()) continue;
          if (bad[abspos]) continue;
          for (unsigned int i{0}; i != bins.size(); ++i) {
            vector<Bin>::const_iterator found{
              upper_bound(bins[i].begin(), bins[i].end(), abspos)};
            if (found-- != bins[i].begin()) {
              if (abspos >= found->abspos_start() &&
                  abspos < found->abspos_start() + found->length()) {
                ++counts[i][found - bins[i].begin()];
              }
            }
          }
        }
        return counts;
      } catch (Error & e) {
        cerr << "paa::Error:" << endl;
        cerr << e.what() << endl;
        throw;
      } catch (exception & e) {
        cerr << "std::exception" << endl;
        cerr << e.what() << endl;
        throw;
      } catch (...) {
        cerr << "unknown exception was caught" << endl;
        throw;
      }
    };
    pool.run(results, block_fun, b_);
  }

  Progress progress{results.size(), 0.01, "Binning"};
  vector<vector<unsigned int>> counts(bins.size());
  for (unsigned int i{0}; i != bins.size(); ++i) {
    counts[i].resize(bins[i].size());
  }
  while (results.size()) {
    const vector<vector<unsigned int>> block_counts{results.get()};
    for (unsigned int i{0}; i != bins.size(); ++i) {
      for (unsigned int b{0}; b != counts[i].size(); ++b) {
        counts[i][b] += block_counts[i][b];
      }
    }
    progress();
  }

  for (unsigned int i{0}; i != bins.size(); ++i) {
    copy_number(ref, bins[i], counts[i],
                title + " " + to_string(bins[i].size()) + " bins",
                create_pdf);
  }

  cerr << "Done" << endl;

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
