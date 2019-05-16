//
// 2a2ds.cpp
//
// Load plink input and compute a2ds measure
//
// Copyright 2019 Peter Andrews @ CSHL
//
// 45678911234567892123456789312345678941234567895123456789612345678971234567898

#include <algorithm>
#include <exception>
#include <iostream>
#include <fstream>
#include <future>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "error.h"
#include "files.h"
#include "mumdex.h"
#include "psplot.h"
#include "threads.h"
#include "utility.h"

using std::cout;
using std::cerr;
using std::endl;
using std::exception;
using std::future;
using std::ifstream;
using std::istream;
using std::istringstream;
using std::map;
using std::move;
using std::ofstream;
using std::ostream;
using std::ostringstream;
using std::string;
using std::to_string;
using std::vector;

using paa::Bounds;
using paa::Error;
using paa::Progress;
using paa::PSDoc;
using paa::PSHSeries;
using paa::Reference;
using paa::ThreadPool;

class a2dsIndividual {
 public:
  class eof {};
  class eol {};
  explicit a2dsIndividual(istream & in) {
    unsigned int sex;
    unsigned int prb;
    if (in.eof()) throw eof();

    in >> family >> id >> dad >> mom >> sex >> prb;
    if (!in) throw eof();
    if (sex == 1) {
      is_male = true;
    } else if (sex == 2) {
      is_male = false;
    } else {
      throw Error("Bad sex code") << sex;
    }
    if (prb == 1) {
      is_proband = false;
    } else if (prb == 2) {
      is_proband = true;
    } else {
      throw Error("Bad proband code") << prb;
    }
    if (mom == "0" && dad == "0") {
      is_child = false;
    } else {
      if (mom != "0" && dad != "0") {
        is_child = true;
      } else {
        throw Error("Unexpected is_child");
      }
    }

    while (in && in.peek() != '\n') {
      int v1, v2;
      in >> v1 >> v2;
      if (!in) throw Error("Parse Error");
      if (v1 == 0 && v2 == 0) {
        data.push_back(0);
      } else if (v1 == 1 && v2 == 1) {
        data.push_back(1);
      } else if (v1 == 1 && v2 == 2) {
        data.push_back(2);
      } else if (v1 == 2 && v2 == 2) {
        data.push_back(3);
      } else {
        throw Error("Unexpected data")
            << v1 << v2 << "at" << data.size() << "for" << id;
      }
    }
    in.get();
    cerr << "Loaded " << id << " with size " << data.size() << endl;
  }

  string family{};
  string id{};
  string dad{};
  string mom{};
  bool is_male{false};
  bool is_proband{false};
  bool is_child{false};
  vector<unsigned char> data{};
};

class a2dsIndividuals {
 public:
  explicit a2dsIndividuals(const unsigned int n_threads,
                           vector<string> file_names) {
    ThreadPool pool{n_threads};
    vector<future<vector<a2dsIndividual>>> futures;
    for (const string & file_name : file_names)
      futures.push_back(pool.run(load, file_name));
    for (future<vector<a2dsIndividual>> & fut : futures)
      for (a2dsIndividual ind : fut.get())
        individuals.push_back(move(ind));
    cerr << "done loading files" << endl;
  }
  static vector<a2dsIndividual> load(const string & file_name) {
    ifstream in{file_name.c_str()};
    if (!in) throw Error("Could not open input file") << file_name;
    vector<a2dsIndividual> inds{};
    while (in) {
      try {
        inds.emplace_back(in);
      } catch (a2dsIndividual::eof) {
        if (!in.eof())
          throw Error("Unexpected input file non-eof error in") << file_name;
      }
    }
    cerr << "Loaded file " << file_name << endl;
    return inds;
  }

  uint64_t size() const { return individuals.size(); }
  const a2dsIndividual & operator[](const uint64_t i) const {
    return individuals[i];
  }
  vector<a2dsIndividual>::const_iterator begin() const {
    return individuals.begin();
  }
  vector<a2dsIndividual>::const_iterator end() const {
    return individuals.end();
  }

 private:
  vector<a2dsIndividual> individuals{};
};

class a2dsLookup {
 public:
  explicit a2dsLookup(const a2dsIndividuals & individuals) {
    for (unsigned int ind{0}; ind != individuals.size(); ++ind) {
      const a2dsIndividual & individual{individuals[ind]};
      auto found = id2index.emplace(individual.id, ind);
      if (found.second != true)
        throw Error("Existing id") << individual.id;
    }
  }
  unsigned int operator[](const string & name) const {
    return id2index.at(name);
  }
  map<string, unsigned int> id2index{};
};

class Individuals {
 public:
  explicit Individuals(const string & file_name,
                       const a2dsLookup & lookup) {
    ifstream in{file_name.c_str()};
    if (!in) throw Error("Problem Reading Individuals") << file_name;
    string name;
    while (in >> name) {
      names.push_back(name);
      ids.push_back(lookup[name]);
    }
    cerr << "Loaded " << names.size() << " names from " << file_name << endl;
  }
  Individuals(const Individuals & cds, const a2dsIndividuals & individuals,
              const bool probands) {
    for (const unsigned int id : cds.ids) {
      const a2dsIndividual & ind{individuals[id]};
      if (ind.is_proband == probands) {
        names.push_back(ind.id);
        ids.push_back(id);
      }
    }
    cerr << "Loaded " << ids.size() << " cds "
         << (probands ? "probands" : "siblings") << endl;
  }
  vector<string> names{};
  vector<unsigned int> ids{};
};

int main(int argc, char * argv[]) try {
  if (--argc < 4) throw Error("a2ds n_threads cds test dataset ...");

  const unsigned int n_threads{static_cast<unsigned int>(atoi(argv[1]))};
  const string cds_name{argv[2]};
  const string test_name{argv[3]};
  argc -= 3;
  argv += 3;
  vector<string> file_names;
  while (argc--) file_names.push_back((argv++)[1]);
  const a2dsIndividuals individuals{n_threads, file_names};
  if (individuals.size() == 0) throw Error("No individuals loaded");
  const a2dsLookup lookup{individuals};
  const Individuals cds{cds_name, lookup};
  const Individuals test{test_name, lookup};
  const Individuals cds_prb{cds, individuals, true};
  const Individuals cds_sib{cds, individuals, false};
  const uint64_t n_positions{individuals[0].data.size()};

  // Get population frequencies from parents
  const vector<double> pop_freqs{[&individuals, n_positions]() {
      vector<uint64_t> allele_seen(n_positions);
      vector<uint64_t> good_pos(n_positions);
      for (const a2dsIndividual & ind : individuals) {
        if (ind.is_child) continue;
        for (uint64_t p{0}; p != n_positions; ++p) {
          if (ind.data[p]) {
            good_pos[p] += 2;
            allele_seen[p] += ind.data[p];
          }
        }
      }
      vector<double> result(n_positions);
      for (uint64_t p{0}; p != n_positions; ++p) {
        if (good_pos[p]) {
          result[p] = 1.0 * allele_seen[p] / good_pos[p];
          if (result[p] > 0.5) result[p] = 1 - result[p];
        }
      }
      return result;
    }()};

  PSDoc plots{"a2ds"};
  PSHSeries<double, uint64_t> pop_freq{plots,
        "Parental allele population frequency;Frequency;N", 100,
        Bounds{0.0, 0.1}};
  for (uint64_t p{0}; p != n_positions; ++p)
    pop_freq.add_point(pop_freqs[p]);

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
