//
// haplothread_sharing.cpp
//
// summarize haplothread sharing data from many chromosomes
//
// Copyright 2019 Peter Andrews @ CSHL
//

#include <exception>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <map>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "error.h"
#include "haha.h"
#include "psplot.h"
#include "pstream.h"
#include "stats.h"
#include "utility.h"

using std::cout;
using std::cerr;
using std::endl;
using std::exception;
using std::ifstream;
using std::istream;
using std::map;
using std::ofstream;
using std::ostream;
using std::ostringstream;
using std::setprecision;
using std::setw;
using std::string;
using std::vector;

using paa::sround;
using paa::Error;
using paa::Marker;
using paa::Matrix;
using paa::NormalParams;
using paa::PSDoc;
using paa::PSGraph;
using paa::PSXYSeries;
using paa::PSXYMSeries;

class IvanMetadata {
 public:
  struct Result {
    string role;
    string sex;
  };
  explicit IvanMetadata(const string & file_name) {
    ifstream file{file_name.c_str()};
    if (!file) throw Error("Problem opening lookup file") << file_name;
    string family;
    string individual;
    string sex_;
    string role_;
    file.ignore(10000, '\n');
    while (file >> family >> individual >> sex_ >> role_) {
      if (lookup.find(individual) != lookup.end())
        throw Error("Duplicate entry for individual") << individual;
      lookup[individual] = role_;
      sex_lookup[individual] = sex_;
    }
  }
  Result operator[](const string & individual) const {
    try {
      return Result{role(individual), sex(individual)};
    } catch (...) {
      cerr << "[] lookup problem for " << individual << endl;
      throw;
    }
  }
  string role(const string & individual) const {
    try {
      return lookup.at(individual);
    } catch (...) {
      cerr << "role lookup problem for " << individual << endl;
      throw;
    }
  }
  string sex(const string & individual) const {
    try {
      return sex_lookup.at(individual);
    } catch (...) {
      cerr << "sex lookup problem for " << individual << endl;
      throw;
    }
  }

 private:
  map<string, string> lookup{};
  map<string, string> sex_lookup{};
};

// 2 dimensional child_odd x father_odd
class Counts {
 public:
  using CMatrix = Matrix<uint64_t, 2, 2>;
  uint64_t size() const { return counts.size(); }
  const CMatrix & operator[](const uint64_t share_parity) const {
    return counts[share_parity];
  }
  CMatrix & operator[](const uint64_t share_parity) {
    return counts[share_parity];
  }
  Counts & operator+=(const Counts & other) {
    for (const bool share_parity : {false, true})
      for (const bool child_odd : {false, true})
        for (const bool father_odd : {false, true})
          counts[share_parity][child_odd][father_odd] +=
              other.counts[share_parity][child_odd][father_odd];
    if (family.empty()) {
      family = other.family;
      individual1 = other.individual1;
      individual2 = other.individual2;
      is_autistic1 = other.is_autistic1;
      is_autistic2 = other.is_autistic2;
      is_male1 = other.is_male1;
      is_male2 = other.is_male2;
    } else {
      if (family != other.family)
        throw Error("Family mismatch") << id() << other.id();
      if (individual1 != other.individual1)
        throw Error("Individual mismatch") << id() << other.id();
      if (individual2 != other.individual2)
        throw Error("Individual mismatch") << id() << other.id();
      if (is_autistic1 != other.is_autistic1)
        throw Error("Is_Autistic mismatch") << id() << other.id();
      if (is_autistic2 != other.is_autistic2)
        throw Error("Is_Autistic mismatch") << id() << other.id();
      if (is_male1 != other.is_male1)
        throw Error("Is_Male mismatch") << id() << other.id();
      if (is_male2 != other.is_male2)
        throw Error("Is_Male mismatch") << id() << other.id();
    }
    return *this;
  }

  string id() const {
    return family +
        " " + individual1 + " " +  individual2 +
        " " + (is_autistic1 ? "P" : "S") + " " + (is_autistic2 ? "P" : "S");
  }

  istream & input(istream & in) {
    string type1;
    string type2;
    in >> family >> individual1 >> individual2 >> type1 >> type2;
    for (const bool share_parity : {false, true})
      for (const bool child_odd : {false, true})
        for (const bool father_odd : {false, true})
          in >> counts[share_parity][child_odd][father_odd];
    static const IvanMetadata lookup{"sids.txt"};
    is_autistic1 = lookup.role(individual1) == "prb";
    is_male1 = lookup.sex(individual1) == "M";
    is_autistic2 = lookup.role(individual2) == "prb";
    is_male2 = lookup.sex(individual2) == "M";
    return in;
  }

  ostream & output(ostream & out) const {
    const char space{'\t'};
    out << family << space << individual1 << space << individual2
        << space << is_autistic1 << space << is_autistic2
        << space << is_male1 << space << is_male2;
    for (const bool share_parity : {false, true})
      for (const bool child_odd : {false, true})
        for (const bool father_odd : {false, true})
          out << space << counts[share_parity][child_odd][father_odd];
    return out;
  }

  string family{};
  string individual1{};
  string individual2{};
  bool is_autistic1{false};
  bool is_autistic2{false};
  bool is_male1{false};
  bool is_male2{false};

 private:
  vector<CMatrix> counts{2, CMatrix{{0, 0}, {0, 0}}};
};
istream & operator>>(istream & in, Counts & counts) {
  return counts.input(in);
}
ostream & operator<<(ostream & out, const Counts & counts) {
  return counts.output(out);
}

int main(int argc, char * argv[]) try {
  // Process command line arguments
  if (--argc < 2) throw Error("haplothread_summary output_name count_files");
  const string output_name{argv[1]};

  // Read count files for all chromosomes
  vector<vector<Counts>> all_counts;
  string header;
  while (--argc) {
    const string in_name{(++argv)[1]};
    // cerr << "Open " << in_name << endl;
    ifstream in_file{in_name.c_str()};
    if (!in_file) throw Error("Problem opening input file") << in_name;
    all_counts.emplace_back(0);
    Counts counts_line;
    getline(in_file, header);
    while (in_file >> counts_line) all_counts.back().push_back(counts_line);
  }

  // Sum counts from all chromosomes
  vector<Counts> counts(all_counts.front().size());
  for (const vector<Counts> & chr_counts : all_counts)
    for (uint64_t c{0}; c != chr_counts.size(); ++c)
      counts[c] += chr_counts[c];

  // Produce output counts file
  ofstream output{(output_name + ".sharing.txt").c_str()};
  if (!output) throw Error("Could not open output file") << output_name;
  output << "family\tindividual1\tindividual2\tprb1\tprb2\tmale1\tmale2"
      "\tNEM\tNEF\tNOM\tNOF\tSEM\tSEF\tSOM\tSOF\n";
  for (const Counts & ind_counts : counts) output << ind_counts << '\n';

  // Randomize order of individuals to even out plots
  std::random_device rd{};
  std::mt19937_64 mersenne{rd()};
  shuffle(counts.begin(), counts.end(), mersenne);
  cerr << "\nLoaded " << counts.size() << " individuals from "
       << output_name << endl;

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
