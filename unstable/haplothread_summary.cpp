//
// haplothread_summary.cpp
//
// summarize haplothread data from many chromosomes
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
    string sex;
    string role;
    file.ignore(10000, '\n');
    while (file >> family >> individual >> sex >> role) {
      if (lookup.find(individual) != lookup.end())
        throw Error("Duplicate entry for individual") << individual;
      lookup[individual] = role;
      sex_lookup[individual] = sex;
    }
  }
  Result operator[](const string & individual) const {
    try {
      return Result{role(individual), sex(individual)};
    } catch (...) {
      cerr << "Lookup problem for " << individual << endl;
      throw;
    }
  }
  string role(const string & individual) const {
    try {
      return lookup.at(individual);
    } catch (...) {
      cerr << "Lookup problem for " << individual << endl;
      throw;
    }
  }
  string sex(const string & individual) const {
    try {
      return sex_lookup.at(individual);
    } catch (...) {
      cerr << "Sex lookup problem for " << individual << endl;
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
  unsigned int size() const { return counts.I; }
  const uint64_t * operator[](const uint64_t child_odd) const {
    return counts[child_odd];
  }
  uint64_t * operator[](const uint64_t child_odd) {
    return counts[child_odd];
  }
  Counts & operator+=(const Counts & other) {
    for (const bool child_odd : {false, true})
      for (const bool father_odd : {false, true})
        counts[child_odd][father_odd] += other.counts[child_odd][father_odd];
    n_positions += other.n_positions;
    if (family.empty()) {
      family = other.family;
      individual = other.individual;
      is_autistic = other.is_autistic;
      is_male = other.is_male;
      n_kids = other.n_kids;
      n_auts = other.n_auts;
    } else {
      if (family != other.family)
        throw Error("Family mismatch") << id() << other.id();
      if (individual != other.individual)
        throw Error("Individual mismatch") << id() << other.id();
      if (is_autistic != other.is_autistic)
        throw Error("Is_Autistic mismatch") << id() << other.id();
      if (is_male != other.is_male)
        throw Error("Is_Male mismatch") << id() << other.id();
      if (n_kids != other.n_kids)
        throw Error("Unexpected n_kids in add") << id() << other.id();
      if (n_auts != other.n_auts)
        throw Error("Unexpected n_auts in add") << id() << other.id();
    }
    return *this;
  }

  string id() const {
    return family + " " + individual + " " + (is_autistic ? "P" : "S");
  }

  istream & input(istream & in) {
    in >> family >> individual >> type >> n_kids >> n_auts >> n_positions;
    for (const bool child_odd : {false, true})
      for (const bool father_odd : {false, true})
        in >> counts[child_odd][father_odd];
    if (false &&
        individual.find_first_of("ps") != string::npos && family[0] == '1') {
      is_autistic = individual.find('p') != string::npos;
    } else {
      static const IvanMetadata lookup{"sids.txt"};
      is_autistic = lookup.role(individual) == "prb";
      is_male = lookup.sex(individual) == "M";
    }
    return in;
  }

  ostream & output(ostream & out) const {
    const char space{'\t'};
    out << family << space << individual
        << space << is_autistic << space << is_male
        << space << n_kids << space << n_auts
        << space << n_positions;
    for (const bool child_odd : {false, true})
      for (const bool father_odd : {false, true})
        out << space << counts[child_odd][father_odd];
    return out;
  }

  string family{};
  string individual{};
  bool type{false};
  bool is_autistic{false};
  bool is_male{false};
  unsigned int n_kids{0};
  unsigned int n_auts{0};
  unsigned int n_positions{0};

 private:
  CMatrix counts{{0, 0}, {0, 0}};
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

  // Randomize order of individuals to even out plots
  std::random_device rd{};
  std::mt19937_64 mersenne{rd()};
  shuffle(counts.begin(), counts.end(), mersenne);
  cerr << "\nLoaded " << counts.size() << " individuals from "
       << output_name << endl;

  // Produce output counts file
  ofstream output{(output_name + ".txt").c_str()};
  if (!output) throw Error("Could not open output file") << output_name;
  output << "family\tindividual\tisproband\tismale\tnkids\tnauts\tnpositions"
         << "\tEM\tEF\tOM\tOF\n";
  for (const Counts & ind_counts : counts) output << ind_counts << '\n';

  // Make plots
  PSDoc plots{output_name};
  plots.pdf(true);
  const string acol{"1 0 0"};
  const string scol{"0 0 1"};
  const double scale{0.6};
  const Marker aut_marker{paa::circle(), scale, acol, 0.1, true, acol};
  const Marker sib_marker{paa::circle(), scale, scol, 0.1, true, scol};
  const Marker parent_marker{paa::circle(), scale, scol, 1, true, "1 1 1"};
  const Marker * const markers[2]{&sib_marker, &aut_marker};
  PSGraph movement_graph{plots, "Hetness Inheritance;"
        "Dad Hetness to Inherited Change;Mom Hetness to Inherited Change"};
  PSXYMSeries movement_series{movement_graph};
  PSXYMSeries mom_dad_inherited_graph{plots, output_name +
        " Auts vs Sibs;Dad Inherited Hetness;Mom Inherited Hetness"};
  PSXYSeries mom_dad_graph{plots, output_name +
        " Auts vs Sibs;Dad Hetness;Mom Hetness"};

  // Fill plots
  ostringstream move_ps;
  move_ps << "0 0 0 c 0.2 lw\n";
  for (const Counts & ind_counts : counts) {
    const bool child_aut{ind_counts.is_autistic};
    const uint64_t dad_inherited_oddness{ind_counts[true][true]};
    const uint64_t mom_inherited_oddness{ind_counts[true][false]};
    const uint64_t dad_oddness{dad_inherited_oddness + ind_counts[false][true]};
    const uint64_t mom_oddness{mom_inherited_oddness +
          ind_counts[false][false]};
    mom_dad_inherited_graph.add_point(
        dad_inherited_oddness, mom_inherited_oddness, *markers[child_aut]);
    mom_dad_graph.add_point(dad_oddness, mom_oddness);
    movement_series.add_point(
        dad_inherited_oddness, mom_inherited_oddness, *markers[child_aut]);
    movement_series.add_point(dad_oddness / 2.0, mom_oddness / 2.0,
                              parent_marker);
    move_ps << "np "
            << dad_inherited_oddness << " " << mom_inherited_oddness << " gc m "
            << dad_oddness / 2.0 << " " << mom_oddness / 2.0 << " gc l sp\n";
  }
  movement_graph.pre_ps(move_ps.str());

  // Statistics
  for (const bool father : {false, true}) {
    vector<double> ratios[2];
    for (const Counts & ind_counts : counts) {
      const bool child_aut{ind_counts.is_autistic};
      if (true ||
          (ind_counts[true][father] > 600000 &&
           ind_counts[true][father] < 800000))
        ratios[child_aut].push_back(
            1.0 * ind_counts[true][father] /
            (ind_counts[true][father] + ind_counts[false][father]));
    }
    for (const bool autistic : {false, true}) {
      const NormalParams norm{ratios[autistic]};
      if (ratios[autistic].size()) {
        cout << "n: " << setw(4) << ratios[autistic].size()
             << ", prb: " << autistic
             << ", dad: "<< father
             << ", mean: " << norm.mean
             << " +/- " << sround(norm.stdev, 3)
             << " (" << sround(norm.seom(), 3) << ")"
             << endl;
      }
    }
  }
  return 0;

  cout << endl;
  for (const bool male : {false, true}) {
    for (const bool father : {false, true}) {
      vector<double> ratios[2];
      for (const Counts & ind_counts : counts) {
        const bool child_aut{ind_counts.is_autistic};
        const bool child_male{ind_counts.is_male};
        if (child_male == male)
          ratios[child_aut].push_back(
              1.0 * ind_counts[true][father] /
              (ind_counts[true][father] + ind_counts[false][father]));
      }
      for (const bool autistic : {false, true}) {
        if (ratios[autistic].size()) {
          const NormalParams norm{ratios[autistic]};
          cout << "n: " << setw(4) << ratios[autistic].size()
               << ", prb: " << autistic
               << ", male: " << male
               << ", dad: "<< father
               << ", mean: " << norm.mean
               << " +/- " << sround(norm.stdev, 3)
               << " ( " << sround(norm.seom(), 3) << ")"
               << endl;
        }
      }
    }
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
