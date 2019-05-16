//
// kenny_positions
//
// Load / process haplothread data
//
// Copyright 2019 Peter Andrews @ CSHL
//

#include <algorithm>
#include <exception>
#include <iostream>
#include <fstream>
#include <functional>
#include <future>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "error.h"
#include "files.h"
#include "ivan_metadata.h"
#include "mumdex.h"
#include "haha.h"
#include "psplot.h"
#include "pstream.h"
#include "threads.h"
#include "utility.h"

using std::cout;
using std::cerr;
using std::endl;
using std::exception;
using std::function;
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

using paa::readable;
using paa::ChromosomeIndexLookup;
using paa::Error;
using paa::IvanMetadata;
using paa::Matrix;
using paa::Progress;
using paa::Reference;
using paa::ThreadPool;

IvanMetadata metadata{"metaData.txt"};

struct PosInfo {
  unsigned int chromosome{0};
  unsigned int pos{0};
  char ref_allele{'0'};
  unsigned int counts[4]{0, 0, 0, 0};
};

const char * const bases{"ACGT"};

const unsigned int max_families{10000};

class AllPosInfo {
 public:
  AllPosInfo(const Reference & ref_,
             const std::string & datadir_,
             const std::string & dataset_,
             const std::string & chr_name_,
             const ChromosomeIndexLookup & chr_lookup) :
      ref{ref_},
      datadir{datadir_},
      dataset{dataset_},
      chr_name{chr_name_} {
        ostringstream input_name;
        input_name << datadir << "/" << dataset;
        if (chr_name != "all") input_name << "-" << chr_name;
        const string pos_name{input_name.str() + "-pos.txt"};
        ifstream pos_file{pos_name.c_str()};
        if (!pos_file) throw Error("Could not open pos file") << pos_name;
        pos_file.ignore(100000, '\n');
        string chrom_name;
        string ref_allele;
        PosInfo info;
        while (pos_file >> chrom_name >> info.pos >> ref_allele
               >> info.counts[0] >> info.counts[1]
               >> info.counts[2] >> info.counts[3]) {
          --info.pos;
          if (chrom_name[0] == 'c') {
            info.chromosome = chr_lookup[chrom_name];
          } else {
            info.chromosome = chr_lookup["chr" + chrom_name];
          }
          info.ref_allele = ref_allele[0];
          data.push_back(info);
          if (info.ref_allele != ref[info.chromosome][info.pos])
            throw Error("Ref allele mismatch")
                << info.ref_allele << ref[info.chromosome][info.pos]
                << pos_name << ref.name(info.chromosome) << info.pos;
        }
        if (data.empty()) throw Error("Empty AllPosInfo");
      }
  const PosInfo & operator[](const uint64_t pos) const { return data[pos]; }
  uint64_t size() const { return data.size(); }

  const Reference & ref;
  std::string datadir;
  std::string dataset;
  std::string chr_name;

 private:
  vector<PosInfo> data{};
};

class AlleleCall {
 public:
  AlleleCall() {}
  explicit AlleleCall(const char data_) : data{data_} {
    switch (data) {
      case 'A':
      case 'C':
      case 'G':
      case 'T':
          break;
      default:
        throw Error("Unknown call character") << data;
    }
  }

 private:
  char data{0};
};

class IndividualInfo {
 public:
  class eof {};
  explicit IndividualInfo(istream & input,
                          function<bool(const string&)> is_aut = ssc_is_aut) {
    string this_family;
    string this_individual;
    string thread_type;
    for (const bool is_father : {false, true}) {
      for (const bool is_transmitted : {true, false}) {
        input >> this_family >> this_individual >> thread_type
              >> data[is_father][is_transmitted];
        if (0) cerr << "Line " << this_family << " " << this_individual << endl;
        if (!input) {
          if (is_father == false && is_transmitted == true) {
            throw eof();
          } else {
            throw Error("Parse Error in IndividualInfo")
                << this_family << this_individual
                << family << individual;
          }
        }
        if (family.size()) {
          if (family != this_family)
            throw Error("Family mismatch") << family << this_family;
        } else {
          family = this_family;
        }
        if (individual.size()) {
          if (individual != this_individual)
            throw Error("Individual mismatch");
        } else {
          individual = this_individual;
          if (is_aut(individual)) is_autistic = true;
        }
      }
    }
  }

  static bool ssc_is_aut(const string & sample_name) {
    return sample_name.find('p') != string::npos;
  }

  char operator()(const uint64_t index,
                  const bool is_father,
                  const bool is_transmitted) const {
    return data[is_father][is_transmitted][index];
  }

  void reduce() {
    for (const bool is_father : {false, true}) {
      for (const bool is_transmitted : {true, false}) {
        data[is_father][is_transmitted].clear();
        data[is_father][is_transmitted].shrink_to_fit();
      }
    }
  }
  uint64_t size() const { return data[0][0].size(); }
  string family{};
  string individual{};
  bool is_autistic{false};
  Matrix<string, 2, 2> data{};
};

// 3 dimensional child_odd x father_odd
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
    for (uint64_t child_odd{0}; child_odd != counts.I; ++child_odd)
      for (uint64_t father_odd{0}; father_odd != counts.J; ++father_odd)
        counts[child_odd][father_odd] += other.counts[child_odd][father_odd];
    return *this;
  }

 private:
  CMatrix counts{{0, 0}, {0, 0}};
};

bool is_good(const char allele) {
  switch (allele) {
    case 'A':
    case 'C':
    case 'G':
    case 'T':
      return true;
    default:
      return false;
  }
}

bool is_odd(const string & genotype) {
  if (genotype.size() != 2) throw Error("Bad genotype") << genotype;
  return genotype[0] != genotype[1];
}

struct KennyPosition {
  string line{};
  unsigned int pos{0};
  char allele1{'0'};
  char allele2{'0'};
};
vector<KennyPosition> kenny_positions;

string genotype_code(const char k1, const char k2,
                     const char i1, const char i2) {
  const string genotype_het{k1 < k2 ? string{k1, k2} : string{k2, k1}};
  const string genotype_1{k1, k1};
  const string genotype_2{k2, k2};
  const string ig{i1 < i2 ? string{i1, i2} : string{i2, i1}};
  if (ig == genotype_het) {
    return "1\t2";
  } else if (ig == genotype_1) {
    return "1\t1";
  } else if (ig == genotype_2) {
    return "2\t2";
  } else {
    return "0\t0";
  }
}

class FamilyInfo {
 public:
  FamilyInfo() {}
  FamilyInfo(const AllPosInfo & pos_info,
             vector<IndividualInfo>::iterator begin,
             vector<IndividualInfo>::iterator end) {
    for (vector<IndividualInfo>::iterator ind{begin}; ind != end; ++ind) {
      if (family.size()) {
        if (ind->family != family)
          throw Error("Family mismatch") << ind->family << family;
      } else {
        family = ind->family;
      }
      if (ind->is_autistic) ++n_auts;
      individuals.push_back(move(*ind));
      counts.emplace_back();
      sharing_counts.emplace_back(vector<vector<Counts>>(
          end - begin, vector<Counts>(2)));
    }

    // process
    for (unsigned int p{0}; p != pos_info.size(); ++p) {
      string last_parental_genotype{""};  // a check
      bool position_good{true};
      string genotype_;
      bool saw_dad_odd{false};
      bool saw_mom_odd{false};
      for (const IndividualInfo & ind : individuals) {
        // good allele
        for (const bool is_dad : {false, true}) {
          for (const bool is_transmitted : {true, false}) {
            const char allele{ind(p, is_dad, is_transmitted)};
            if (!is_good(allele)) {
              position_good = false;
              break;
            }
            genotype_ += allele;
          }
          if (!position_good) break;
        }

        // genotypes
        const string parent_genotypes[2]{
          string() + ind(p, false, true) + ind(p, false, false),
              string() + ind(p, true, true) + ind(p, true, false)};
        const string ind_genotype{
          string() + ind(p, false, true) + ind(p, true, true)};
        const string & mom_genotype{parent_genotypes[false]};
        const string & dad_genotype{parent_genotypes[true]};
        if (0) cerr << ind.family
                    << " " << ind.individual << " [" << ind_genotype
                    << "] [" << mom_genotype << "] [" << dad_genotype << "]"
                    << endl;
        // oddness
        const bool dad_odd{is_odd(dad_genotype)};
        if (dad_odd) saw_dad_odd = true;
        const bool mom_odd{is_odd(mom_genotype)};
        if (mom_odd) saw_mom_odd = true;

        // Only one odd parent
        const unsigned int n_parents_odd{0u + mom_odd + dad_odd};
        if (n_parents_odd != 1) {
          position_good = false;
          break;
        }
      }
      // Check only two alleles in family
      sort(genotype_.begin(), genotype_.end());
      auto uend = unique(genotype_.begin(), genotype_.end());
      if (uend - genotype_.begin() != 2) position_good = false;

      if (!position_good) continue;
      if (saw_mom_odd == saw_dad_odd) throw Error("Unexpected boolean value");
      positions.push_back(p);
      genotype += genotype_ + ';';

      // increment counts!
      for (uint64_t i{0}; i != individuals.size(); ++i) {
        const IndividualInfo & ind{individuals[i]};
        const bool ind_odd{ind(p, false, true) != ind(p, true, true)};
        ++counts[i][ind_odd][saw_dad_odd];
        for (uint64_t i2{i + 1}; i2 != individuals.size(); ++i2) {
          const IndividualInfo & ind2{individuals[i2]};
          const bool ind2_odd{ind2(p, false, true) != ind2(p, true, true)};
          const bool share_parity{ind_odd == ind2_odd};
          ++sharing_counts[i][i2][share_parity][ind_odd][saw_dad_odd];
        }
      }
    }
  }

  void output_kenny(ostream & out, ostream & pos_out,
                    const AllPosInfo & pos_info) const {
    const bool output_pos{family == metadata.first_family()};
    // output genotypes for family for specific positions
    for (uint64_t i{0}; i != individuals.size(); ++i) {
      const IndividualInfo & ind{individuals[i]};
      // Output parent info
      if (i == 0) {
        for (const bool is_dad : {false, true}) {
          const string id{is_dad ? metadata.dad(ind.individual) :
                metadata.mom(ind.individual)};
          out << ind.family
              << '\t' << id
              << '\t' << metadata.dad(id)
              << '\t' << metadata.mom(id)
              << '\t' << (is_dad ? 1 : 2)
              << '\t' << 1;
          unsigned int kp{0};
          for (unsigned int p{0}; p != pos_info.size(); ++p) {
            if (kp == kenny_positions.size()) break;
            if (kenny_positions[kp].pos > pos_info[p].pos) continue;
            if (kenny_positions[kp].pos == pos_info[p].pos) {
              out << '\t' << genotype_code(kenny_positions[kp].allele1,
                                           kenny_positions[kp].allele2,
                                           ind.data[is_dad][true][p],
                                           ind.data[is_dad][false][p]);
            }
            ++kp;
          }
          out << '\n';
        }
      }
      // Output child info
      out << ind.family
          << '\t' << ind.individual
          << '\t' << metadata.dad(ind.individual)
          << '\t' << metadata.mom(ind.individual)
          << '\t' << (metadata[ind.individual].is_male() ? 1 : 2)
          << '\t' << (metadata[ind.individual].is_proband() ? 2 : 1);
      unsigned int kp{0};
      unsigned int n_mismatch{0};
      for (unsigned int p{0}; p != pos_info.size(); ++p) {
        if (kp == kenny_positions.size()) break;
        if (kenny_positions[kp].pos > pos_info[p].pos) continue;
        if (kenny_positions[kp].pos != pos_info[p].pos) {
          ++n_mismatch;
        } else {
          if (output_pos && i == 0) pos_out << kenny_positions[kp].line << '\n';
          out << '\t' << genotype_code(kenny_positions[kp].allele1,
                                       kenny_positions[kp].allele2,
                                       ind.data[false][true][p],
                                       ind.data[true][true][p]);
        }
        ++kp;
      }
      out << '\n';
      if (n_mismatch != 0 && output_pos && i == 0) {
        cerr << "Position mismatches " << n_mismatch
             << " " << 1.0 * n_mismatch / pos_info.size() << endl;
      }
    }
  }

  void reduce() {
    // clear space
    for (auto ind : individuals) ind.reduce();
    genotype.clear();
    genotype.shrink_to_fit();
  }

  static FamilyInfo create(const string & dataset,
                           const string & hpth_name,
                           const string & chr_name,
                           const string & family,
                           const AllPosInfo & pos_info) try {
    const string alt_file_name{"data_new/" + chr_name + "/split." +
          family + ".txt"};
    const string command{readable(alt_file_name) ?
          "cat " + alt_file_name : "tabix -f " + hpth_name + " " + family +
          " 2>&1 | grep -v -e Warning -e tabix -e older"};
    if (0) cerr << command << endl;
    redi::ipstream hpth_file{command.c_str()};
    if (!hpth_file) throw Error("Could not open hpth file for family")
                        << hpth_name << family;
    function<bool(string)> proband_function{};
    proband_function = [](const string & sample_name) {
      return metadata.role(sample_name) == "prb";
    };
    vector<IndividualInfo> individuals;
    while (hpth_file) {
      try {
        individuals.emplace_back(hpth_file, proband_function);
      } catch (IndividualInfo::eof) {
        if (!hpth_file.eof())
          throw Error("Unexpected hpth file non-eof error")
              << dataset << hpth_name << family;
      }
    }

    if (individuals.empty()) throw Error("Empty individual") << family;
    return FamilyInfo{pos_info, individuals.begin(), individuals.end()};
  } catch (Error & e) {
    cerr << "paa::Error:" << endl;
    cerr << e.what() << endl;
    exit(1);
  } catch (exception & e) {
    cerr << "std::exception" << endl;
    cerr << e.what() << endl;
    exit(1);
  } catch (...) {
    cerr << "unknown exception was caught" << endl;
    exit(1);
  }
  uint64_t size() const { return genotype.size(); }

  void output_counts(ostream & out, const char space) const {
    for (uint64_t i{0}; i != individuals.size(); ++i) {
      const IndividualInfo & ind{individuals[i]};
      out << family << space << ind.individual
          << space << (ind.is_autistic ? 1 : 0)
          << space << individuals.size() << space << n_auts
          << space << positions.size();
      for (const bool child_odd : {false, true})
        for (const bool father_odd : {false, true})
          out << space << counts[i][child_odd][father_odd];
      out << '\n';
    }
  }
  void output_sharing(ostream & out, const char space) const {
    for (uint64_t i{0}; i != individuals.size(); ++i) {
      const IndividualInfo & ind{individuals[i]};
      for (uint64_t i2{i+1}; i2 != individuals.size(); ++i2) {
        const IndividualInfo & ind2{individuals[i2]};
        out << family
            << space << ind.individual
            << space << ind2.individual
            << space << (ind.is_autistic ? 1 : 0)
            << space << (ind2.is_autistic ? 1 : 0);
        for (const bool share_parity : {false, true})
          for (const bool child_odd : {false, true})
            for (const bool father_odd : {false, true})
              out << space
                  << sharing_counts[i][i2][share_parity][child_odd][father_odd];
        out << '\n';
      }
    }
  }

  string family{};
  vector<IndividualInfo> individuals{};
  uint64_t n_auts{0};
  vector<unsigned int> positions{};
  string genotype{};
  vector<Counts> counts{};
  vector<vector<vector<Counts>>> sharing_counts{};
};

class AllThreadInfo {
 public:
  explicit AllThreadInfo(const AllPosInfo & pos_info_,
                         ThreadPool & pool) :
      pos_info{pos_info_},
      chr_name{pos_info.chr_name} {
      ostringstream input_name;
      input_name << pos_info.datadir
                 << "/" << pos_info.dataset;
      if (chr_name != "all") input_name << "-" << chr_name;
      const string hpth_name{input_name.str() + "-hpth.txt.gz"};
      const string command{"tabix -l " + hpth_name +
            " 2>&1 | grep -v -e Warning -e tabix -e older | head -n 50000"};
      // cerr << command << endl;
      redi::ipstream family_file{command.c_str()};
      if (!family_file)
        throw Error("Problem opening family file from") << hpth_name;
      family_file.ignore(100000, '\n');

      // Read in families in parallel
      ThreadPool::Results<FamilyInfo> results;
      uint64_t f{0};
      for (string family; family_file >> family;) {
        pool.run(results, FamilyInfo::create, pos_info.dataset, hpth_name,
                 chr_name, family, std::ref(pos_info));
        if (++f == max_families) break;
      }
      Progress progress{results.size(), 1.0 / results.size(), chr_name};
      cerr << "Processing " << results.size() << " results" << endl;
      ofstream out{(pos_info.chr_name + ".txt").c_str()};
      if (!out) throw Error("Problem opening kenny output file");
      ofstream pos_out{(pos_info.chr_name + ".pos.txt").c_str()};
      if (!out) throw Error("Problem opening kenny positions output file");
      for (unsigned int r{0}; results.size(); ++r) {
        FamilyInfo result{results.get()};
        progress();
        cerr << endl;
        result.output_counts(cerr, '\t');
        result.output_kenny(out, pos_out, pos_info);
        result.reduce();
        auto place = data.emplace(result.family, move(result));
        if (place.second == false) throw Error("Duplicate family");
        // place.first->second.add(move(result));
      }
    }

  uint64_t size() const { return data.size(); }

  void output_counts(ostream & out) const {
    const char * const odd{"EO"};
    const char * const parent{"MF"};

    out << "family\tindividual\ttype\tkids\tauts\tpositions";
    auto counts = data.begin()->second.counts;
    for (uint64_t child_odd{0}; child_odd != 2; ++child_odd)
        for (uint64_t father_odd{0}; father_odd != 2; ++father_odd)
          out << '\t' << odd[child_odd] << parent[father_odd];
    out << '\n';
    for (const auto & info : data) info.second.output_counts(out, '\t');
  }
  void output_sharing(ostream & out) const {
    const char * const odd{"EO"};
    const char * const parent{"MF"};
    const char * const parity{"NS"};

    out << "family\tindividual1\tindividual2\ttype1\ttype2";
    auto counts = data.begin()->second.sharing_counts;
    for (const bool share_parity : {false, true})
      for (const bool child_odd : {false, true})
        for (const bool father_odd : {false, true})
          out << '\t'
              << parity[share_parity]
              << odd[child_odd]
              << parent[father_odd];
    out << '\n';
    for (const auto & info : data) info.second.output_sharing(out, '\t');
  }

 private:
  const AllPosInfo & pos_info;
  const std::string chr_name;
  map<string, FamilyInfo> data{};
};

int main(int argc, char * argv[]) try {
  if (--argc != 5) throw Error("kenny_positions ref dir dataset n_threads chr");

  const Reference ref{argv[1]};
  const ChromosomeIndexLookup chr_lookup{ref};
  const string datadir{argv[2]};
  const string dataset{argv[3]};
  const unsigned int n_threads{static_cast<unsigned int>(atoi(argv[4]))};
  const string chromosome_name{argv[5]};

  ThreadPool pool{n_threads};

  const AllPosInfo pos_info{ref, datadir, dataset, chromosome_name, chr_lookup};
  cout << "Loaded " << pos_info.size() << " positions" << endl;

  ifstream kenny{"SPARK.SNP.clean.bim"};
  if (!kenny) throw Error("Problem opening kenny file");
  unsigned int chrno;
  string dummy;
  unsigned int dummy2;
  KennyPosition kp;
  while (getline(kenny, kp.line)) {
    istringstream line_stream{kp.line.c_str()};
    line_stream >> chrno >> dummy >> dummy2
                >> kp.pos >> kp.allele1 >> kp.allele2;
    if (!line_stream) {
      if (line_stream.eof()) break;
      throw Error("Kenny position parse error") << kp.line;
    }
    const string chr_name{chrno < 23 ? "chr" + to_string(chrno) :
          string(chrno == 23 ? "chrX" : (chrno == 24 ? "chrY" : "chrNONE"))};
    --kp.pos;
    if (chr_name == chromosome_name) kenny_positions.push_back(kp);
  }

  const AllThreadInfo haplo_info{pos_info, pool};
  cout << "Loaded " << haplo_info.size() << " families" << endl;

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