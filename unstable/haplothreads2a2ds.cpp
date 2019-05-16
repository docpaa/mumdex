//
// haplothreads2a2ds
//
// Load / process haplothread data to prepare for a2ds processing
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

const IvanMetadata metadata{"metaData.txt"};

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

class AlleleCall {  // unused
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
  explicit IndividualInfo(istream & input) {
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
          is_autistic = metadata[individual].is_proband();
        }
      }
    }
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

struct A2dsPosition {
  string line{};
  unsigned int pos{0};
  char allele1{'0'};
  char allele2{'0'};
};
vector<A2dsPosition> a2ds_positions;

string genotype_code(const char k1, const char k2,
                     char i1, char i2) {
  for (const char allele : {i1, i2}) {
    switch (allele) {
      case 'A':
      case 'C':
      case 'G':
      case 'T':
        break;
      default:
        if (i1 != i2) cerr << "Strange genotype "
                           << k1 << " " << k2 << " "<< i1 << " " << i2 << endl;
        // throw Error("Quitting");
    }
  }
  switch (i1) {
    case 'A':
    case 'C':
    case 'G':
    case 'T':
    case 'c':
    case 'u':
    case 'd':
    case 'o':
        break;
    case 'M':
      i1 = 'A';
      i2 = 'C';
      break;
    case 'R':
      i1 = 'A';
      i2 = 'G';
      break;
    case 'W':
      i1 = 'A';
      i2 = 'T';
      break;
    case 'S':
      i1 = 'C';
      i2 = 'G';
      break;
    case 'Y':
      i1 = 'C';
      i2 = 'T';
      break;
    case 'K':
      i1 = 'G';
      i2 = 'T';
      break;
    default:
      cerr << "Must handle " << i1 << endl;
      exit(1);
  }
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
    // cerr << string{k1, k2, i1, i2} << std::endl;
    return "0\t0";
  }
}

class FamilyInfo {
 public:
  FamilyInfo() {}
  FamilyInfo(vector<IndividualInfo>::iterator begin,
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
    }
  }

  void output_a2ds(ostream & out, ostream & pos_out,
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
            if (kp == a2ds_positions.size()) break;
            if (a2ds_positions[kp].pos > pos_info[p].pos) continue;
            if (a2ds_positions[kp].pos == pos_info[p].pos) {
              out << '\t' << genotype_code(a2ds_positions[kp].allele1,
                                           a2ds_positions[kp].allele2,
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
        if (kp == a2ds_positions.size()) break;
        if (a2ds_positions[kp].pos > pos_info[p].pos) continue;
        if (a2ds_positions[kp].pos != pos_info[p].pos) {
          ++n_mismatch;
        } else {
          if (output_pos && i == 0) pos_out << a2ds_positions[kp].line << '\n';
          out << '\t' << genotype_code(a2ds_positions[kp].allele1,
                                       a2ds_positions[kp].allele2,
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
  }

  static FamilyInfo create(const string & dataset,
                           const string & hpth_name,
                           const string & chr_name,
                           const string & family) try {
    const string alt_file_name{"data_new/data_zipped/" + chr_name + "/split." +
          family + ".txt.gz"};
    const string command{readable(alt_file_name) ?
          "zcat " + alt_file_name : "tabix -f " + hpth_name + " " + family +
          " 2>&1 | grep -v -e Warning -e tabix -e older"};
    if (0) cerr << command << endl;
    redi::ipstream hpth_file{command.c_str()};
    if (!hpth_file) throw Error("Could not open hpth file for family")
                        << hpth_name << family;
    vector<IndividualInfo> individuals;
    while (hpth_file) {
      try {
        individuals.emplace_back(hpth_file);
      } catch (IndividualInfo::eof) {
        if (!hpth_file.eof())
          throw Error("Unexpected hpth file non-eof error")
              << dataset << hpth_name << family;
      }
    }

    if (individuals.empty()) throw Error("Empty individual") << family;
    return FamilyInfo{individuals.begin(), individuals.end()};
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

  string family{};
  vector<IndividualInfo> individuals{};
  uint64_t n_auts{0};
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
      const string command{readable("families.txt") ?
            "cat families.txt" :
            "tabix -l " + hpth_name +
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
                 chr_name, family);
        if (++f == max_families) break;
      }
      Progress progress{results.size(), 10.0 / results.size(),
            pos_info.dataset + " " + chr_name};
      cerr << "Processing " << results.size() << " results" << endl;
      ofstream out{(pos_info.chr_name + ".txt").c_str()};
      if (!out) throw Error("Problem opening a2ds output file");
      ofstream pos_out{(pos_info.chr_name + ".pos.txt").c_str()};
      if (!out) throw Error("Problem opening a2ds positions output file");
      for (unsigned int r{0}; results.size(); ++r) {
        FamilyInfo result{results.get()};
        progress();
        cerr << endl;
        result.output_a2ds(out, pos_out, pos_info);
        result.reduce();
        auto place = data.emplace(result.family, move(result));
        if (place.second == false) throw Error("Duplicate family");
        // place.first->second.add(move(result));
      }
    }

  uint64_t size() const { return data.size(); }

 private:
  const AllPosInfo & pos_info;
  const std::string chr_name;
  map<string, FamilyInfo> data{};
};

int main(int argc, char * argv[]) try {
  if (--argc != 5)
    throw Error("haplothreads2a2ds ref dir dataset n_threads chr");

  const Reference ref{argv[1]};
  const ChromosomeIndexLookup chr_lookup{ref};
  const string datadir{argv[2]};
  const string dataset{argv[3]};
  const unsigned int n_threads{static_cast<unsigned int>(atoi(argv[4]))};
  const string chromosome_name{argv[5]};

  ThreadPool pool{n_threads};

  const AllPosInfo pos_info{ref, datadir, dataset, chromosome_name, chr_lookup};
  cout << "Loaded " << pos_info.size() << " positions" << endl;

  ifstream a2ds{"SPARK.SNP.clean.bim"};
  if (!a2ds) throw Error("Problem opening a2ds file");
  unsigned int chrno;
  string dummy;
  unsigned int dummy2;
  A2dsPosition ap;
  while (getline(a2ds, ap.line)) {
    istringstream line_stream{ap.line.c_str()};
    line_stream >> chrno >> dummy >> dummy2
                >> ap.pos >> ap.allele1 >> ap.allele2;
    if (!line_stream) {
      if (line_stream.eof()) break;
      throw Error("A2ds position parse error") << ap.line;
    }
    const string chr_name{chrno < 23 ? "chr" + to_string(chrno) :
          string(chrno == 23 ? "chrX" : (chrno == 24 ? "chrY" : "chrNONE"))};
    --ap.pos;
    if (chr_name == chromosome_name) a2ds_positions.push_back(ap);
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
