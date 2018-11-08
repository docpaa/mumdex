//
// haha
//
// Process HAHA data to find genome phase
//
// Copyright 2018 Peter Andrews @ CSHL
//

#include <exception>
#include <iostream>
#include <fstream>
#include <functional>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "error.h"
#include "files.h"
#include "haha.h"
#include "mumdex.h"
#include "psplot.h"
#include "threads.h"

using std::cout;
using std::cerr;
using std::endl;
using std::exception;
using std::ofstream;
using std::ostringstream;
using std::set;
using std::stoul;
using std::string;
using std::vector;

using paa::mkdir;
using paa::ChromosomeIndexLookup;
using paa::CoverageHMM;
using paa::Error;
using paa::FragmentInfo;
using paa::FragmentTable;
using paa::HaHaHMM;
using paa::HetTable;
using paa::PSDoc;
using paa::Reference;
using paa::ThreadPool;

int main(int argc, char * argv[]) try {
  const std::string usage{R"xxx(usage:
     haha sample ref in_dir out_dir n_samples n_versions
     haha copy ref in_dir out_dir init0 bkg flen n_threads chr ...
     haha haha ref in_dir out_dir init0 bkg flen n_threads chr
)xxx"};

  const string analysis_type{argv[1]};
  --argc;
  if (!(analysis_type == "sample" && argc == 6) &&
      !(analysis_type == "copy" && argc >= 9) &&
      !(analysis_type == "haha" && argc == 9)) throw Error(usage);

  const Reference ref{argv[2]};
  const ChromosomeIndexLookup chr_lookup{ref};
  const string data_dir{argv[3]};
  const string output_prefix{argv[4]};
  argv += 4; argc -= 4;

  // Chromosomes to process
  const vector<string> chromosome_names{
    [&ref, &chr_lookup, &analysis_type, argc, argv]() {
      vector<string> result;
      if (analysis_type == "sample") {
        for (unsigned int c{0}; c != ref.n_chromosomes(); ++c)
          if (ref.size(c) > 40000000 &&
              ref.name(c).find('X') == string::npos &&
              ref.name(c).find('Y') == string::npos)
            result.push_back(ref.name(c));
      } else {
        for (int a{4}; a != argc; ++a) {
          const std::string name{argv[a + 1]};
          chr_lookup[name];  // throws if name not found
          result.push_back(name);
        }
      }
      return result;
    }()};
  std::cerr << "Running " << analysis_type << " on chromosome"
            << (chromosome_names.size() > 1 ? "s" : "");
  for (const string name : chromosome_names) std::cerr << " " << name;
  std::cerr << endl;

  // Fragment input data
  const vector<FragmentTable> fragment_tables{
    [&chromosome_names, &chr_lookup, &data_dir]() {
      vector<FragmentTable> result;
      result.reserve(chromosome_names.size());
      for (const string name : chromosome_names) {
        const string table_name{data_dir + "/coverage_" + name + ".txt"};
        result.emplace_back(table_name, chr_lookup);
      }
      return result;
    }()};

  // Het input data
  const vector<HetTable> het_tables{
    [&chromosome_names, &chr_lookup, &analysis_type, &data_dir]() {
      vector<HetTable> result;
      if (analysis_type == "copy") return result;  // not needed for "copy"
      result.reserve(chromosome_names.size());
      for (const string name : chromosome_names) {
        const string table_name{data_dir + "/hets_" + name + ".txt"};
        result.emplace_back(table_name, chr_lookup);
      }
      return result;
    }()};

  if (analysis_type == "sample") {
    const uint64_t n_samples{stoul(argv[1])};
    const uint64_t n_versions{stoul(argv[2])};
    for (unsigned int v{0}; v != n_versions; ++v) {
      cerr << "Version " << v + 1 << endl;
      const std::set<uint64_t> samples{[&fragment_tables, n_samples]() {
          std::random_device rd{};
          std::mt19937_64 mersenne{rd()};
          std::function<uint64_t()> gen{std::bind(
              std::uniform_int_distribution<uint64_t>(
                  0, fragment_tables.front().n_samples() - 1), mersenne)};
          std::set<uint64_t> result;
          while (result.size() != n_samples) result.insert(gen());
          return result;
        }()};
      ostringstream out_dir;
      out_dir << output_prefix << "." << n_samples << "." << v;
      mkdir(out_dir.str());
      for (unsigned int c{0}; c != chromosome_names.size(); ++c) {
        const string & chr{chromosome_names[c]};
        const string cover_name{out_dir.str() + "/coverage_" + chr + ".txt"};
        ofstream cover_out{cover_name.c_str()};
        if (!cover_out)
          throw Error("Problem opening coverage out file") << cover_name;
        for (uint64_t f{0}; f != fragment_tables[c].n_fragments(); ++f) {
          cover_out << chr << '\t' << f
                    << '\t' << fragment_tables[c][f].start()
                    << '\t' << fragment_tables[c][f].stop();
          for (const uint64_t s : samples)
            cover_out << '\t' << paa::counts[fragment_tables[c](s, f)];
          cover_out << '\n';
        }
        const string hets_name{out_dir.str() + "/hets_" + chr + ".txt"};
        ofstream hets_out{hets_name.c_str()};
        if (!hets_out)
          throw Error("Problem opening hets out file") << hets_name;
        for (uint64_t h{0}; h != het_tables[c].n_hets(); ++h) {
          hets_out << chr << '\t' << het_tables[c][h].position()
                   << '\t' << h
                   << '\t' << het_tables[c][h].phase1()
                   << '\t' << het_tables[c][h].phase2();
          for (const uint64_t s : samples)
            hets_out << '\t'
                     << paa::alleles_symbols[het_tables[c].unflipped(s, h)];
          hets_out << '\n';
        }
        cover_out.close();
        hets_out.close();
        const HetTable hets{hets_name, chr_lookup};
        const FragmentTable fragments{cover_name, chr_lookup};
      }
    }
    return 0;
  }

  const double init0{atof(argv[1])};
  const double background{atof(argv[2])};
  const double fragment_length{atof(argv[3])};
  const unsigned int n_threads{static_cast<unsigned int>(stoul(argv[4]))};

  ThreadPool pool{n_threads};
  const string plots_name{data_dir + "/" + analysis_type + "." +
        (chromosome_names.size() == 1 ? chromosome_names.front() : "all")};
  PSDoc plots{plots_name, plots_name};

  // Coverage HMMs for each chromosome
  const vector<CoverageHMM> coverage_hmms{
    [&chromosome_names, &fragment_tables, &pool, &plots,
     init0, background, fragment_length, &analysis_type, &data_dir]() {
      vector<CoverageHMM> coverages;
      coverages.reserve(chromosome_names.size());
      for (unsigned int c{0}; c != chromosome_names.size(); ++c) {
        const string name{data_dir + "/" + chromosome_names[c]};
        coverages.emplace_back(name, fragment_tables[c], pool, plots,
                            init0, background, fragment_length);
      }

      if (analysis_type == "copy") {
        // Run coverage HMM on all chromosomes, if needed
        const unsigned int max_iterations{100};
        unsigned int iteration{0};
        bool converged{false};
        while (!converged && ++iteration < max_iterations) {
          converged = true;
          CoverageHMM::update_counts all_counts{};
          for (unsigned int c{0}; c != chromosome_names.size(); ++c) {
            CoverageHMM & hmm{coverages[c]};
            hmm.viterbi();
            hmm.show_iteration();
            const CoverageHMM::update_counts chr_counts{
              hmm.update_parameters()};
            all_counts += chr_counts;
            if (!hmm.converged) converged = false;
          }
          for (unsigned int c{0}; c != chromosome_names.size(); ++c)
            coverages[c].recalculate_parameters(all_counts);
        }
        for (unsigned int c{0}; c != chromosome_names.size(); ++c) {
          coverages[c].plot_states();
          const string name{chromosome_names[c]};
          ofstream all_out{(name + ".all.txt").c_str()};
          if (!all_out) throw Error("problem opening all.txt");
          coverages[c].show_all(all_out);
        }
        cout << "Coverage HMM final parameters:"
             << " " << coverages[0].init0
             << " " << coverages[0].background
             << " " << coverages[0].fragment_length
             << endl;
      }
      return coverages;
    }()};

  if (analysis_type == "copy") return 0;

  // Initial phase
  vector<double> phase{simple_push(het_tables.front())};

  // HaHa HMM - just for one chromosome
  for (unsigned int i{0}; i != 3 ; ++i) {
    const HaHaHMM haha_hmm{
      [&chromosome_names, &het_tables, &coverage_hmms,
       &pool, &plots, &data_dir, &phase]() {
        const string name{data_dir + "/" + chromosome_names.front()};
        return HaHaHMM{name, pool, plots,
              het_tables.front(), coverage_hmms.front(), phase};
      }()};

    // HaHa afterburner run to fix phase at points of uncertainty
    const double tolerance{-100};
    phase = haha_hmm.afterburner(chromosome_names.front(), ref,
                                 chr_lookup, tolerance);
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
