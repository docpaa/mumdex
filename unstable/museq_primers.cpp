//
// museq_primers
//
// find primers for museq
//
// Copyright 2016 Peter Andrews @ CSHL
//

#include <algorithm>
#include <exception>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "error.h"
#include "mumdex.h"
#include "pstream.h"
#include "utility.h"

using std::cerr;
using std::cout;
using std::endl;
using std::exception;
using std::ifstream;
using std::istringstream;
using std::min;
using std::ofstream;
using std::ostringstream;
using std::string;
using std::vector;

using paa::serr;
using paa::sout;
using paa::ChromosomeIndexLookup;
using paa::Error;
using paa::Mappability;
using paa::Reference;

using redi::ipstream;

unsigned int gcat_score(const char base) {
  switch (base) {
    case 'A':
    case 'T':
        return 2;
    case 'C':
    case 'G':
        return 4;
    case 'N':
      return 10000;
    default:
      throw Error("Unknown base seen") << base;
  }
}

int main(int argc, char* argv[])  try {
  if (--argc != 2) throw Error("usage: museq_primers ref_fasta input");

  const string ref_fasta{argv[1]};
  const Reference reference{ref_fasta};
  const Mappability mappability{ref_fasta, true};
  const ChromosomeIndexLookup lookup{reference};

  const string input_name{argv[2]};
  ifstream input{input_name.c_str()};
  if (!input) throw Error("Cannot open input") << input_name;

  unsigned int chr;
  unsigned int chr_offset{0};
  string last_chr;
  string chr_name;
  unsigned int start;
  unsigned int stop;
  unsigned int distance;
  unsigned int length;
  string motif;
  unsigned int n_motif{0};
  unsigned int n_unique_motif{0};
  unsigned int n_subseq{0};
  unsigned int n_good_score{0};
  unsigned int n_good_unique{0};
  const unsigned int notify{10000};

  time_t start_time{time(nullptr)};

  vector<string> saved;

  const string sequences_name{input_name + ".temp"};
  ofstream sequence_out(sequences_name.c_str());

  input.ignore(1000, '\n');
  while (input >> chr_name >> start >> stop >> distance >> length >> motif) {
    ++n_motif;

    if (chr_name != last_chr) {
      chr = lookup[chr_name];
      last_chr = chr_name;
      chr_offset = reference.offset(chr);
    }

    // Don't bother if motif does not have unique sequence
    const unsigned int abspos{chr_offset + start};
    if (mappability.low(abspos) > min(length, 253U)) {
      continue;
    }

    ++n_unique_motif;

    // Find sequences with good GCAT score
    const unsigned int target_score{60};
    const unsigned int target_slop{2};
    for (unsigned int start_base{0};
         start_base + 10 < motif.size(); ++start_base) {
      unsigned int score{0};
      string good_sequence;
      ++n_subseq;
      for (unsigned int b{start_base}; b != motif.size(); ++b) {
        score += gcat_score(motif[b]);
        if (score + target_slop < target_score) continue;
        if (score > target_score + target_slop) break;
        const unsigned int slength{b - start_base + 1};
        if (score + target_slop == target_score &&
            b + 1 != motif.size() &&
            slength != 30 &&
            gcat_score(motif[b + 1]) == 2) continue;
        if (slength < 20) continue;
        if (slength > 30) break;
        good_sequence = motif.substr(start_base, slength);
        break;
      }

      // Is sequence good
      if (good_sequence.size()) {
        ++n_good_score;

        // Is sequence sufficiently unique
        const unsigned int sequence_abspos{abspos + start_base};
        if (mappability.low(sequence_abspos + 2) >
            min(good_sequence.size() - 4, 253UL)) {
          continue;
        }
        ++n_good_unique;

        const unsigned int n_expected{20427871/8};
        if ((n_good_unique % notify) == 0) {
          const time_t current_time{time(nullptr)};
          const uint64_t elapsed(current_time - start_time);
          const double minutes_elapsed{elapsed / 60.0};
          const double fraction_done{1.0 * n_good_unique / n_expected};
          serr << "processed" << n_good_unique
               << "of expected" << n_expected
               << "or" << fraction_done << "of total"
               << "in" << minutes_elapsed << "minutes"
               << minutes_elapsed * (1 - fraction_done) / fraction_done
               << "minutes remaining"
               << endl;
        }

        ostringstream info;
        const unsigned int position{start + start_base};
        info << good_sequence << " "
             << good_sequence.size() << " "
             << score << " "
             << chr_name << " "
             << position << " "
             << position + good_sequence.size() << " "
             << distance;
        saved.push_back(info.str());
        sequence_out << good_sequence << endl;
      }
    }
  }

  serr << n_unique_motif << "of" << n_motif << "motifs were unique" << endl;
  serr << n_good_score << "of" << n_subseq
       << "subsequences tried had a good score" << endl;
  serr << n_good_unique << "of" << n_good_score
       << "good scores were unique" << endl;

  sequence_out.close();
  const std::string bowtie_command{
    string("/data/software/bowtie/bowtie-0.12.8/bowtie ") +
        "-a --best -v 3 hg19 --suppress 2,3,4,6,7 --quiet -r -p 16 " +
        sequences_name + " | perl -pe 's/,/ /g'"};
  cerr << "running " << bowtie_command << endl;
  ipstream bowtie{bowtie_command};

  start_time = time(nullptr);
  string alignment;
  std::vector<unsigned int> mismatches(4);
  const unsigned int bad_n{static_cast<unsigned int>(-1)};
  unsigned int last_n{bad_n};
  unsigned int n_processed{0};
  const unsigned int n_expected{static_cast<unsigned int>(saved.size())};
  while (getline(bowtie, alignment)) {
    istringstream line_to_parse{alignment};
    unsigned int n;
    string sequence;
    string mismatch;
    unsigned int n_mismatches{0};
    line_to_parse >> n >> sequence;
    if (n != last_n) {
      if (last_n != bad_n) {
        ++n_processed;
        if ((n_processed % notify) == 0) {
          const time_t current_time{time(nullptr)};
          const uint64_t elapsed(current_time - start_time);
          const double minutes_elapsed{elapsed / 60.0};
          const double fraction_done{1.0 * n_processed / n_expected};
          serr << "processed" << n_processed
               << "of expected" << n_expected
               << "or" << fraction_done << "of total"
               << "in" << minutes_elapsed << "minutes"
               << minutes_elapsed * (1 - fraction_done) / fraction_done
               << "minutes remaining"
               << endl;
        }

        const string & out = saved[n];
        cout << out << " "
             << mismatches[0] << " "
             << mismatches[1] << " "
             << mismatches[2] << " "
             << mismatches[3] << endl;
      }
      mismatches[0] = 0;
      mismatches[1] = 0;
      mismatches[2] = 0;
      mismatches[3] = 0;
      last_n = n;
    }
    while (line_to_parse >> mismatch) {
      ++n_mismatches;
    }
    ++mismatches[n_mismatches];
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



