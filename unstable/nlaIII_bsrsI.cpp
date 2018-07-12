//
// nlaIII_bsrsI
//
// Get info on restriction fragments
//
// Copyright 2018 Peter Andrews @ CSHL
//

#include <exception>
#include <iostream>
#include <string>
#include <vector>

#include "error.h"
#include "mumdex.h"
#include "utility.h"

using std::cerr;
using std::cin;
using std::cout;
using std::endl;
using std::exception;
using std::string;
using std::vector;

using paa::sout;
using paa::ChromosomeIndexLookup;
using paa::Error;
using paa::Mappability;
using paa::Reference;

int main(int argc, char* argv[])  try {
  if (--argc != 1) throw Error("usage: nlaIII_bsrsI ref_fasta");

  const string ref_fasta{argv[1]};
  const Reference ref{ref_fasta};
  const Mappability mappability{ref};
  const ChromosomeIndexLookup chr_lookup{ref};
  const string nlaIII{"CATG"};
  const unsigned int nlaIII_cut{4};
  const string bsrsI{"CCGCTC"};
  const string bsrsI_rc{paa::reverse_complement(bsrsI)};

  cout << "chr\tstart\tstop\tid\tlength\t"
       << "broken\tlow_map\thigh_map\thigh_map2\n";

  unsigned int fragment_id{0};
  for (const string & cc :
    {"1", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "2",
          "20", "21", "22", "3", "4", "5", "6", "7", "8", "9", "M", "X", "Y"}) {
    const unsigned int c{chr_lookup[string("chr") + cc]};
    unsigned int last_position{0};
    for (unsigned int b{0}; b != ref.size(c); ++b) {
      bool match{true};
      for (unsigned int b2{0}; b2 != nlaIII.size(); ++b2) {
        if (ref[c][b + b2] != nlaIII[b2]) {
          match = false;
          break;
        }
      }
      if (match) {
        const unsigned int position{b + nlaIII_cut};
        if (last_position != 0) {
          bool broken{false};
          for (const string & seq : {bsrsI, bsrsI_rc}) {
            for (unsigned int b3{last_position}; b3 != position + 1; ++b3) {
              bool match2{true};
              for (unsigned int b4{0}; b4 != seq.size(); ++b4) {
                if (ref[c][b3 + b4] != seq[b4]) {
                  match2 = false;
                  break;
                }
              }
              if (match2) {
                broken = true;
                break;
              }
            }
            if (broken) break;
         }
          const unsigned int fragment_end{position - nlaIII_cut};
          cout << ref.name(c)
               << '\t' << last_position + 1
               << '\t' << position - nlaIII_cut
               << '\t' << fragment_id++
               << '\t' << fragment_end - last_position
               << '\t' << broken
               << '\t' << mappability.low(ref.abspos(c, last_position))
               << '\t' << mappability.high(ref.abspos(c, fragment_end -1))
               << '\t' << mappability.high(ref.abspos(c, position -1))
               << '\n';
        }
        last_position = position;
      }
    }
    ++fragment_id;
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



