//
// reference_sequences.cpp
//
// Get sequences from the reference
//
// Copyright 2020 Peter Andrews @ CSHL
//

#include <exception>
#include <iostream>
#include <sstream>
#include <string>

#include "error.h"
#include "mumdex.h"

using std::cout;
using std::cerr;
using std::endl;
using std::exception;
using std::istringstream;
using std::string;

using paa::reverse_complement;
using paa::ChromosomeIndexLookup;
using paa::Error;
using paa::Mappability;
using paa::Reference;

int main(int argc, char * argv[]) try {
  // Check initial command line arguments
  const string usage{"usage: reference_sequences [-rc] reference "
        "chr:start-stop|chr:start:n ..."};
  if (--argc < 2) throw Error(usage);
  bool rc{false};
  if (string(argv[1]) == "-rc") {
    rc = true;
    --argc;
    ++argv;
  }
  const Reference reference{argv[1]};
  const ChromosomeIndexLookup chr_lookup{reference};
  const Mappability mappability{reference};
  --argc;
  argv += 2;

  // Process all sequence specifications
  string chr_name;
  unsigned int start;
  char type;
  unsigned int stop_or_n;
  unsigned int stop;
  unsigned int n;
  cout << "region\tsequence\tsize\tunique_length\n";
  for (int a{0}; a != argc; ++a) {
    istringstream specification{argv[a]};
    getline(specification, chr_name, ':');
    const unsigned int chr{chr_lookup[chr_name]};
    specification >> start >> type >> stop_or_n;
    if (!specification) throw Error("Specification parse error for") << argv[a];
    --start;
    if (type == ':') {
      stop = start + stop_or_n;
      n = stop_or_n;
    } else if (type == '-') {
      stop = stop_or_n - 1;
      n = stop - start;
    } else {
      throw Error(usage + "\nBad specification") << argv[a];
    }
    if (stop <= start) throw Error("Start does not precede stop") << argv[a];
    if (stop > reference.size(chr) || start >= reference.size(chr))
      throw Error("Region goes off end of chromosome") << argv[a];
    const unsigned int abspos{reference.abspos(chr, start)};
    const string seq{reference.subseq(chr, start, stop)};
    cout << argv[a] << '\t' << (rc ? reverse_complement(seq) : seq)
         << '\t' << n << '\t' << mappability.min(abspos, n) << '\n';
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
