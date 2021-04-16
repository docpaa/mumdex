//
// lookup_bridge
//
// lookup a specific bridge in bridges files
//
// Copyright 2016 Peter Andrews @ CSHL
//

#include <algorithm>
#include <exception>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>

#include "bridges.h"
#include "error.h"
#include "files.h"
#include "mumdex.h"
#include "population.h"
#include "utility.h"

using std::cerr;
using std::cout;
using std::endl;
using std::exception;
using std::lower_bound;
using std::map;
using std::ostringstream;
using std::right;
using std::setw;
using std::string;

using paa::readable;
using paa::sout;
using paa::BridgeInfo;
using paa::ChromosomeIndexLookup;
using paa::Error;
using paa::Family;
using paa::FileVector;
using paa::MUM;
using paa::Population;
using paa::Reference;
using paa::Sample;

int main(int argc, char* argv[])  try {
  // Check usage
  if (--argc < 11) {
    throw Error("usage: lookup_bridge ref bridges_dir pop_file "
                "chr1 pos1 high1 chr2 pos2 high2 inv [sample|family] ...");
  }

  // Process command line arguments
  const Reference ref{argv[1]};
  const ChromosomeIndexLookup chr_lookup{ref};
  const string bridges_dir{argv[2]};
  const Population pop{argv[3]};
  const string chr1s{argv[4]};
  const unsigned int pos1{static_cast<unsigned int>(atoi(argv[5]))};
  const bool high1{static_cast<bool>(atoi(argv[6]))};
  const string chr2s{argv[7]};
  const unsigned int pos2{static_cast<unsigned int>(atoi(argv[8]))};
  const bool high2{static_cast<bool>(atoi(argv[9]))};
  const int64_t invariant{atol(argv[10])};

  argc -= 10;
  argv += 10;

  const unsigned int chr1{chr_lookup[chr1s]};
  const unsigned int chr2{chr_lookup[chr2s]};

  while (argc--) {
    // Get samples for sample or family name
    const string sample_or_family{argv++[1]};
    for (const Sample & sample : pop.samples(sample_or_family)) {
      // Get sample info
      const Family family{pop.family(sample)};
      const string sample_name{pop.sample(sample)};
      const string family_name{pop.family(family)};

      cerr << "checking " << sample_name << endl;

      // Open bridges file
      ostringstream bridges_name;
      bridges_name << bridges_dir << "/" << sample_name << "/"
                   << get_bridge_file_name(ref, chr1);
      const FileVector<BridgeInfo> bridges{bridges_name.str()};
      if (!readable(bridges_name.str()))
        throw Error("Could not open bridges file")
            << bridges_name.str() << paa::bridges_bad_message();

      // Advance to start position
      FileVector<BridgeInfo>::const_iterator start{
        lower_bound(bridges.begin(), bridges.end(), pos1,
                    [](const BridgeInfo & lhs, const unsigned int lpos) {
                      return lhs.pos1() < lpos;
                    })};

      // Check for matching bridge
      for (FileVector<BridgeInfo>::const_iterator test_iter{start};
           test_iter != bridges.end(); ++test_iter) {
        const BridgeInfo test{*test_iter};
        if (test.pos1() > pos1) break;
        if (test.pos1() == pos1 &&
            test.high1() == high1 &&
            test.chr2() == chr2 &&
            test.pos2() == pos2 &&
            test.high2() == high2 &&
            test.invariant() == invariant) {
          sout << family_name << sample_name << pop.member(sample)
               << ref.name(test.chr1()) << test.pos1() << test.high1()
               << ref.name(test.chr2()) << test.pos2() << test.high2()
               << test.invariant() << test.offset() << test.bridge_count()
               << test.anchor1_length() << test.anchor2_length()
               << test.mate_anchor1_length() << test.mate_anchor2_length()
               << test.mate_anchor1_count() << test.mate_anchor1_count()
               << endl;
        }
      }
    }
  }

  cerr << "done" << endl;

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


