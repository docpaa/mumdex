//
// sssa
//
// species specific suffix array
//
// Copyright 2016 Peter Andrews @ CSHL
//

#include "sssa.h"

#include <exception>
#include <iostream>
#include <string>

#include "error.h"

using std::cout;
using std::cerr;
using std::endl;
using std::exception;
using std::string;

using paa::Error;
using SSSA = paa::SSSA<24>;

int main(int argc, char* argv[]) try {
  using Species = paa::Species_t<15>;
  if (0) cout << "max_n " << Species::max_n() << " "
              << "n_values " << Species::n_values() << " "
              << "n_bits " << Species::n_bits() << " "
              << "mask " << static_cast<unsigned int>(Species::mask()) << " "
              << endl;
  if (0) cout << "Kmer size " << sizeof(SSSA::Kmer) << " and "
              << "Species size " << sizeof(SSSA::Species) << endl;

  cout << paa::SSSA<20>::info() << endl;
  cout << paa::SSSA<22>::info() << endl;
  cout << paa::SSSA<24>::info() << endl;
  cout << paa::SSSA<26>::info() << endl;
  cout << paa::SSSA<28>::info() << endl;
  cout << paa::SSSA<30>::info() << endl;
  if (!--argc) throw Error("usage: sssa genomes ...");
  SSSA sssa;
  while (argc--) {
    const string genome{argv++[1]};
    const SSSA sssa2{genome};
    sssa.merge(sssa2);
    cout << "Created sssa of size " << sssa2.size()
         << " for " << genome << " and merged size is " << sssa.size()
         << " with kmer density of " << sssa2.kmer_density() << endl;
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
