//
// mumdex_sequences
//
// extract read sequences from a mumdex
//
// Copyright 2015 Peter Andrews @ CSHL
//

#include <array>
#include <exception>
#include <iostream>
#include <string>

#include "error.h"
#include "mumdex.h"
#include "utility.h"

using std::array;
using std::cout;
using std::cerr;
using std::endl;
using std::exception;
using std::string;

using paa::Error;
using paa::MUMdex;

int main(int argc, char * argv[]) try {
  paa::exit_on_pipe_close();

  if (--argc != 1) throw Error("usage: mumdex_sequences mumdex_name");
  const string mumdex_name{argv[1]};

  cerr << "Output mumdex sequences for " << mumdex_name << endl;
  const MUMdex mumdex(mumdex_name);
  for (unsigned int i = 0; i != mumdex.n_pairs(); ++i) {
    const array<string, 2> sequences(mumdex.sequences(i));
    for (const auto & sequence : sequences) {
      cout << sequence << '\n';
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
