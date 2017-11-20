//
// bridges2txt
//
// convert bridge information to text format
//
// Copyright 2016 Peter Andrews CSHL
//

#include <exception>
#include <iostream>
#include <string>

#include "bridges.h"
#include "error.h"
#include "files.h"
#include "mumdex.h"
#include "utility.h"

using std::cerr;
using std::cout;
using std::endl;
using std::exception;

using paa::ref_ptr;
using paa::sout;
using paa::BridgeInfo;
using paa::Error;
using paa::MappedVector;
using paa::Reference;

const Reference * paa::ref_ptr;

int main(int argc, char* argv[])  try {
  paa::exit_on_pipe_close();  // Handle pipe close signal properly

  if (argc != 3) {
    throw Error("usage: bridges2txt ref bridge_file");
  }

  const Reference ref{argv[1]};
  ref_ptr = &ref;

  const MappedVector<BridgeInfo> bridges{argv[2]};

  sout << "chr" << "pos" << "high"
       << "chr2" << "pos2" << "high2"
       << "it" << "inv" << "ioff"
       << "al" << "bl"
       << "aml" << "bml"
       << "bc" << "amc" << "bmc" << endl;
  for (const BridgeInfo & bridge : bridges) {
    bridge.output(sout);
    sout << endl;
    // sout << bridge << endl;  // Problem on mac...
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

