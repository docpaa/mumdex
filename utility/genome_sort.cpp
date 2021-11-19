//
// genome_sort.cpp
//
// sort input by genome coordinates
// First columns are chr, chrpos
//
// Copyright 2021 Peter Andrews @ CSHL
//

#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#include <climits>

#include "error.h"
#include "mumdex.h"
#include "utility.h"

using std::cerr;
using std::cin;
using std::cout;
using std::endl;
using std::exception;
using std::ifstream;
using std::istream;
using std::istringstream;
using std::make_unique;
using std::map;
using std::pair;
using std::string;
using std::unique_ptr;
using std::vector;

using paa::ChromosomeIndexLookup;
using paa::Error;
using paa::Reference;

int main(int argc, char ** argv) try {
  paa::exit_on_pipe_close();
  std::ios_base::sync_with_stdio(false);
  cin.tie(nullptr);

  const string usage{
    "usage: genome_sort [-h] [-d D] reference [input_files ...]"};

  // Process optional arguments
  --argc;
  char delimeter_{0};
  bool header_{false};
  while (argc) {
    bool acted{false};
    if (argc > 1 && argv[1] == string("-d")) {
      delimeter_ = argv[2][0];
      argc -= 2;
      argv += 2;
      acted = true;
    }
    if (argc >= 1 && argv[1] == string("-h")) {
      --argc;
      ++argv;
      acted = header_ = true;
    }
    if (!acted) break;
  }
  const bool header{header_};

  // Load reference
  if (argc < 1) throw Error(usage);
  const Reference reference{argv[1]};
  const ChromosomeIndexLookup chr_lookup{reference};
  --argc;
  ++argv;

  // Open inputs
  vector<unique_ptr<ifstream>> input_files;
  vector<istream *> inputs;
  if (argc) {
    while (argc--) {
      const string file_name{(argv++)[1]};
      input_files.push_back(make_unique<ifstream>(file_name.c_str()));
      if (!*input_files.back())
        throw Error("Could not open file for input") << file_name;
      inputs.push_back(input_files.back().get());
    }
  } else {
    inputs.push_back(&cin);
  }

  // Determine delimeter automatically, if one was not specified on command line
  unique_ptr<istringstream> first_line;
  const char delimeter{[delimeter_, &inputs, &first_line]() {
      if (delimeter_) {
        return delimeter_;
      } else {
        string line;
        getline(*inputs.front(), line);
        map<char, size_t> spaces;
        for (const char c : line) {
          if (isspace(c) || c == ',') {
            ++spaces[c];
          }
        }
        line += '\n';
        first_line = make_unique<istringstream>(line.c_str());
        inputs.insert(inputs.begin(), first_line.get());
        if (spaces.size() == 1) {
          return spaces.begin()->first;
        } else if (spaces.find('\t') != spaces.end()) {
          return '\t';
        } else if (spaces.find(',') != spaces.end()) {
          return ',';
        } else if (spaces.find(' ') != spaces.end()) {
          return ' ';
        } else {
          return '\t';
        }
      }
    }()};

  // Output header if exists
  string line;
  if (header) {
    getline(*inputs.front(), line);
    cout << line << endl;
  }

  vector<pair<unsigned int, string>> lines;
  string chr_name;
  unsigned int pos;
  for (istream * in : inputs) {
    while (getline(*in, line)) {
      istringstream line_stream{line};
      getline(line_stream, chr_name, delimeter);
      line_stream >> pos;
      const unsigned int chr{chr_lookup[chr_name]};
      const unsigned int abspos{reference.abspos(chr, pos)};
      lines.emplace_back(abspos, line);
    }
  }
  sort(lines.begin(), lines.end());

  // Output sorted lines
  for (const auto & info : lines) cout << info.second << endl;

  return 0;
} catch (Error & e) {
  cerr << "paa::Error:" << endl;
  cerr << e.what() << endl;
  return 1;
} catch(exception & e) {
  cerr << e.what() << endl;
  return 1;
} catch(...) {
  cerr << "Some exception was caught." << endl;
  return 1;
}
