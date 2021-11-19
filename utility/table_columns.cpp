//
// table_columns.cpp
//
// extract columns from a table
// auto-determines the delimeter where feasible
//
// Copyright 2020 Peter Andrews @ CSHL
//

#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "error.h"
#include "utility.h"

using std::cerr;
using std::cin;
using std::cout;
using std::endl;
using std::exception;
using std::istringstream;
using std::map;
using std::string;
using std::vector;

using paa::Error;

int main(int argc, char ** argv) try {
  paa::exit_on_pipe_close();
  std::ios_base::sync_with_stdio(false);
  cin.tie(nullptr);

  // Process optional arguments
  --argc;
  char delimeter_{0};
  while (argc) {
    bool acted{false};
    if (argc > 2 && argv[1] == string("-d")) {
      delimeter_ = argv[2][0];
      argc -= 2;
      argv += 2;
      acted = true;
    }
    if (!acted) break;
  }

  // Check for proper usage
  const string usage{"usage: table_columns [-d D] column_name ..."};
  if (argc < 1) throw Error(usage);

  // Table header line
  const string header{[]() {
      string result;
      getline(cin, result);
      return result;
    }()};

  // Determine delimeter automatically, if one was not specified on command line
  const char delimeter{[delimeter_, &header]() {
      if (delimeter_) {
        return delimeter_;
      } else {
        map<char, size_t> spaces;
        for (const char c : header) {
          if (isspace(c) || c == ',') {
            ++spaces[c];
          }
        }
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

  // Column name - column index lookup table
  const map<string, uint64_t> column_indexes{[&header, delimeter]() {
      map<string, uint64_t> result;
      string column;
      istringstream header_stream{header.c_str()};
      while (getline(header_stream, column, delimeter))
        result.emplace(column, result.size());
      return result;
    }()};

  // Column names to output
  const vector<string> column_names{[argc, argv]() {
      vector<string> result;
      for (int c{0}; c != argc; ++c) result.push_back(argv[c + 1]);
      return result;
    }()};

  // Output column indexes
  const vector<uint64_t> output_columns{
    [&header, &column_names, &column_indexes]() {
      vector<uint64_t> result;
      for (const string & column : column_names)
        try {
          result.push_back(column_indexes.at(column));
        } catch (...) {
          throw Error("Column name") << column << "not found in" << header;
        }
      return result;
    }()};

  // Process table
  string line;
  vector<string> row(column_indexes.size());
  cout << column_names[0];
  for (uint64_t c{1}; c != column_names.size(); ++c)
    cout << delimeter << column_names[c];
  cout << "\n";
  while (getline(cin, line)) {
    istringstream line_stream{line};
    for (string & column : row) getline(line_stream, column, delimeter);
    cout << row[output_columns[0]];
    for (uint64_t c{1}; c != output_columns.size(); ++c)
      cout << delimeter << row[output_columns[c]];
    cout << "\n";
  }

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
