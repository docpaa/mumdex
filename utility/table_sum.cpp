//
// table_sum.cpp
//
// Sum cells of a table
//
// Copyright 2021 Peter Andrews @ CSHL
//

#include <fstream>
#include <iostream>
#include <map>
#include <memory>
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
using std::ifstream;
using std::istringstream;
using std::make_unique;
using std::map;
using std::string;
using std::unique_ptr;
using std::vector;

using paa::Error;

int main(int argc, char ** argv) try {
  paa::exit_on_pipe_close();
  std::ios_base::sync_with_stdio(false);
  cin.tie(nullptr);

  // Process optional arguments
  --argc;
  char delimeter{0};
  vector<string> protect_column_names;
  bool not_negative{false};
  bool average{false};
  while (argc) {
    bool acted{false};
    if (argc > 2 && argv[1] == string("-p")) {
      protect_column_names.push_back(argv[2]);
      argc -= 2;
      argv += 2;
      acted = true;
    }
    if (argc > 2 && argv[1] == string("-d")) {
      delimeter = argv[2][0];
      argc -= 2;
      argv += 2;
      acted = true;
    }
    if (argc > 1 && argv[1] == string("-n")) {
      not_negative = true;
      argc -= 1;
      argv += 1;
      acted = true;
    }
    if (argc > 1 && argv[1] == string("-a")) {
      average = true;
      argc -= 1;
      argv += 1;
      acted = true;
    }
    if (!acted) break;
  }

  // Check for proper usage
  const string usage{"usage: table_sum [-d delim] [-p col] ... table ..."};
  if (argc < 1) throw Error(usage);

  // Open input files
  vector<unique_ptr<ifstream>> tables;
  string header{""};
  while (argc--) {
    const string file_name{argv++[1]};
    tables.push_back(make_unique<ifstream>(file_name));
    if (!tables.back()) throw Error("Problem opening input") << file_name;
    string header_;
    getline(*tables.back(), header_);
    if (header.size()) {
      if (header != header_) throw Error("Header mismatch in") << file_name;
    } else {
      header = header_;
      if (!delimeter) {
        map<char, size_t> spaces;
        for (const char c : header) {
          if (isspace(c) || c == ',') {
            ++spaces[c];
          }
        }
        if (spaces.size() == 1) {
          delimeter = spaces.begin()->first;
        } else if (spaces.find('\t') != spaces.end()) {
          delimeter = '\t';
        } else if (spaces.find(',') != spaces.end()) {
          delimeter = ',';
        } else if (spaces.find(' ') != spaces.end()) {
          delimeter = ' ';
        } else {
          delimeter = '\t';
        }
      }
    }
  }

  // Columns to protect - just keep value of last table cell
  const vector<uint64_t> protect_column{
    [&header, &protect_column_names, delimeter]() {
      vector<uint64_t> result;
      string column;
      istringstream header_stream{header.c_str()};
      while (getline(header_stream, column, delimeter))
        result.push_back(find(protect_column_names.begin(),
                              protect_column_names.end(), column) !=
                         protect_column_names.end());
      return result;
    }()};

  // Process table
  string line;
  cout << header << "\n";
  string sval;
  double dval{0};
  int64_t ival{0};
  bool done{false};
  while (true) {
    for (uint64_t c{0}; c != protect_column.size(); ++c) {
      double dsum{0};
      int64_t isum{0};
      bool is_int{true};
      const char this_delimeter{
        c + 1 == protect_column.size() ? '\n' : delimeter};
      const bool protect{static_cast<bool>(protect_column[c])};
      uint64_t n{0};
      for (uint64_t t{0}; t != tables.size(); ++t) {
        getline(*tables[t], sval, this_delimeter);
        if (!*tables[t]) {
          if (c == 0 && t == 0) {
            done = true;
            break;
          }
          throw Error("Problem getting cell in table") << t << "col" << c;
        }
        if (protect) continue;
        istringstream dval_stream{sval};
        dval_stream >> dval;
        if (!dval_stream) throw Error("Problem parsing") << sval << dval;
        if (!not_negative || dval >= 0) {
          ++n;
          dsum += dval;
        }
        if (is_int) {
          if (sval.find('.') != string::npos) is_int = false;
          istringstream ival_stream{sval};
          ival_stream >> ival;
          if (!not_negative || ival >= 0)
            isum += ival;
        }
      }
      if (done) break;
      if (c) cout << delimeter;
      if (protect) {
        cout << sval;
      } else {
        if (average) {
          if (is_int) {
            cout << (n ? 1.0 * isum / n : -1);
          } else {
            cout << (n ? dsum / n : -1);
          }
        } else {
          if (is_int) {
            cout << isum;
          } else {
            cout << dsum;
          }
        }
      }
    }
    if (done) break;
    cout << "\n";
    if (!*tables.front()) break;
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
