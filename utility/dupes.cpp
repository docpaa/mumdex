//
// dupes.cpp
//
// Look for duplicate text
//
// Copyright 2021 Peter Andrews @ CSHL
//

#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <set>
#include <string>
#include <utility>
#include <vector>
#include <climits>

#include "error.h"
#include "files.h"
#include "paastrings.h"
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
using std::set;
using std::string;
using std::unique_ptr;
using std::vector;

using paa::Error;
using paa::file_string;
using paa::replace;
using paa::replace_all;

string reverse(string input) {
  std::reverse(input.begin(), input.end());
  return input;
}

int main(int argc, char ** argv) try {
  paa::exit_on_pipe_close();
  std::ios_base::sync_with_stdio(false);
  cin.tie(nullptr);

  const string usage{"usage: dupes [input_files ...]"};
  vector<string> files;
  while (--argc) {
    const string file_name{(++argv)[0]};
    cerr << "Load file " << file_name << endl;
    string contents{replace(file_string(file_name), '\n', ' ')};
    for (uint64_t i{0}; i != 5; ++i) {
      contents = replace_all(contents, "      ", " ");
      contents = replace_all(contents, "  ", " ");
    }
    files.push_back(contents);
  }
  map<string, uint64_t> counts;
  set<string> seen;
  for (const string & file : files) {
    cerr << "file size " << file.size() << endl;
    for (uint64_t b{0}; b != file.size(); ++b) {
      uint64_t e{b};
      uint64_t last_saw{0};
      uint64_t saw{0};
      for (; e != file.size(); ++e) {
        saw = 0;
        const string sub{file.substr(b, e - b + 1)};
        if (seen.insert(sub).second == false ||
            seen.insert(reverse(sub)).second == false) continue;
        for (const string & file2 : files) {
          uint64_t pos{0};
          while (pos != string::npos &&
                 (pos = file2.find(sub, pos)) != string::npos) {
            ++pos;
            ++saw;
          }
        }
        if (saw <= 1) break;
        last_saw = saw;
      }
      if (e > b && last_saw > 1) {
        counts[file.substr(b, e - b)] = last_saw;
        for (uint64_t a{b}; a != e; ++a) {
          const string sub{file.substr(a, e - a)};
          seen.insert(sub);
          seen.insert(reverse(sub));
        }
      }
    }
  }
  vector<pair<string, uint64_t>> repeats;
  for (const auto & item : counts)
    if (item.second > 1) repeats.push_back(item);
  sort(repeats.begin(), repeats.end(),
       [](const pair<string, uint64_t> & l,
          const pair<string, uint64_t> & r) {
         return r.first.size() < l.first.size();
       });
  set<string> ending;
  for (const auto & item : repeats) {
    if (ending.count(item.first) == 0)
      cout << "Repeat: " << item.second << " " << item.first.size() << " "
           << item.second * item.first.size() << " " << item.first << endl;
    for (uint64_t b{0}; b != item.first.size(); ++b)
      ending.insert(item.first.substr(b, item.first.size() - b + 1));
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
