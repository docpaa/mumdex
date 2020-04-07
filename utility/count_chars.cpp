//
// count_characters.cpp
//
// Count occurrences of characters, separated by whitespace
//
// Copyright 2020 Peter Andrews @ CSHL
//

#include <fstream>
#include <future>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "error.h"
#include "threads.h"
#include "utility.h"

using std::cerr;
using std::cin;
using std::cout;
using std::endl;
using std::exception;
using std::future;
using std::ifstream;
using std::istream;
using std::make_unique;
using std::map;
using std::pair;
using std::string;
using std::unique_ptr;
using std::vector;

using paa::Error;
using paa::ThreadPool;

int main(int argc, char ** argv) try {
  paa::exit_on_pipe_close();
  std::ios_base::sync_with_stdio(false);
  cin.tie(nullptr);

  const string usage{"usage: count_characters [input_files ...]"};
  --argc;

  // Open inputs
  vector<unique_ptr<ifstream>> input_files;
  vector<istream *> inputs;
  if (argc) {
    while (argc--) {
      const string file_name{(argv++)[1]};
      input_files.push_back(make_unique<ifstream>(file_name.c_str()));
      if (!*input_files.back())
        throw Error(usage + "\nCould not open file for input") << file_name;
      inputs.push_back(input_files.back().get());
    }
  } else {
    inputs.push_back(&cin);
  }

  using CharCounts = map<char, uint64_t>;
  const CharCounts char_counts{[&inputs]() {
      auto count_chars = [](istream & in) {
        CharCounts result;
        char character;
        while (in >> character) ++result[character];
        return result;
      };
      vector<future<CharCounts>> futures;
      ThreadPool pool{32};
      for (istream * in : inputs)
        futures.push_back(pool.run(count_chars, ref(*in)));
      CharCounts result;
      for (future<CharCounts> & future : futures)
        for (const auto & item : future.get())
          result[item.first] += item.second;
      return result;
    }()};

  for (const pair<char, uint64_t> & char_count : char_counts)
    cout << char_count.first << '\t' << char_count.second << '\n';

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
