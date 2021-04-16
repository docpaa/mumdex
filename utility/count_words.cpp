//
// count_words.cpp
//
// Count occurrences of words, separated by whitespace
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

  const string usage{"usage: count_words [input_files ...]"};
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

  using WordCounts = map<string, uint64_t>;
  const WordCounts word_counts{[&inputs]() {
      auto count_words = [](istream & in) {
        WordCounts result;
        string word;
        while (in >> word) ++result[word];
        return result;
      };
      vector<future<WordCounts>> futures;
      ThreadPool pool{32};
      for (istream * in : inputs)
        futures.push_back(pool.run(count_words, ref(*in)));
      WordCounts result;
      for (future<WordCounts> & future : futures)
        for (const auto & item : future.get())
          result[item.first] += item.second;
      return result;
    }()};

  for (const auto & word_count : word_counts)
    cout << word_count.first << '\t' << word_count.second << '\n';

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
