//
// uniqc.cpp
//
// Count occurrences of lines
//
// Copyright 2021 Peter Andrews @ CSHL
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

  const string usage{"usage: uniqc [input_files ...]"};
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

  using LineCounts = map<string, uint64_t>;
  const LineCounts line_counts{[&inputs]() {
      auto count_lines = [](istream & in) {
        LineCounts result;
        string line;
        while (getline(in, line)) ++result[line];
        return result;
      };
      vector<future<LineCounts>> futures;
      ThreadPool pool{32};
      for (istream * in : inputs)
        futures.push_back(pool.run(count_lines, ref(*in)));
      LineCounts result;
      for (future<LineCounts> & future : futures)
        for (const auto & item : future.get())
          result[item.first] += item.second;
      return result;
    }()};

  for (const auto & line_count : line_counts)
    cout << line_count.second << '\t' << line_count.first << '\n';

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
