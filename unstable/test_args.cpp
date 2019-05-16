//
// test_args.cpp
//
// Fancy constructors with "optional" arguments in any order
//
// Copyright 2019 Peter Andrews
//

#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <climits>

#include "error.h"
#include "utility.h"

using std::cerr;
using std::cin;
using std::cout;
using std::endl;
using std::exception;
using std::ifstream;
using std::istream;
using std::make_unique;
using std::string;
using std::unique_ptr;
using std::vector;

using paa::Error;

struct Color {
  string value;
};
struct Size {
  double value;
};

class FlexArgs {
 public:
  template <class ... Args>
  explicit FlexArgs(Args && ... args) {
    process(std::forward<Args>(args)...);
    cout << color.value << " " << size.value << " " << truth << endl;
  }
  template <class Arg1, class ... Args>
  void process(Arg1 && arg1, Args && ... args) {
    process(arg1);
    process(std::forward<Args>(args)...);
  }
  void process(const Color & color_) {
    color = color_;
    return;
  }
  void process(const Size & size_) {
    size = size_;
    return;
  }
  void process(const bool truth_) {
    truth = truth_;
    return;
  }
  void process() const { }

 private:
  Color color{"blue"};
  Size size{3.0};
  bool truth{false};
};

int main(int argc, char ** argv) try {
  paa::exit_on_pipe_close();
  std::ios_base::sync_with_stdio(false);
  cin.tie(nullptr);

  if (--argc != 0) throw Error("usage: test_args") << argv[1];

  const Color red{"red"};
  const Size big{20.0};
  const FlexArgs test1{};
  const FlexArgs test2(big, red);
  const FlexArgs test3(true, big);
  const FlexArgs test4{red, big, true};

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
