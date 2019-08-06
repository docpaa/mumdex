//
// test_args.cpp
//
// Fancy constructors with "optional" arguments in any order
//
// Copyright 2019 Peter Andrews @ CSHL
//

#include <iostream>
#include <string>

#include "error.h"
#include "utility.h"

using std::cerr;
using std::cin;
using std::cout;
using std::endl;
using std::exception;
using std::string;

using paa::Error;

struct Color {
  string value;
};
struct Size {
  double value;
};

// Probably a bad idea - inefficient to create members and then set them later
class FlexArgs {
 public:
  template <class ... Args>
  explicit FlexArgs(Args && ... args) {
    cout << endl <<  __PRETTY_FUNCTION__ << endl;
    process_args(std::forward<Args>(args)...);
    cout << color.value << " " << size.value
         << " " << (truth ? "true" : "false") << endl;
  }
#if 0
  template <class Arg1>
  void process_args(Arg1 && arg1) {
    cout << __PRETTY_FUNCTION__ << endl;
    process_arg(std::forward<Arg1>(arg1));
  }
#endif
  template <class Arg1, class ... Args>
  void process_args(Arg1 && arg1, Args && ... args) {
    cout << __PRETTY_FUNCTION__ << endl;
    process_arg(std::forward<Arg1>(arg1));
    process_args(std::forward<Args>(args)...);
  }
  void process_arg(const Color & color_) {
    cout << __PRETTY_FUNCTION__ << endl;
    color = color_;
    return;
  }
  // template <class SIZE>
  void process_arg(const Size & size_) {
    cout << __PRETTY_FUNCTION__ << endl;
    size = size_;
    return;
  }
  void process_arg(Size && size_) {
    cout << __PRETTY_FUNCTION__ << endl;
    size = std::move(size_);
    return;
  }
  void process_arg(const bool truth_) {
    cout << __PRETTY_FUNCTION__ << endl;
    truth = truth_;
    return;
  }
  void process_args() const {
    cout << __PRETTY_FUNCTION__ << endl;
  }

 private:
  Color color{"blue"};
  Size size{1.0};
  bool truth{false};
};

int main(int argc, char ** argv) try {
  if (--argc != 0) throw Error("usage: test_args") << argv[1];

  const Color red{"red"};
  Size big{2.0};
  const Size bigger{3.0};
  const FlexArgs test1{};
  const FlexArgs test2(big, red);
  const FlexArgs test3(true, big);
  const FlexArgs test4{red, big, true, bigger};
  const FlexArgs test5{red, big, true, Size{0.5}};

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
