//
// video
//
// play with X11 video
// just a thought - unimplemented
//
// Copyright 2018 Peter Andrews @ CSHL
//

#include <exception>
#include <iostream>
#include <fstream>
#include <vector>

#include "error.h"
#include "tsv.h"
#include "x11plot.h"

namespace paa {

class FrameInfo {
 public:
  FrameInfo() = default;

 private:
  double duration{1.0};
  function<void(void)> render{[]() {}};
};

class Video {
 public:
  Video(FrameInfo * begin_, FrameInfo * end_) :
      begin{begin_}, end{end_} { }

 private:
  FrameInfo * const begin;
  FrameInfo * const end;
}

}  // namespace paa

using std::cout;
using std::cerr;
using std::endl;
using std::exception;
using std::ofstream;
using std::vector;


using paa::Error;


int main(int argc, char** argv) try {
  if (--argc != 0) throw Error("usage: video");
  if (0) cout << argv[1];

  const vector<FrameInfo> frames;
  const Video video{frames.begin(), frames.end()};


  return 0;
} catch (Error & e) {
  cerr << "Error:" << endl;
  cerr << e.what() << endl;
  return 1;
} catch (exception & e) {
  cerr << "exception" << endl;
  cerr << e.what() << endl;
  return 1;
} catch (...) {
  cerr << "unknown exception was caught" << endl;
  return 1;
}
