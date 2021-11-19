//
// error.h
//
// exception class for reporting errors
//
// Copyright 2012 Peter Andrews @ CSHL
//

#ifndef PAA_UTILITY_ERROR_H_
#define PAA_UTILITY_ERROR_H_

#include <exception>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace paa {

// Output a vector
template <class Element>
std::ostream & operator<<(std::ostream & out,
                          const std::vector<Element> & elements) {
  for (uint64_t e{0}; e != elements.size(); ++e) {
    if (e) out << " ";
    out << elements[e];
  }
  return out;
}

// Class to handle error reporting
class Error: public std::exception {
 public:
  // Initialize error with a message
  explicit Error(const std::string & m, const bool immediate = false) :
      std::exception(), message(m) {
    if (immediate) {
      std::cerr << "paa::Error: " << message << std::endl;
      exit(1);
    }
  }

  // Destructor
  virtual ~Error() throw() { }

  // Make a copy
  Error(const Error & t) : std::exception(), message(t.what()) {}

  // Returns the message
  virtual const char * what() const throw() {
    return message.c_str();
  }

  // Allows you to write to the error like it was a stream
  template <class IN> Error & operator << (const IN & in) {
    std::ostringstream out;
    out << " " << in;
    message += out.str();
    return *this;
  }

 private:
  std::string message;
};

class UsageError : public Error {
 public:
  using Error::Error;
  virtual ~UsageError() throw() { }
};

#if 0
class warn {
 public:
  explicit warn(const std::string & message) {
    std::cerr << message << std::endl;
  }
};
#endif

}  // namespace paa

#endif  // PAA_UTILITY_ERROR_H_
