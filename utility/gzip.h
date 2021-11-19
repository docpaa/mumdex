//
// gzip.h
//
// ostream interface to gzip
//
// copyright 2021 Peter Andrews
//

#ifndef PAA_GZIP_H_
#define PAA_GZIP_H_

#include <zlib.h>
#include <mutex>
#include <set>
#include <sstream>
#include <string>
#include <typeinfo>
#include <vector>
#include "error.h"

namespace paa {

// Buffer size for gzip operations
constexpr uint64_t gz_buffer{16 * 8192};

// Write output in gzip format.  Should just work.
struct Gzip {
  // Construction and destruction
  explicit Gzip(const std::string & file_name,
                const unsigned int compression = 5,
                const bool complain_ = true) :
      file{gzopen(file_name.c_str(),
                  ("wb" + std::to_string(compression)).c_str())},
      complain{complain_} {
        if (compression > 9)
          throw Error("Bad compression in Gzip") << compression;
        gzbuffer(file, gz_buffer);
      }
  ~Gzip() { close(); }
  Gzip(const Gzip &) = delete;
  Gzip & operator=(const Gzip &) = delete;
  Gzip(Gzip && other) : file{other.file}, complain{other.complain} {
    other.file = nullptr;
  }

  // User operations
  operator bool() const { return file; }
  void close() {
    if (file) {
      gzclose(file);
      file = nullptr;
    }
  }
  template <class Value>
  Gzip & operator<<(const Value & value) {
    if (complain) {
      static std::set<size_t> seen;
      const std::type_info & type{typeid(Value)};
      static std::mutex mutex;
      std::lock_guard<std::mutex> lock(mutex);
      if (seen.insert(type.hash_code()).second)
        std::cerr << "Writing gzip type " << type.name() << std::endl;
    }
    std::ostringstream stream;
    stream << value;
    gzfwrite(stream.str().c_str(), stream.str().size(), 1, file);
    return *this;
  }
  Gzip & operator<<(const bool value) {
    printf("%d", value);
    return *this;
  }
  Gzip & operator<<(const char value) {
    printf("%c", value);
    return *this;
  }
  Gzip & operator<<(const int value) {
    printf("%d", value);
    return *this;
  }
  Gzip & operator<<(const unsigned int value) {
    printf("%u", value);
    return *this;
  }
  Gzip & operator<<(const int64_t value) {
    printf("%ld", value);
    return *this;
  }
  Gzip & operator<<(const uint64_t value) {
    printf("%lu", value);
    return *this;
  }
  Gzip & operator<<(const double value) {
    printf("%f", value);
    return *this;
  }
  Gzip & operator<<(const char * value) {
    printf("%s", value);
    return *this;
  }
  Gzip & operator<<(const std::string & value) {
    printf("%s", value.c_str());
    return *this;
  }

 private:
  template <class Value>
  void printf(const char * const format, const Value value) {
    gzprintf(file, format, value);
  }
  gzFile file;
  bool complain;
};

// Read input in gzip format.  Must specify maximum expected line length.
// Input via getline or stream operations must not be interspersed badly.
struct Zcat {
  // Construction and destruction
  explicit Zcat(const std::string & file_name,
                const uint64_t max_line_length = 151) :
      file{gzopen(file_name.c_str(), "rb")},
      buffer(max_line_length + 3) {  // \n + 0 + overflow checker
        if (!file) throw Error("Zcat open failed") << file_name;
        gzbuffer(file, gz_buffer);
      }
  ~Zcat() { close(); }
  Zcat(const Zcat &) = delete;
  Zcat & operator=(const Zcat &) = delete;
  Zcat(Zcat && other) :
      file{other.file}, buffer{std::move(other.buffer)},
      stream{std::move(other.stream)} {
    other.file = nullptr;
  }

  // User operations
  operator bool() const { return file; }
  void close() {
    if (file) {
      gzclose(file);
      file = nullptr;
    }
  }
  Zcat & getline(std::string & line) {
    const char * gets_result{gzgets(file, &buffer[0], buffer.size())};
    if (!gets_result) {
      close();
    } else {
      line = gets_result;
      if (line.size() == buffer.size() - 1)
        throw Error("Max line size reached in Zcat");
      line.resize(line.size() - 1);  // Remove newline character
    }
    return *this;
  }
  template<class Value>
  Zcat & operator>>(Value & value) {
    if (!stream || !stream.rdbuf()->in_avail()) {
      std::string line;
      getline(line);
      if (!*this) return *this;
      stream.str(line);
      stream.clear();
    }
    stream >> value;
    return *this;
  }

 private:
  Zcat & getc(char & value) {
    const int result{gzgetc(file)};
    if (result < 0) {
      close();
    } else {
      value = result;
    }
    return *this;
  }
  gzFile file;
  std::vector<char> buffer;
  std::istringstream stream{};
};
Zcat & getline(Zcat & zcat, std::string & line) { return zcat.getline(line); }

}  // namespace paa

#endif  // PAA_GZIP_H_
