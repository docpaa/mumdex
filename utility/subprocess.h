//
// subprocess.h
//
// Read and write from subprocesses
//
// opopen works much better than opstream, which hangs sometimes
//
// copyright 2021 Peter Andrews
//

#ifndef PAA_EMPTY_H_
#define PAA_EMPTY_H_

#include <sys/wait.h>
#include <atomic>
#include <functional>
#include <future>
#include <iostream>
#include <set>
#include <string>
#include <typeinfo>
#include <vector>

#include "error.h"
#include "files.h"
#include "pstream.h"
#include "threads.h"
#include "utility.h"

namespace paa {

static std::mutex pmutex;

// Opens a named pipe (fifo) and deletes it on destruction
struct Fifo {
  explicit Fifo(const std::string & name__ = get_name()) : name_{name__} {
    unlink(name());
    mkfifo(name());
  }
  Fifo(Fifo && other) : name_{std::move(other.name_)} { other.name_ = ""; }
  Fifo(const Fifo &) = delete;
  Fifo & operator=(const Fifo &) = delete;
  ~Fifo() { if (name().size()) unlink(name()); }
  operator const std::string &() const { return name(); }
  const std::string & name() const { return name_; }
 private:
  std::string name_{get_name()};
  static std::atomic<uint64_t> n;
  static std::string base;
  static std::string get_name() { return base + std::to_string(++n); }
};
std::atomic<uint64_t> Fifo::n{0};
std::string Fifo::base("fifo." + std::to_string(getpid()) + ".");
struct Fifos {
  Fifo in_fifo{};
  Fifo out_fifo{};
};

// Runs a system command in the background
struct System {
  explicit System(const std::string & command_) :
      command{command_}, future{pasync([](const std::string & command__) {
        if (command__.empty()) throw Error("No command passed to system");
        return system(command__.c_str());
      }, command)} {}
  System(System && other) :
      command{std::move(other.command)}, future{std::move(other.future)} {
        other.command = "";
      }
  ~System() {
    if (command.size()) {
      const int status{future.get()};
      if (status)
        std::cerr << "System command error " << status << " for " << command
                  << std::endl;
    }
  }

 private:
  std::string command;
  std::future<int> future;
};

// Input stream from a process
struct isystem : private Fifo, private System, public std::ifstream {
  explicit isystem(const std::string & command_) :
      Fifo{}, System{"(" + command_ + ") >" + name()}, std::ifstream{name()} {}
  isystem(isystem && other) :
      Fifo{std::move(other)}, System{std::move(other)},
      std::ifstream{std::move(static_cast<std::ifstream&&>(other))} {}
};

// Output stream to a process
struct osystem : private Fifo, private System, public std::ofstream {
  explicit osystem(const std::string & command_) :
      Fifo{}, System{"(" + command_ + ") <" + name()}, std::ofstream{name()} {}
  osystem(osystem && other) :
      Fifo{std::move(other)}, System{std::move(other)},
      std::ofstream{std::move(other)} {}
  ~osystem() { close(); }
  void close() { std::ofstream::close(); }
};

// Bidirectional writing to and reading from a process
struct iosystem : private Fifos, private System, public std::ifstream {
  explicit iosystem(const std::string & command_) :
      System{"(" + command_ + ") >" + in_fifo.name() + " <" + out_fifo.name()},
      std::ifstream{in_fifo.name()}, out_{out_fifo} {}
  iosystem(iosystem && other) :
      Fifos{std::move(other)}, System{std::move(other)},
      std::ifstream{std::move(other)}, out_{std::move(other.out_)} {}
  std::ifstream & in() { return *this; }
  std::ofstream & out() { return out_; }
  void close() { out_.close(); }

 private:
  std::ofstream out_;
};
template <class Value>
std::ostream & operator<<(iosystem & stream, const Value & value) {
  return stream.out() << value;
}

struct ofile {
  ofile(const std::string & name_, const bool complain_) :
      name{name_}, complain{complain_} {}
  ~ofile() { close(); }
  ofile(const ofile &) = delete;
  ofile & operator=(const ofile &) = delete;
  ofile(ofile && other) :
      file{other.file}, name{std::move(other.name)}, complain{other.complain},
      closer{std::move(other.closer)} {
    other.file = nullptr;
  }
  operator bool() const { return file; }
  void open(int fd) {
    file = fdopen(fd, "wb");
    if (file == nullptr) {
      std::cerr << "fdopen problem in " << name << std::endl;
      ::close(fd);
      throw Error("Problem getting FILE * in") << name;
    } else {
      closer = [](const std::string & name_, FILE * file_) {
        if (::fclose(file_) == -1)
          std::cerr << "fclose failure in " << name_ << std::endl;
      };
    }
  }
  void open(const std::string & command) {
    file = popen(command.c_str(), "w");
    if (!file) throw Error("popen failed for") << command << "in" << name;
    closer = [](const std::string & name_, FILE * file_) {
      if (pclose(file_) == -1)
        std::cerr << "pclose failure in " << name_ << std::endl;
    };
  }
  void close() {
    if (file) {
      closer(name, file);
      file = nullptr;
    }
  }
  template <class Value>
  ofile & operator<<(const Value & value) {
    if (complain) {
      static std::set<size_t> seen;
      const std::type_info & type{typeid(Value)};
      static std::mutex mutex;
      std::lock_guard<std::mutex> lock(mutex);
      if (seen.insert(type.hash_code()).second)
        std::cerr << "Writing type " << type.name() << " in " << name
                  << std::endl;
    }
    std::ostringstream stream;
    stream << value;
    if (fwrite(stream.str().c_str(), stream.str().size(), 1, file) != 1)
      throw Error("Write error in") << name;
    return *this;
  }
  ofile & operator<<(const bool value) {
    printf("%d", value);
    return *this;
  }
  ofile & operator<<(const char value) {
    printf("%c", value);
    return *this;
  }
  ofile & operator<<(const int value) {
    printf("%d", value);
    return *this;
  }
  ofile & operator<<(const unsigned int value) {
    printf("%u", value);
    return *this;
  }
  ofile & operator<<(const int64_t value) {
    printf("%ld", value);
    return *this;
  }
  ofile & operator<<(const uint64_t value) {
    printf("%lu", value);
    return *this;
  }
  ofile & operator<<(const double value) {
    printf("%f", value);
    return *this;
  }
  ofile & operator<<(const char * value) {
    printf("%s", value);
    return *this;
  }
  ofile & operator<<(const std::string & value) {
    printf("%s", value.c_str());
    return *this;
  }

 private:
  template <class Value>
  void printf(const char * const format, const Value value) {
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wformat-nonliteral"
    if (fprintf(file, format, value) < 0)
      throw Error("fprintf error in") << name;
#pragma GCC diagnostic pop
  }
  FILE * file{nullptr};
  std::string name;
  bool complain{true};
  std::function<void(const std::string &, FILE *)> closer{};
};

struct ifile {
  ifile(const std::string & name_, const bool complain_) :
      name{name_}, complain{complain_} {}
  ~ifile() { close(); }
  ifile(const ifile &) = delete;
  ifile & operator=(const ifile &) = delete;
  ifile(ifile && other) :
      file{other.file}, name{std::move(other.name)}, complain{other.complain},
      closer{std::move(other.closer)},
      getline_line{other.getline_line}, getline_n{other.getline_n} {
        other.file = nullptr;
      }
  operator bool() const { return file; }
  void open(int fd) {
    file = fdopen(fd, "rb");
    if (file == nullptr) {
      std::cerr << "fdopen problem in " << name << std::endl;
      ::close(fd);
      throw Error("Problem getting FILE * in") << name;
    } else {
      closer = [](const std::string & name_, FILE * file_) {
        if (::fclose(file_) == -1)
          std::cerr << "Close failure in " << name_ << std::endl;
      };
    }
  }
  void open(const std::string & command) {
    file = popen(command.c_str(), "r");
    if (!file) throw Error("popen failed for") << command << "in" << name;
    closer = [](const std::string & name_, FILE * file_) {
      if (pclose(file_) == -1)
        std::cerr << "pclose failure in " << name_ << std::endl;
    };
  }
  void close() {
    if (file) {
      closer(name, file);
      file = nullptr;
      if (getline_line) free(getline_line);
    }
  }
  ifile & operator>>(bool & value) {
    int val;
    scanf("%d", val);
    value = val;
    return *this;
  }
  ifile & operator>>(char & value) {
    scanf("%c", value);
    return *this;
  }
  ifile & operator>>(int & value) {
    scanf("%d", value);
    return *this;
  }
  ifile & operator>>(unsigned int & value) {
    scanf("%u", value);
    return *this;
  }
  ifile & operator>>(int64_t & value) {
    scanf("%ld", value);
    return *this;
  }
  ifile & operator>>(uint64_t & value) {
    scanf("%lu", value);
    return *this;
  }
  ifile & operator>>(double & value) {
    scanf("%f", value);
    return *this;
  }
  ifile & operator>>(std::string & value) {
    if (!file) return *this;
    static thread_local char buffer[256];
    value.clear();
    while (true) {
      const int n_read{fscanf(file, "%255s", buffer)};
      if (n_read <= 0) {
        close();
        break;
      }
      value += buffer;
      if (n_read < 255) break;
      const int c{fgetc(file)};
      ungetc(c, file);
      if (c == EOF) close();
      if (isspace(c)) break;
    }
    return *this;
  }
  int peek() {
    const int c{fgetc(file)};
    ungetc(c, file);
    return c;
  }
  void ignore(const uint64_t n, const char delimeter) {
    for (uint64_t i{0}; i != n; ++i) {
      const int c{fgetc(file)};
      if (c == delimeter) return;
      if (c == EOF) {
        close();
        return;
      }
    }
  }
  void getline(std::string & line, const char delimeter = '\n') {
    if (file) {
      if (getdelim(&getline_line, &getline_n, delimeter, file) < 0) {
        close();
      } else {
        line = getline_line;  // Will not get embedded nulls
        if (line.back() == delimeter) line.pop_back();
      }
    }
  }

 private:
  template <class Value>
  void scanf(const char * const format, Value & value) {
    if (file) {
      if (peek() == EOF) {
        close();
      } else {
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wformat-nonliteral"
        if (fscanf(file, format, &value) == 0) close();
#pragma GCC diagnostic pop
      }
    }
  }
  FILE * file{nullptr};
  std::string name;
  bool complain{true};
  std::function<void(std::string &, FILE *)> closer{};
  char * getline_line{nullptr};
  size_t getline_n{0};
};
bool getline(ifile & in, std::string & line, const char delimeter = '\n') {
  in.getline(line, delimeter);
  return in;
}

// Output stream to a process
struct opopen : public ofile {
  // Constructors and destructors
  explicit opopen(const std::string & command, const bool complain_ = true) :
      ofile{"opopen", complain_} { open(command); }
  opopen(const opopen &) = delete;
  opopen & operator=(const opopen &) = delete;
  opopen(opopen && other) = default;
};

// Output stream from a process
struct ipopen : public ifile {
  // Constructors and destructors
  explicit ipopen(const std::string & command, const bool complain_ = true) :
      ifile{"ipopen", complain_} { open(command); }
  ipopen(const ipopen &) = delete;
  ipopen & operator=(const ipopen &) = delete;
  ipopen(ipopen && other) = default;
};
bool getline(ipopen & in, std::string & line) {
  in.getline(line);
  return in;
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wnon-virtual-dtor"
struct iopopen : private Fifo, public ipopen, public std::ofstream {
  explicit iopopen(const std::string & command_) :
      ipopen{"(" + command_ + ") < " + Fifo::name()},
      std::ofstream{Fifo::name()} {}
  ~iopopen() { close(); }
  iopopen(const iopopen &) = delete;
  iopopen & operator=(const iopopen &) = delete;
  iopopen(iopopen && other) :
      Fifo{std::move(other)}, ipopen{std::move(other)},
      std::ofstream{std::move(other)} {}
  ipopen & in() { return *this; }
  std::ofstream & out() { return *this; }
  void close() { std::ofstream::close(); }
};
#pragma GCC diagnostic pop


// These both do not work - mysterious hang occurs (deadlock?)
#if 0
struct iopopen {
  explicit iopopen(const std::string & command_) :
      command{"(" + command_ + ") < " + fifo.name()} {}
  void open_in() { in_ = std::make_unique<ipopen>(command); }
  void open_out() { out_ = std::make_unique<std::ofstream>(fifo.name()); }
  ipopen & in() const { return *in_; }
  std::ofstream & out() const { return *out_; }

 private:
  Fifo fifo{};
  std::unique_ptr<ipopen> in_{nullptr};
  std::unique_ptr<std::ofstream> out_{nullptr};
  std::string command;
};
struct oipopen {
  explicit oipopen(const std::string & command_) :
      command{"(" + command_ + ") > " + fifo.name()} {}
  void open_in() { in_ = std::make_unique<std::ifstream>(fifo.name()); }
  void open_out() { out_ = std::make_unique<opopen>(command); }
  std::ifstream & in() const { return *in_; }
  opopen & out() const { return *out_; }

 private:
  Fifo fifo{};
  std::unique_ptr<std::ifstream> in_{nullptr};
  std::unique_ptr<opopen> out_{nullptr};
  std::string command;
};
#endif

constexpr uint64_t InputEnd{1};
constexpr uint64_t OutputEnd{0};
struct opipe : public ofile {
  explicit opipe(const std::string & command, const bool complain_ = true) :
      ofile{"opipe", complain_} {
    static std::mutex mutex;
    std::lock_guard<std::mutex> lock(mutex);  // Necessary
    int out_fd[2]{0, 0};
    if (pipe(out_fd) < 0) throw Error("Pipe open error in opipe");
    pid = fork();
    if (pid == 0) {  // Child
      ::close(out_fd[InputEnd]);
      dup2(out_fd[OutputEnd], STDIN_FILENO);
      ::close(out_fd[OutputEnd]);
      execlp("/bin/bash", "/bin/bash", "-c", command.c_str(), nullptr);
      perror("opipe");  // Only reaches here on exec error
      _exit(1);
    } else {  // Parent
      ::close(out_fd[OutputEnd]);
      ofile::open(out_fd[InputEnd]);
    }
  }
  ~opipe() { close(); }
  opipe(const opipe &) = delete;
  opipe & operator=(const opipe &) = delete;
  opipe(opipe && other) : ofile{std::move(other)}, pid{other.pid} {
    other.pid = 0;
  }
  void close() {
    ofile::close();
    if (pid) {
      int status{0};
      pid = waitpid(pid, &status, 0);
      pid = 0;
    }
  }

 private:
  pid_t pid{0};
};

struct ipipe : public ifile {
  explicit ipipe(const std::string & command, const bool complain_ = true) :
      ifile{"ipipe", complain_} {
    static std::mutex mutex;
    std::lock_guard<std::mutex> lock(mutex);
    int in_fd[2]{0, 0};
    if (pipe(in_fd) < 0) throw Error("Pipe open error in ipipe");
    pid = fork();
    if (pid == 0) {
      // Child
      ::close(in_fd[OutputEnd]);
      dup2(in_fd[InputEnd], STDOUT_FILENO);
      ::close(in_fd[InputEnd]);
      execlp("/bin/bash", "/bin/bash", "-c", command.c_str(), nullptr);
      // Only reaches here on exec error
      perror("ipipe");
      _exit(1);
    } else {
      // Parent
      ::close(in_fd[InputEnd]);
      open(in_fd[OutputEnd]);
    }
  }
  ~ipipe() { close(); }
  ipipe(const ipipe &) = delete;
  ipipe & operator=(const ipipe &) = delete;
  ipipe(ipipe && other) : ifile{std::move(other)}, pid{other.pid} {
    other.pid = 0;
  }
  // operator bool() const { return file; }
  void close() {
    ifile::close();
    if (pid) {
      int status{0};
      pid = waitpid(pid, &status, 0);
      pid = 0;
    }
  }

 private:
  pid_t pid{0};
};
ipipe & getline(ipipe & in, std::string & line, const char delimeter = '\n') {
  in.getline(line, delimeter);
  return in;
}
#if 0
inline std::string getline(ipipe & input) {
  std::string line;
  getline(input, line);
  return line;
}
#endif
struct iopipe : public ifile, public ofile {
  explicit iopipe(const std::string & command, const bool complain_ = true) :
      ifile{"iopipe", complain_}, ofile{"iopipe", complain_} {
    static std::mutex mutex;
    std::lock_guard<std::mutex> lock(mutex);
    int in_fd[2]{0, 0};
    if (pipe(in_fd) < 0) throw Error("Pipe in open error in iopipe");
    int out_fd[2]{0, 0};
    if (pipe(out_fd) < 0) throw Error("Pipe out open error in iopipe");
    pid = fork();
    if (pid == 0) {
      // Child
      ::close(in_fd[OutputEnd]);
      dup2(in_fd[InputEnd], STDOUT_FILENO);
      ::close(in_fd[InputEnd]);
      ::close(out_fd[InputEnd]);
      dup2(out_fd[OutputEnd], STDIN_FILENO);
      ::close(out_fd[OutputEnd]);
      execlp("/bin/bash", "/bin/bash", "-c", command.c_str(), nullptr);
      // Only reaches here on exec error
      perror("iopipe");
      _exit(1);
    } else {
      // Parent
      ::close(in_fd[InputEnd]);
      ::close(out_fd[OutputEnd]);
      ifile::open(in_fd[OutputEnd]);
      ofile::open(out_fd[InputEnd]);
    }
  }
  ifile & in() { return *this; }
  ofile & out() { return *this; }
  ~iopipe() {
    close();
    ifile::close();
  }
  iopipe(const iopipe &) = delete;
  iopipe & operator=(const iopipe &) = delete;
  iopipe(iopipe && other) :
      ifile{std::move(other)}, ofile{std::move(other)}, pid{other.pid} {
    other.pid = 0;
  }
  operator bool() const { return ifile::operator bool(); }
  void close_out() {
    ofile::close();
  }
  void close() {
    ofile::close();
    if (pid) {
      int status{0};
      pid = waitpid(pid, &status, 0);
      pid = 0;
    }
  }

 private:
  pid_t pid{0};
};

}  // namespace paa

#endif  // PAA_EMPTY_H_
