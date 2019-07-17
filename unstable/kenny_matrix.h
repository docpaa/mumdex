//
// kenny_matrix.h
//
// simple matrix computations
//
// Copyright 2016 Peter Andrews @ CSHL
//
//

#ifndef PAA_KENNY_MATRIX_H
#define PAA_KENNY_MATRIX_H

#include <algorithm>
#include <future>
#include <iostream>
#include <memory>
#include <fstream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "error.h"
#include "files.h"
#include "threads.h"
#include "utility.h"

namespace paa {

template <class Type>
class DiagMatrix {
 public:
  DiagMatrix(DiagMatrix &&) = default;
  DiagMatrix(const DiagMatrix &) = delete;
  DiagMatrix & operator=(const DiagMatrix &) = delete;
  DiagMatrix & operator=(DiagMatrix &&) = delete;

  explicit DiagMatrix(const uint64_t n__) : n_{n__}, data_(n_) { }

  explicit DiagMatrix(const std::string file_name) {
    const uint64_t found{file_name.find("_matrix")};
    if (found == std::string::npos) {
      throw Error("Filename bad format in DiagMatrix") << file_name;
    }
    const std::string matrix_string{file_name.substr(found)};
    std::istringstream file_stream{matrix_string.c_str()};
    std::string extension;
    char sep;
    file_stream.ignore(10000, '.');
    file_stream >> n_ >> sep >> extension;
    if (extension == "txt") {
      data_.clear_and_resize(n_);
      read(file_name);
    } else if (extension == "bin") {
      data_.read_mapped(file_name);
    } else {
      throw Error("Bad diag matrix file name extension") << extension;
    }
  }

  void save(const std::string & base_name) const {
    std::ostringstream file_name;
    file_name << base_name << "_matrix." << n_ << ".txt";
    std::ofstream file{file_name.str().c_str()};
    if (!file) throw Error("Problem opening matrix file for writing")
                   << file_name;
    output(file);
  }
  void read(const std::string & file_name) {
    std::ifstream file{file_name.c_str()};
    if (!file) throw Error("Problem opening matrix file for reading")
                   << file_name;
    for (uint64_t x{0}; x != n_; ++x) {
      file >> (*this)(x);
    }
    if (!file) throw Error("DiagMatrix parse error") << file_name;
  }
  void write(const std::string & base_name) const {
    std::ostringstream file_name;
    file_name << base_name << "_matrix." << n_ << ".bin";
    data_.save(file_name.str());
  }
  std::ostream & output(std::ostream & out) const {
    for (uint64_t x{0}; x != n_; ++x) {
      if (x) out << '\t';
      out << (*this)(x);
    }
    return out << '\n';
  }
  uint64_t nx() const { return n_; }
  uint64_t ny() const { return n_; }

  // Type operator()(const uint64_t x, const uint64_t y) const {
  //  return x == y ? data_[x] : Type(0);
  // }
  Type & operator()(const uint64_t x) {
    return data_[x];
  }
  Type operator()(const uint64_t x) const {
    return data_[x];
  }

 private:
  uint64_t n_{};
  FlexVector<Type> data_{};
};

template <class Type>
inline std::ostream & operator<<(std::ostream & out,
                                 const DiagMatrix<Type> & matrix) {
  return matrix.output(out);
}

template <class Type>
class BMatrix {
 public:
  BMatrix(BMatrix &&) = default;
  BMatrix(const BMatrix &) = delete;
  BMatrix & operator=(const BMatrix &) = delete;
  BMatrix & operator=(BMatrix &&) = delete;

  BMatrix(const uint64_t nx__, const uint64_t ny__, bool row_major_ = true) :
      nx_{nx__}, ny_{ny__}, data_(nx_ * ny_), row_major{row_major_} { }
  BMatrix(const uint64_t nx__, const uint64_t ny__, const Type val,
         bool row_major_ = true) :
      BMatrix{nx__, ny__, row_major_} {
    set(val);
  }
  explicit BMatrix(const std::string file_name, bool row_major_ = true) :
      row_major{row_major_} {
    const uint64_t found{file_name.find("_matrix")};
    if (found == std::string::npos) {
      throw Error("Filename bad format in BMatrix") << file_name;
    }
    const std::string matrix_string{file_name.substr(found)};
    std::istringstream file_stream{matrix_string.c_str()};
    std::string extension;
    char sep;
    file_stream.ignore(10000, '.');
    file_stream >> nx_ >> sep >> ny_ >> sep >> extension;
    if (extension == "txt") {
      data_.clear_and_resize(nx_ * ny_);
      read(file_name);
    } else if (extension == "bin") {
      data_.read_mapped(file_name, true);
    } else {
      throw Error("Bad matrix file name extension") << extension;
    }
  }
  // Create from rows of other matrix
  BMatrix(const BMatrix & other, const uint64_t y_start,
          const uint64_t n_rows) :
      BMatrix{other.nx_, n_rows} {
    for (uint64_t y{0}; y != ny_; ++y) {
      for (uint64_t x{0}; x != nx_; ++x) {
        (*this)(x, y) = other(x, y + y_start);
      }
    }
  }
  ~BMatrix() { }

  void save(const std::string & base_name) const {
    std::ostringstream file_name;
    file_name << base_name << "_matrix." << nx_ << "." << ny_ << ".txt";
    std::ofstream file{file_name.str().c_str()};
    if (!file) throw Error("Problem opening matrix file for writing")
                   << file_name.str();
    output(file);
  }
  std::ostream & output(std::ostream & out) const {
    for (uint64_t y{0}; y != ny_; ++y) {
      for (uint64_t x{0}; x != nx_; ++x) {
        if (x) out << '\t';
        out << (*this)(x, y);
      }
      out << '\n';
    }
    return out;
  }
  void read(const std::string & file_name) {
    std::ostringstream message;
    message << "Reading BMatrix " << nx_ << " x " << ny_;
    Progress progress(ny_, 0.1, message.str());
    std::ifstream file{file_name.c_str()};
    if (!file) throw Error("Problem opening matrix file for reading")
                   << file_name;
    for (uint64_t y{0}; y != ny_; ++y) {
      progress();
      for (uint64_t x{0}; x != nx_; ++x) {
        file >> (*this)(x, y);
      }
    }
    if (!file) throw Error("BMatrix parse error") << file_name;
  }
  std::string write(const std::string & base_name) const {
    std::ostringstream file_name;
    file_name << base_name << "_matrix." << nx_ << "." << ny_ << ".bin";
    data_.save(file_name.str());
    return file_name.str();
  }

  uint64_t nx() const { return nx_; }
  uint64_t ny() const { return ny_; }

  Type & operator()(const uint64_t x, const uint64_t y) {
    return row_major ? data_[x + y * nx_] : data_[x * ny_ + y];
  }
  Type operator()(const uint64_t x, const uint64_t y) const {
    return row_major ? data_[x + y * nx_] : data_[x * ny_ + y];
  }

  std::vector<std::pair<std::string, unsigned int> >
  write_rows(const std::string & base) const {
    std::vector<std::pair<std::string, unsigned int>> result;
    const uint64_t n_rows_goal{30};
    for (uint64_t y{0}; y < ny_; y += n_rows_goal) {
      const uint64_t n_rows{std::min(n_rows_goal, ny_ - y)};
      std::ostringstream out_name;
      out_name << base << "." << y;
      const BMatrix rows{*this, y, n_rows};
      result.push_back(std::pair<std::string, unsigned int>(
          rows.write(out_name.str()), n_rows));
    }
    return result;
  }
  BMatrix dist_mul(const BMatrix & rhs, const std::string & out_dir) const {
    throw Error("Set up distributed runs with Peter first");
    std::cerr << "Writing matrices for distributed computation" << std::endl;
    const std::string rhs_name{rhs.write("rhs")};
    const std::vector<std::pair<std::string, unsigned int>>
        names{write_rows("lhs")};
    const unsigned int n_threads{8};
    const double memory{std::max(
        1.0 * sizeof(Type) *
        (rhs.nx_ * ny_) / 1024 / 1024 / 1024,
        static_cast<double>(n_threads)) / n_threads};
    std::vector<std::string> result_names;
    Progress progress{names.size(), 0.01, "Distributed run"};
    std::ofstream sge_commands{"sge_commands.txt"};
    std::cerr << "Running jobs" << std::endl;
    for (unsigned int j{0}; j != names.size(); ++j) {
      std::ostringstream job_name;
      job_name << "mm" << j;
      const std::string jn{job_name.str()};
      const std::string ob{out_dir + "/" + jn};
      std::ostringstream job_stream;
      job_stream << "qsub -l vf=" << memory << "G -N " << jn
                 << " -pe threads " << n_threads << " -wd " << out_dir
                 << " -o " << ob << ".out.txt -e " << ob << ".err.txt"
                 << " -shell no -b yes /data/unsafe/paa/mums/mumdex/matrix row "
                 << names[j].first << " " << rhs_name
                 << " " << ob << " dummy " << n_threads << " > /dev/null";
      sge_commands << "qdel -j " << jn << "; " << job_stream.str() << std::endl;
      if (system(job_stream.str().c_str()) == -1) {
        std::cerr << "Problem running jobs" << std::endl;
      }
      std::ostringstream result_name;
      result_name << ob << "_matrix."
                  << rhs.nx_ << "." << names[j].second << ".bin";
      result_names.push_back(result_name.str());
    }
    for (unsigned int j{0}; j != result_names.size(); ++j) {
      while (!readable(result_names[j])) {
        sleep(1);
      }
      progress();
    }
    BMatrix result{rhs.nx_, ny_};
    for (unsigned int j{0}; j != names.size(); ++j) {
      const BMatrix input{result_names[j]};
      for (unsigned int x{0}; x != result.nx_; ++x) {
        result(x, j) = input(x, 0);
      }
    }
    return result;
  }

  std::string info() const {
    std::ostringstream out;
    out << nx_ << " x " << ny_ << " " << data_.size();
    return out.str();
  }

  // Slow algorithm
  BMatrix operator*(const BMatrix & rhs) const {
    if (nx_ != rhs.ny_) {
      throw Error("BMatrix dimensions bad for multiplication")
          << "lhs nx" << nx_ << "!=" << "rhs ny" << rhs.ny_;
    }
    BMatrix result(rhs.nx_, ny_, 0.0);
    Progress progress{rhs.nx_ * ny_, 0.01, "BMatrix multiplication"};
    for (uint64_t x2{0}; x2 != rhs.nx_; ++x2) {
      for (uint64_t y1{0}; y1 != ny_; ++y1) {
        progress();
        for (uint64_t x1{0}; x1 != nx_; ++x1) {
          result(x2, y1) += (*this)(x1, y1) * rhs(x2, x1);
        }
      }
    }
    return result;
  }

  void set(const Type val) {
    if (row_major) {
      for (unsigned int y{0}; y != ny_; ++y) {
        for (unsigned int x{0}; x != nx_; ++x) {
          (*this)(x, y) = val;
        }
      }
    } else {
      for (unsigned int x{0}; x != nx_; ++x) {
        for (unsigned int y{0}; y != ny_; ++y) {
          (*this)(x, y) = val;
        }
      }
    }
  }

  // Parallel algorithm
  BMatrix thread_mul(const BMatrix & rhs, const unsigned int n_threads) const {
    if (nx_ != rhs.ny_) {
      throw Error("BMatrix dimensions bad for multiplication")
          << "lhs nx" << nx_ << "!=" << "rhs ny" << rhs.ny_;
    }
    if (rhs.nx_ == 1 || n_threads == 1) { return (*this) * rhs; }
    BMatrix result(rhs.nx_, ny_, 0.0);
    ThreadPool pool{n_threads};
    ThreadPool::Results<void> results;
    for (uint64_t x2{0}; x2 != rhs.nx_; ++x2) {
      pool.run(results, [this, x2, &rhs](BMatrix & result_) {
          for (uint64_t y1{0}; y1 != ny_; ++y1) {
            for (uint64_t x1{0}; x1 != nx_; ++x1) {
              result_(x2, y1) += (*this)(x1, y1) * rhs(x2, x1);
            }
          }}, std::ref(result));
    }

    Progress progress{rhs.nx_, 0.01, "BMatrix multiplication"};
    while (results.size()) {
      results.get();
      progress();
    }
    return result;
  }

  BMatrix operator*(const DiagMatrix<Type> & rhs) const {
    if (nx_ != rhs.ny()) {
      throw Error("BMatrix dimensions bad for multiplication")
          << "lhs nx" << nx_ << "!=" << "rhs ny" << rhs.ny();
    }
    BMatrix result(rhs.nx(), ny_, 0.0);
    Progress progress{ny_, 0.1, "Diagonal BMatrix multiplication"};
    for (uint64_t y1{0}; y1 != ny_; ++y1) {
      progress();
      for (uint64_t x2{0}; x2 != rhs.nx(); ++x2) {
        result(x2, y1) += (*this)(x2, y1) * rhs(x2);
      }
    }
    return result;
  }

  BMatrix diag_mul(const DiagMatrix<Type> & rhs,
                  const unsigned int n_threads) const {
    if (nx_ != rhs.ny()) {
      throw Error("BMatrix dimensions bad for multiplication")
          << "lhs nx" << nx_ << "!=" << "rhs ny" << rhs.ny();
    }
    BMatrix result(rhs.nx(), ny_, 0.0);
    ThreadPool pool{n_threads};
    ThreadPool::Results<void> results;
    Progress progress{ny_, 0.1, "Diagonal BMatrix multiplication"};
    for (uint64_t y1{0}; y1 != ny_; ++y1) {
      pool.run(results, [this, y1, &rhs](BMatrix & result_) {
          for (uint64_t x2{0}; x2 != rhs.nx(); ++x2) {
            result_(x2, y1) += (*this)(x2, y1) * rhs(x2);
          }
        }, std::ref(result));
    }
    while (results.size()) {
      results.get();
      progress();
    }
    return result;
  }

 private:
  Type * data() { return &data_[0]; }
  Type * row_start(const uint64_t y) { return &data_[y * nx_]; }
  Type * col_start(const uint64_t x) { return &data_[x]; }

  uint64_t nx_{};
  uint64_t ny_{};
  FlexVector<Type> data_{};
  bool row_major{true};
};

template <class Type>
inline std::ostream & operator<<(std::ostream & out,
                                 const BMatrix<Type> & matrix) {
  return matrix.output(out);
}

}  // namespace paa

#endif  // PAA_KENNY_MATRIX_H
