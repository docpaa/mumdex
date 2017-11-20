//
// sequence.h
//
// various sequence classes and functions
//
// Copyright 2016 Peter Andrews
//

#ifndef PAA_SEQUENCE_H
#define PAA_SEQUENCE_H

#include <algorithm>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "pstream.h"

namespace paa {

class Similarity {
 public:
  Similarity(const std::string & seq1, const std::string & seq2) {
    if (seq1.size() != seq2.size())
      throw Error("Bad sequence size in Similarity");
    for (unsigned int b{0}; b != seq1.size(); ++b) {
      if (seq1[b] != 'x' && seq1[b] != 'N' &&
          seq1[b] != 'x' && seq2[b] != 'N') {
        if (seq1[b] == seq2[b]) {
          ++n_matches_;
        } else {
          ++n_mismatches_;
        }
      }
    }
  }

  unsigned int n_matches() const { return n_matches_; }
  unsigned int n_testable() const { return n_matches_ + n_mismatches_; }
  unsigned int n_mismatches() const { return n_mismatches_; }

 private:
  unsigned int n_matches_{0};
  unsigned int n_mismatches_{0};
};

class ConsensusSequence {
 public:
  class BaseId {
   public:
    explicit BaseId(const char base) {
      switch (base) {
        case 'A':
          index = 0;
          break;
        case 'C':
          index = 1;
          break;
        case 'G':
          index = 2;
          break;
        case 'T':
          index = 3;
          break;
        case 'N':
          index = 4;
          break;
        case 'x':
          index = 5;
          break;
        default:
          throw Error("Unexpected base seen") << base;
      }
    }
    operator unsigned int() const {
      return index;
    }
    unsigned int value() const {
      return index;
    }

   private:
    unsigned int index{};
  };

  class Base {
   public:
    explicit Base(const unsigned int index) {
      switch (index) {
        case 0:
          base = 'A';
          break;
        case 1:
          base = 'C';
          break;
        case 2:
          base = 'G';
          break;
        case 3:
          base = 'T';
          break;
        case 4:
          base = 'N';
          break;
        case 5:
          base = 'x';
          break;
        default:
          throw Error("Unexpected index seen") << index;
      }
    }
    operator char() const {
      return base;
    }
    char value() const {
      return base;
    }

    static constexpr unsigned int n_bases{6};
    static constexpr unsigned int n_good_bases{4};

   private:
    char base{};
  };

  ConsensusSequence(const unsigned int length_arg) :  // NOLINT
      base_counts(length_arg, std::array<unsigned int, Base::n_bases>{
          {0, 0, 0, 0, 0, 0}}) { }
  void add(const std::string & seq) {
    if (seq.size() != base_counts.size())
      throw Error("Bad sequence size in ConsensusSequence::add");
    for (unsigned int b{0}; b != seq.size(); ++b) {
      ++base_counts[b][BaseId{seq[b]}];
    }
  }
  // Need to call sequences() before similarity() will work
  const std::string & sequence() const {
    determine();
    return consensus;
  }

  Similarity similarity(const std::string & seq) const {
    return Similarity{seq, consensus};
  }
  unsigned int size() const {
    return static_cast<unsigned int>(base_counts.size());
  }
  unsigned int count(const unsigned int b, const char base_) const {
    return base_counts[b][BaseId{base_}];
  }
  unsigned int count(const unsigned int b, const unsigned int index) const {
    return base_counts[b][index];
  }

  void determine() const {
    consensus.clear();
    for (const auto counts : base_counts) {
      using entry = std::pair<unsigned int, unsigned int>;
      std::array<entry, Base::n_good_bases> sorted;
      for (unsigned int b{0}; b != Base::n_good_bases; ++b) {
        sorted[b] = {b, counts[b]};
      }
      sort(sorted.begin(), sorted.end(),
           [](const entry left, const entry right) {
             return left.second > right.second;
           });
      const char base_{[&sorted] {
          if (sorted.front().first < Base::n_good_bases &&
              sorted.front().second &&
              sorted[0].second != sorted[1].second) {
            return static_cast<char>(Base{sorted.front().first});
          } else {
            return 'x';
          }
        }()};
      consensus += base_;
    }
  }

 private:
  std::vector<std::array<unsigned int, Base::n_bases>> base_counts{};
  mutable std::string consensus{};
};
using Base = ConsensusSequence::Base;

class Repeatness {
 public:
  template <class longSA>
  explicit Repeatness(const std::string & sequence_arg, const longSA & sa,
                      const Reference & ref, const Mappability & mappability) :
      sequence_{sequence_arg}, mismatches_(4) {
    const auto mams = sa.find_mams(sequence_);
    n_mams_ = static_cast<unsigned int>(mams.size());
    for (const auto mam : mams) {
      if (mam.len > max_length_) max_length_ = mam.len;
    }
    for (const auto mam : mams) {
      const unsigned int chr{mam.chr};
      const unsigned int pos{mam.pos};
      const unsigned int abspos_start{ref.offset(chr) + pos};
      for (unsigned int b{0}; b != mam.len; ++b) {
        const unsigned int low_map{mappability.low(abspos_start + b)};
        if (low_map <= b + 1) {
          if (min_mappability_ > mam.len - low_map) {
            min_mappability_ = mam.len - low_map;
          }
          if (max_mappability_ < mam.len - low_map) {
            max_mappability_ = mam.len - low_map;
          }
        }
        const unsigned int high_map{mappability.high(abspos_start + b)};
        if (high_map <= mam.len - b) {
          if (min_mappability_ > mam.len - high_map) {
            min_mappability_ = mam.len - high_map;
          }
          if (max_mappability_ < mam.len - high_map) {
            max_mappability_ = mam.len - high_map;
          }
        }
      }
    }
    const std::string bowtie_command{
      "/data/software/bowtie/bowtie-0.12.8/bowtie "
          "--mm -a --best -v 3 hg19 --suppress 1,2,3,4,5,6,7 --quiet -c "};
    const std::string filter{"| perl -pe 's/,/ /g'"};
    redi::ipstream bowtie{bowtie_command + sequence_ + filter};
    std::string alignment;
    while (getline(bowtie, alignment)) {
      std::istringstream line_to_parse{alignment};
      std::string mismatch;
      unsigned int n_mismatches{0};
      while (line_to_parse >> mismatch) {
        ++n_mismatches;
      }
      ++mismatches_[n_mismatches];
    }
  }

  const std::string & sequence() const {
    return sequence_;
  }
  unsigned int n_mams() const {
    return n_mams_;
  }
  const std::vector<unsigned int> & mismatches() const {
    return mismatches_;
  }
  unsigned int max_length() const {
    return max_length_;
  }
  unsigned int min_mappability() const {
    return min_mappability_;
  }
  unsigned int max_mappability() const {
    return max_mappability_;
  }

  std::string summary() const {
    std::ostringstream out;
    out << "sequence " << sequence() << " of length " << sequence().size()
        << " has " << n_mams() << " mum" << (n_mams() == 1 ? "" : "s");
    if (n_mams()) {
      out << " with max length " << max_length()
          << " and excess mappability "
          << "min: " << min_mappability() << ", "
          << "max: " << max_mappability();
    }
    const unsigned int alignments{mismatches()[0] + mismatches()[1] +
          mismatches()[2] + mismatches()[3]};
    out << " and "
        << alignments
        << " bowtie alignments";
    if (alignments) {
      out << " with mismatches "
          << 0 << ":" << mismatches()[0] << ", "
          << 1 << ":" << mismatches()[1] << ", "
          << 2 << ":" << mismatches()[2] << ", "
          << 3 << ":" <<  mismatches()[3];
    }
    return out.str();
  }

  std::string minimal() const {
    std::ostringstream out;
    out << n_mams() << " ";
    if (n_mams()) {
      out << max_mappability();
    } else {
      out << 0;
    }
    out << " " << mismatches()[0]
        << " " << mismatches()[1]
        << " " << mismatches()[2]
        << " " << mismatches()[3];
    return out.str();
  }

 private:
  std::string sequence_{};
  unsigned int n_mams_{};
  std::vector<unsigned int> mismatches_{};
  unsigned int max_length_{0};
  unsigned int min_mappability_{255};
  unsigned int max_mappability_{0};
};

class GCATscore {
 public:
  explicit GCATscore(const std::string & sequence) {
    for (const char base : sequence) {
      switch (base) {
        case 'A':
        case 'T':
        case 'N':
            score += 2;
        break;
        case 'C':
        case 'G':
            score += 4;
        break;
        default:
          throw Error("Unknown base seen") << base;
      }
    }
  }
  operator unsigned int() const { return score; }
  unsigned int value() const { return score; }

 private:
  unsigned int score{0};
};


}  // namespace paa

#endif  // PAA_SEQUENCE_H
