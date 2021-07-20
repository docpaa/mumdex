//
// sam.h
//
// sam info
//
// Copyright 2016 Peter Andrews @ CSHL
//

#ifndef PAA_SAM_H
#define PAA_SAM_H

#include <algorithm>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "error.h"

namespace sam {

enum MapFlag {
  is_paired = 1 << 0,
  is_proper = 1 << 1,
  is_unmapped = 1 << 2,
  is_mate_unmapped = 1 << 3,
  is_reversed = 1 << 4,
  is_mate_reversed = 1 << 5,
  is_first = 1 << 6,
  is_second = 1 << 7,
  is_not_primary = 1 << 8,
  is_bad_vendor_quality = 1 << 9,
  is_a_duplicate = 1 << 10,
  is_supplemental = 1 << 11
};

const std::vector<std::string> descriptions = {
  "paired",
  "proper",
  "unmapped",
  "mate_unmapped",
  "reversed",
  "mate_reversed",
  "first",
  "second",
  "not_primary",
  "bad_vendor_quality",
  "a_duplicate",
  "supplementary"
};

std::string flag_text(unsigned int flag) {
  std::string output;
  unsigned int index = 0;
  while (flag) {
    if (flag & 1) {
      if (output.size()) output += ' ';
      output += descriptions[index];
    }
    ++index;
    flag >>= 1;
  }
  return output;
}

using RefPos = std::pair<unsigned int, bool>;
RefPos get_ref_pos(const std::string & cigar, const unsigned int read_pos,
                   const unsigned int ref_pos_start) {
  if (ref_pos_start == 0 || cigar == "*") return RefPos{0, false};
  std::istringstream cigar_stream{cigar.c_str()};
  unsigned int cigar_bases;
  char cigar_char;
  unsigned int base{0};
  unsigned int pos = ref_pos_start;
  while (cigar_stream >> cigar_bases >> cigar_char) {
    switch (cigar_char) {
      case 'S':
        {
          base += cigar_bases;
          if (read_pos < base) return RefPos{pos, false};
          break;
        }
      case 'M':
      case '=':
      case 'X':
        {
          if (read_pos < base + cigar_bases)
            return RefPos{pos + read_pos - base, true};
          base += cigar_bases;
          pos += cigar_bases;
          break;
        }
      case 'D':
      case 'N':
        {
          pos += cigar_bases;
          break;
        }
      case 'I':
        {
          if (read_pos < base + cigar_bases)
            return RefPos{pos, false};
          base += cigar_bases;
          break;
        }
      case 'H':
        {
          break;
        }
      default:
        {
          throw paa::Error("Unknown cigar character") << cigar_char;
        }
    }
  }
  throw paa::Error("Reached end of cigar string in get_ref_pos")
      << cigar << ref_pos_start << read_pos;
}

void fix_clipped(unsigned int n_clip,
                 std::string & cigar, unsigned int & ref_pos) {
  if (ref_pos == 0 || cigar == "*") return;
  std::istringstream cigar_stream{cigar.c_str()};
  std::ostringstream new_cigar;
  unsigned int cigar_bases;
  char cigar_char;
  while (cigar_stream >> cigar_bases >> cigar_char) {
    const unsigned int to_clip{std::min(cigar_bases, n_clip)};
    switch (cigar_char) {
      case 'S':
        {
          cigar_bases -= to_clip;
          n_clip -= to_clip;
          break;
        }
      case 'M':
      case '=':
      case 'X':
        {
          cigar_bases -= to_clip;
          n_clip -= to_clip;
          ref_pos += to_clip;
          break;
        }
      case 'D':
      case 'N':
        {
          if (to_clip) {
            ref_pos += cigar_bases;
            cigar_bases = 0;
          }
          break;
        }
      case 'I':
        {
          cigar_bases -= to_clip;
          n_clip -= to_clip;
          break;
        }
      case 'H':
        {
          break;
        }
      default:
        {
          throw paa::Error("Unknown cigar character") << cigar_char;
        }
    }
    if (cigar_bases) new_cigar << cigar_bases << cigar_char;
  }
  cigar = new_cigar.str();
}

}  // namespace sam

#endif  // PAA_SAM_H

