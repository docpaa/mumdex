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
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "error.h"
#include "mumdex.h"
#include "named_ints.h"

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

namespace paa {

const uint64_t good_alignment_score{40};
const uint64_t good_base_quality{30};
const double good_quality_fraction{0.666};

struct Sam {
  Sam() = default;
  template<class Stream> explicit Sam(Stream & sam) : Sam{getline(sam)} { }
  explicit Sam(const std::string & sam_line) {
    std::istringstream line_to_parse{sam_line};
    std::string name_plus;
    line_to_parse >> name_plus;
    std::istringstream name_stream{name_plus.c_str()};
    getline(name_stream, name, '_');
    name_stream >> read_conversion;
    if (!line_to_parse || name.empty())
      throw Error("sam empty line") << sam_line;
    line_to_parse >> flag >> chr >> pos >> score >> cigar
                  >> mchr >> mpos >> tlen >> bases >> qual;
    if (!line_to_parse)
      throw Error("sam parse error") << ('"' + sam_line + '"');
    line_to_parse.get();
    getline(line_to_parse, optional);
    // cerr << "optional (" << optional << ")" << endl;
    std::istringstream optional_stream{optional.c_str()};
    std::string opt;
    while (optional_stream >> opt) {
      std::istringstream optional_part{opt};
      std::string oname;
      getline(optional_part, oname, ':');
      std::string val_type;
      getline(optional_part, val_type, ':');
      if (val_type == "i") {
        int val;
        optional_part >> val;
        if (oname == "AS") {
          alignment_score = val;
        } else if (oname == "XS") {
          next_alignment_score = val;
        } else if (oname == "NM") {
          edit_distance = val;
        } else if (oname == "XM") {
          n_mismatch = val;
        } else if (oname == "XO") {
          n_gap_open = val;
        } else if (oname == "XG") {
          n_gap_extension = val;
        }
      } else if (val_type == "Z") {
        std::string val;
        optional_part >> val;
        if (oname == "SB") {
          sample_barcode = val;
        } else if (oname == "GC") {
          genome_conversion = val;
        } else if (oname == "RC") {
          read_conversion = val;
        }
      }
    }
  }
  bool is_reversed() const { return flag & sam::is_reversed; }
  bool is_mapped() const { return pos && (flag & sam::is_unmapped) == 0; }
  bool is_unique() const {
    return is_mapped() &&
        (next_alignment_score == 10000 ||
         next_alignment_score < alignment_score);
  }
  bool is_proper() const {
    return flag & sam::is_proper;
  }
  void replace_sequence(std::string new_sequence) {
    if (false && bases.size() != new_sequence.size())
      throw Error("Sequence size mismatch in replace_sequence")
          << bases.size() << new_sequence.size();
    if (is_reversed()) new_sequence = reverse_complement(new_sequence);
    bases = new_sequence;
  }
  bool good_base_qual() const {
    // Are most qualities good?
    const char good_qual{33 + good_base_quality};
    uint64_t n_good{0};
    for (const char q : qual) n_good += q >= good_qual;
    return n_good > good_quality_fraction * qual.size();
  }
  void output(std::ostream & out) const {
    out << name
        << '\t' << flag
        << '\t' << chr
        << '\t' << pos
        << '\t' << score
        << '\t' << cigar
        << '\t' << mchr
        << '\t' << mpos
        << '\t' << tlen
        << '\t' << bases
        << '\t' << qual;
    if (sample_barcode.size())
      out << '\t' << "SB:Z:" << sample_barcode;
    if (template_barcode.size())
      out << '\t' << "TB:Z:" << template_barcode;
    if (linear_barcode.size())
      out << '\t' << "LB:Z:" << linear_barcode;
    if (optional.size())
      out << '\t' << optional;
    out << '\t' << "RC:Z:" << read_conversion;
  }
  std::string name{};
  std::string read_conversion{"*"};
  std::string genome_conversion{"*"};
  std::string sample_barcode{""};
  std::string template_barcode{""};
  std::string linear_barcode{""};
  uint64_t flag{0};
  std::string chr{"*"};
  unsigned int pos{0};
  unsigned int score{0};
  std::string cigar{"*"};
  std::string mchr{"*"};
  unsigned int mpos{0};
  int tlen{0};
  std::string bases{};
  std::string qual{};
  std::string optional{};
  int alignment_score{10000};
  int next_alignment_score{0};
  int edit_distance{0};
  int n_mismatch{0};
  int n_gap_open{0};
  int n_gap_extension{0};
};
std::ostream & operator<<(std::ostream & out, const Sam & sam) {
  sam.output(out);
  return out;
}

}  // namespace paa

#endif  // PAA_SAM_H

