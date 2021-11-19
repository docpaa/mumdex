//
// named_ints.h
//
// Small integers named and type safe with other properties too
//
// copyright 2021 Peter Andrews
//

#ifndef PAA_NAMED_INTS_H_
#define PAA_NAMED_INTS_H_

#include <array>
#include <string>

#include "error.h"

namespace paa {

// Named bases
struct BaseType {
  const char * name() const {
    switch (index) {
      case 0:
        return "A";
      case 1:
        return "C";
      case 2:
        return "G";
      case 3:
        return "T";
      case 4:
        return "N";
      default:
        throw Error("Bad BaseType value") << index;
    }
  }
  const char * color() const {
    switch (index) {
      case 0:
        return "0.243 0.188 1";
      case 1:
        return "0.839 0 0.024";
      case 2:
        return "0.09 0.722 0.12";
      case 3:
        return "0.882 0.894 0.031";
      case 4:
        return "0 0 0";
      default:
        throw Error("Bad PanelCountType value") << index;
    }
  }
  uint64_t index;
  constexpr operator uint64_t() const { return index; }
  inline static const BaseType & char2base(const char c);
};
constexpr BaseType A{0};
constexpr BaseType C{1};
constexpr BaseType G{2};
constexpr BaseType T{3};
constexpr BaseType N{4};
constexpr BaseType ACGTN{5};
constexpr std::array<BaseType, ACGTN> acgtn_types{{A, C, G, T, N}};
const BaseType & BaseType::char2base(const char c) {
  switch (c) {
    case 'A':
      return A;
    case 'C':
      return C;
    case 'G':
      return G;
    case 'T':
      return T;
    case 'N':
      return N;
    default:
      throw Error("Unknown base") << c;
  }
}

// Named sequence / quality
struct SeqQual {
  const char * name() const {
    switch (index) {
      case 0:
        return "Sequence";
      case 1:
        return "Quality";
      case 2:
        return "Pattern";
      default:
        throw Error("Bad SeqQual value") << index;
    }
  }
  uint64_t index;
  constexpr operator uint64_t() const { return index; }
};
constexpr SeqQual Seq{0};
constexpr SeqQual Qual{1};
constexpr SeqQual Pattern{2};
constexpr SeqQual SQTs{2};
constexpr SeqQual SQPTs{3};
constexpr std::array<SeqQual, SQTs> sq_types{{Seq, Qual}};

// Named read indices
struct Read {
  std::string sname() const { return name(); }
  const char * name() const {
    switch (index) {
      case 0:
        return "Read 1";
      case 1:
        return "Read 2";
      case 2:
        return "Pair";
      case 3:
        return "Agree";
      default:
        throw Error("Bad Read value") << index;
    }
  }
  uint64_t index;
  constexpr operator uint64_t() const { return index; }
};
constexpr Read R1{0};
constexpr Read R2{1};
constexpr Read RP{2};
constexpr Read Agree{3};
constexpr Read RRP{3};
constexpr Read RRPA{4};
constexpr std::array<Read, RP> r12_types{{R1, R2}};

// Named begin / end limits
struct BeginEnd {
  const char * name() const {
    switch (index) {
      case 0:
        return "Begin";
      case 1:
        return "End";
      default:
        throw Error("Bad BeginEnd value") << index;
    }
  }
  uint64_t index;
  constexpr operator uint64_t() const { return index; }
};
constexpr BeginEnd B{0};
constexpr BeginEnd E{1};
constexpr BeginEnd BETs{2};
constexpr std::array<BeginEnd, BETs> be_types{{B, E}};

// Named parts of a locus
struct Part {
  uint64_t index;
  const char * name() const {
    switch (index) {
      case 0:
        return "Left";
      case 1:
        return "Right";
      case 2:
        return "MS";
      default:
        throw Error("Bad Part value") << index;
    }
  }
  const char * short_name() const {
    switch (index) {
      case 0:
        return "L";
      case 1:
        return "R";
      case 2:
        return "MS";
      default:
        throw Error("Bad Part value") << index;
    }
  }
  const char * color() const {
    switch (index) {
      case 0:
        return "0 0.8 0";
      case 1:
        return "1 0 0";
      case 2:
        return "0 0 0";
      default:
        throw Error("Bad Part value") << index;
    }
  }
  constexpr operator uint64_t() const { return index; }
};
constexpr Part L{0};
constexpr Part R{1};
constexpr Part LRTs{2};
constexpr Part MS{2};
constexpr Part LRMTs{3};
constexpr Part UNSET{3};
constexpr std::array<Part, LRTs> lr_types{{L, R}};
constexpr std::array<Part, LRMTs> lrm_types{{L, R, MS}};
constexpr std::array<Part, LRMTs> lmr_types{{L, MS, R}};

// Named panel hist type indices
struct PanelCountType {
  const char * name() const {
    switch (index) {
      case 0:
        return "Mapping Hits";
      case 1:
        return "Assigned Pair";
      case 2:
        return "Pair Matching Assignment";
      default:
        throw Error("Bad PanelCountType value") << index;
    }
  }
  const char * color() const {
    switch (index) {
      case 0:
        return "1 0 0";
      case 1:
        return "0 0 1";
      case 2:
        return "0 0 0";
      default:
        throw Error("Bad PanelCountType value") << index;
    }
  }
  unsigned int index;
  constexpr operator unsigned int() const { return index; }
};
static constexpr PanelCountType PH{0};
static constexpr PanelCountType PA{1};
static constexpr PanelCountType PM{2};
static constexpr PanelCountType PCTs{3};
static constexpr std::array<PanelCountType, PCTs> pc_types{{PH, PA, PM}};

// Named assigned / not assigned
struct Assigned {
  uint64_t index;
  const char * name() const {
    switch (index) {
      case 0:
        return "Not Assigned";
      case 1:
        return "Assigned";
      default:
        throw Error("Bad Assigned value") << index;
    }
  }
  bool filled() const { return index; }
  constexpr operator uint64_t() const { return index; }
};
constexpr Assigned NASS{0};
constexpr Assigned ASS{1};
constexpr Assigned ATs{2};
constexpr std::array<Assigned, ATs> na_types{{NASS, ASS}};

// Named C or AC microsatellite
struct MSType {
  std::string sname() const { return name(); }
  const char * name() const {
    switch (index) {
      case 0:
        return "C";
      case 1:
        return "AC";
      default:
        throw Error("Bad MSType value") << index;
    }
  }
  uint64_t index;
  constexpr operator uint64_t() const { return index; }
};
constexpr MSType MSC{0};
constexpr MSType MSAC{1};
constexpr MSType MSTs{2};
constexpr std::array<MSType, MSTs> ms_types{{MSC, MSAC}};

// Named Disrupted or Not Disrupted
struct DisType {
  const char * name() const {
    switch (index) {
      case 0:
        return "Not Disrupted";
      case 1:
        return "Disrupted";
      default:
        throw Error("Bad DisType value") << index;
    }
  }
  const char * color() const {
    switch (index) {
      case 0:
        return "1 0 0";
      case 1:
        return "0 0.8 0";
      default:
        throw Error("Bad DisType value") << index;
    }
  }
  uint64_t index;
  constexpr operator uint64_t() const { return index; }
};
constexpr DisType NotDis{0};
constexpr DisType Dis{1};
constexpr DisType DTs{2};
constexpr std::array<DisType, DTs> dis_types{{NotDis, Dis}};

// Named reference type
struct RefType {
  const char * name() const {
    switch (index) {
      case 0:
        return "4N";
      case 1:
        return "3N";
      case 2:
        return "C2T";
      default:
        throw Error("Bad RefType value") << index;
    }
  }
  uint64_t index;
  constexpr operator uint64_t() const { return index; }
};
constexpr RefType R4N{0};
constexpr RefType R3N{1};
constexpr RefType RC2T{2};
constexpr RefType R43CTs{3};
constexpr std::array<RefType, R43CTs> ref_types{{R4N, R3N, RC2T}};

// Named min / max
struct MinMax {
  const char * name() const {
    switch (index) {
      case 0:
        return "Min";
      case 1:
        return "Max";
      default:
        throw Error("Bad MinMax value") << index;
    }
  }
  uint64_t index;
  constexpr operator uint64_t() const { return index; }
};
constexpr MinMax Min{0};
constexpr MinMax Max{1};
constexpr MinMax MMTs{2};
constexpr std::array<MinMax, MMTs> mm_types{{Min, Max}};

// Named worse / best
struct WorstBest {
  const std::string sname() const { return name(); }
  const char * name() const {
    switch (index) {
      case 0:
        return "Worst";
      case 1:
        return "Best";
      default:
        throw Error("Bad WorstBest value") << index;
    }
  }
  const char * color() const {
    switch (index) {
      case 0:
        return "1 0 0";
      case 1:
        return "0 0.8 0";
      default:
        throw Error("Bad BadGood value") << index;
    }
  }
  constexpr operator uint64_t() const { return index; }
  uint64_t index;
};
constexpr WorstBest Worst{0};
constexpr WorstBest Best{1};
constexpr WorstBest WBTs{2};
constexpr std::array<WorstBest, WBTs> wb_types{{Worst, Best}};
constexpr std::array<WorstBest, WBTs> bw_types{{Best, Worst}};

// Named worse / good
struct BadGood {
  const std::string sname() const { return name(); }
  const char * name() const {
    switch (index) {
      case 0:
        return "Bad";
      case 1:
        return "Good";
      default:
        throw Error("Bad BadGood value") << index;
    }
  }
  const char * color() const {
    switch (index) {
      case 0:
        return "1 0 0";
      case 1:
        return "0 0.8 0";
      default:
        throw Error("Bad BadGood value") << index;
    }
  }
  constexpr operator uint64_t() const { return index; }
  uint64_t index;
};
constexpr BadGood Bad{0};
constexpr BadGood Good{1};
constexpr BadGood BGTs{2};
constexpr std::array<BadGood, BGTs> bg_types{{Bad, Good}};
constexpr std::array<BadGood, BGTs> gb_types{{Good, Bad}};

}  // namespace paa

#endif  // PAA_NAMED_INTS_H_
