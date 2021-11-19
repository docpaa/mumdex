//
// sequencing_reads.h
//
// Template for a header file
//
// copyright 2021 Peter Andrews
//

#ifndef PAA_SEQUENCING_READS_H_
#define PAA_SEQUENCING_READS_H_

#include <algorithm>
#include <atomic>
#include <map>
#include <set>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "error.h"
#include "gzip.h"
#include "histogram.h"
#include "named_ints.h"
#include "psplot.h"

namespace paa {

template <class Stream>
void fastq_record(Stream & out, const std::string & name,
                  const std::string & sequence1, const std::string & quality1,
                  const std::string & sequence2, const std::string & quality2) {
  out << '@' << name << '\n' << sequence1 << "\n+\n" << quality1 << '\n'
      << '@' << name << '\n' << sequence2 << "\n+\n" << quality2 << '\n';
}

inline bool is_mismatch(const char c1, const char c2) {
  switch (c1) {
    case 'N':
      return false;
    case 'Y':
      return c2 != 'C' && c2 != 'T';
    case 'R':
      return c2 != 'G' && c2 != 'A';
    case 'B':
      return c2 == 'A';
    case 'D':
      return c2 == 'C';
    case 'H':
      return c2 == 'G';
    case 'V':
      return c2 == 'T';
    default:
      return c1 != c2;
  }
}

// Count mismatches
inline uint64_t n_mismatches(const std::string & seq1,
                             const std::string & seq2) {
  uint64_t mismatches{0};
  const uint64_t size{std::min(seq1.size(), seq2.size())};
  for (uint64_t b{0}; b != size; ++b)
    mismatches += is_mismatch(seq1[b], seq2[b]);
  return mismatches;
}

// Holds a raw pair
struct RawPair {
  RawPair() = default;
  RawPair(std::string && pair_name_,
          std::string && read_1, std::string && quality_1,
          std::string && read_2, std::string && quality_2) :
      pair_name{std::move(pair_name_)},
      read1_seq{std::move(read_1)}, read1_qual{std::move(quality_1)},
      read2_seq{std::move(read_2)}, read2_qual{std::move(quality_2)} {
        if (read1_seq.size() != read1_qual.size() ||
            read2_seq.size() != read2_qual.size())
          throw Error("Read sequence - quality length mismatch");
      }
  RawPair(RawPair &&) = default;
  RawPair & operator=(RawPair &&) = default;
  RawPair(const RawPair &) = delete;
  RawPair & operator=(const RawPair &) = delete;
  const std::string * data(const Read & read) const {
    return (read ? &read2_seq : &read1_seq);
  }
  const std::string & data(const Read & read, const SeqQual & sq) const {
    return data(read)[sq];
  }
  std::string pair_name{};
  std::string read1_seq{};
  std::string read1_qual{};
  std::string read2_seq{};
  std::string read2_qual{};
};
using Pairs = std::vector<RawPair>;

class FastqsReader : public RawPair {
 public:
  FastqsReader(const std::string & name1, const std::string & name2,
               const uint64_t n_max_ = 0, const uint64_t n_notify_ = 10000) :
      reads1{name1}, reads2{name2},
      n_max{n_max_}, n_notify{n_notify_} {
        if (!reads1) throw Error("Problem opening fasta file") << name1;
        if (!reads2) throw Error("Problem opening fasta file") << name2;
      }
  bool get_pair() {
    reads1 >> at_char >> pair_.pair_name >> extra;
    if (0) std::cerr << "Read 1 " << at_char << " " << pair_.pair_name
                     << " " << extra << " " << !!reads1 << std::endl;
    if (!reads1 || (n == n_max && n)) return finish();
    reads1 >> pair_.read1_seq >> plus >> pair_.read1_qual;
    if (!reads1) throw Error("Read 1 parse error");
    reads2 >> at_char >> read2_name >> extra;
    if (0) std::cerr << "Read 2 " << at_char << " " << read2_name
                     << " " << extra << " " << !!reads2 << std::endl;
    reads2 >> pair_.read2_seq >> plus >> pair_.read2_qual;
    if (!reads2) throw Error("Read 2 parse error");
    if (pair_.pair_name != read2_name) throw Error("Read name mismatch");
    if (!(++n % n_notify)) notify();
    return true;
  }
  RawPair & pair() { return pair_; }
  uint64_t n_pairs() const { return n; }

 private:
  bool finish() const {
    std::cerr << "\r                             " << std::endl;
    notify(true, std::cout);
    return false;
  }
  void notify(const bool end = false, std::ostream & out = std::cerr) const {
    out << "\rN pairs input: " << n_pairs();
    if (end) out << std::endl;
  }
  Zcat reads1;
  Zcat reads2;
  uint64_t n_max;
  uint64_t n_notify;
  uint64_t n{0};
  RawPair pair_{};
  std::string read2_name{};
  std::string extra{};
  std::string plus{};
  char at_char{};
};

inline std::string seq_or_star(const std::string & seq) {
  return seq.size() ? seq : "*";
}

// Show differences between sequences
template <class Out>
void show_diffs(Out & out, const std::string & pattern, const uint64_t start,
                const std::string & seq) {
  for (uint64_t b{0}; b != seq.size(); ++b) {
    const uint64_t sb{start + b};
    if (sb >= pattern.size()) {
      out << '-';
      continue;
    }
    const char s{seq[b]};
    const char p{pattern[sb]};
    switch (p) {
      case 'N':
        out << ' ';
        break;
      case 'D':
        out << (s == 'C' ? 'D' : ' ');
        break;
      case 'H':
        out << (s == 'G' ? 'H' : ' ');
        break;
      case 'B':
        out << ((s == 'C' || s == 'T') ? ' ' : 'B');
        break;
      case 'F':
        out << ((s == 'G' || s == 'A') ? ' ' : 'F');
        break;
      default:
        out << (p == s ? ' ' : 'X');
        break;
    }
  }
  for (uint64_t b{seq.size()}; b < pattern.size(); ++b) out << '-';
}

// PairStats - standard plots for read pairs
const Marker big_dot{paa::circle(), 0.75};
constexpr uint64_t max_read_length{151};
using ReadCounts = std::vector<uint64_t>;

struct PairStats {
  explicit PairStats(const uint64_t interval_ = 1,
                     const uint64_t max_length_ = max_read_length) :
      interval{interval_}, max_length{max_length_} {}

  // Update stats for each pair seen
  void update(const RawPair & pair) {
    if (++n % interval) return;
    for (const Read & read : r12_types) {
      const std::string & seq{pair.data(read, Seq)};
      const std::string & qual{pair.data(read, Qual)};
      ++read_length[read](qual.size());
      uint64_t quality_sum{0};
      for (uint64_t b{0}; b != seq.size(); ++b) {
        ++bases[read][BaseType::char2base(seq[b])][b];
        const int q{qual[b] - 33};
        ++position_quality[read][0][b];
        position_quality[read][1][b] += q;
        ++base_quality[read](q);
        quality_sum += q;
      }
      if (qual.size()) ++read_quality[read](quality_sum / qual.size());
    }
  }

  // Make final plots from histograms
  void plot(PSDoc & plots, const std::string & method_) const {
    const std::string method{" for " + method_};
    // Pretty base colors by read position
    PSPage & bases_page{*PSPage::create(plots, "Bases Seen" + method, "1 2")};
    for (const Read & read : r12_types) {
      PSGraph & bases_graph{*PSGraph::create(
          bases_page, read.sname() + ";Read Position;Bases Seen",
          Bounds{0, 1.0 * max_length, 0, 1})};
      bases_graph.hist(true);
      std::ostringstream bases_ps[ACGTN];
      for (const BaseType & base : acgtn_types)
        bases_ps[base] << base.color() << " c\n";
      for (uint64_t b{0}; b != max_length; ++b) {
        const uint64_t total{bases[read][0][b] + bases[read][1][b] +
              bases[read][2][b] + bases[read][3][b] + bases[read][4][b]};
        double low{0};
        for (const BaseType & base : acgtn_types) {
          const double high{low + 1.0 * bases[read][base][b] / total};
          bases_ps[base] << " np " << b << " " << low << " gc m "
                         << b + 1 << " " << low << " gc l "
                         << b + 1 << " " << high << " gc l "
                         << b << " " << high << " gc l cp fp\n";
          low = high;
        }
      }
      if (read) {
        std::ostringstream legend;
        legend << lll_font << " sf 0 xfc 0 yfc " << lll_y << " sub m 0 0 0 c "
               << "(Base Colors:) s";
        for (const BaseType & base : acgtn_types) {
          if (base) legend << " 0 0 0 c (,) s";
          legend << " " << base.color() << " c ( " << base.name() << ") s";
        }
        bases_graph.ps(legend.str());
      }
      bases_graph.mid_ps(bases_ps[0].str() + bases_ps[1].str() +
                         bases_ps[2].str() + bases_ps[3].str() +
                         bases_ps[4].str());
    }

    // Quality scores by read position
    PSPage & position_quality_page{
      *PSPage::create(plots, "Base Quality Scores" + method, "1 2")};
    for (const Read & read : r12_types) {
      PSXYSeries & position_quality_series{*PSXYSeries::create(
          position_quality_page,
          read.sname() + ";Read Position;Average Quality", big_dot)};
      for (uint64_t b{0}; b != max_length; ++b)
        if (position_quality[read][0][b])
          position_quality_series.add_point(
              b, 1.0 * position_quality[read][1][b] /
              position_quality[read][0][b]);
    }

    // Histograms of read length, quality scores
    PSPage & read_length_page{
      *PSPage::create(plots, "Read Sizes" + method, "1 2")};
    PSPage & base_quality_page{
      *PSPage::create(plots, "Base Quality Scores" + method, "1 2")};
    PSPage & read_quality_page{
      *PSPage::create(plots, "Read Average Quality Scores" + method, "1 2")};
    for (const Read & read : r12_types) {
      auto log_hist = [&read](PSPage & page, const UIntHistogram hist[RP],
                              const std::string & title) {
        PSHSeries<uint64_t, uint64_t>::create(
            page, hist[read], read.name() + title)->graph().log_y(true);
      };
      log_hist(read_length_page, read_length, ";Read Length;N");
      log_hist(base_quality_page, base_quality, ";Quality Score;N");
      log_hist(read_quality_page, read_quality, ";Average Quality Score;N");
    }
  }

  // Add stats together
  PairStats & operator+=(const PairStats & other) {
    for (const Read & read : r12_types) {
      for (const BaseType & base : acgtn_types)
        for (uint64_t b{0}; b != max_length; ++b)
          bases[read][base][b] += other.bases[read][base][b];
      for (const bool numerator : {false, true})
        for (uint64_t b{0}; b != max_length; ++b)
          position_quality[read][numerator][b] +=
              other.position_quality[read][numerator][b];
      read_length[read] += other.read_length[read];
      base_quality[read] += other.base_quality[read];
      read_quality[read] += other.read_quality[read];
    }
    return *this;
  }

 private:
  ReadCounts empty() const { return ReadCounts(max_length); }
  uint64_t n{0};
  uint64_t interval;
  uint64_t max_length;
  ReadCounts bases[RP][ACGTN]{{empty(), empty(), empty(), empty(), empty()},
    {empty(), empty(), empty(), empty(), empty()}};
  ReadCounts position_quality[RP][2]{{empty(), empty()}, {empty(), empty()}};
  UIntHistogram read_length[RP]{{max_length + 1, 0, max_length + 1},
    {max_length + 1, 0, max_length + 1}};
  UIntHistogram base_quality[RP]{{42, 0, 42}, {42, 0, 42}};
  UIntHistogram read_quality[RP]{{42, 0, 42}, {42, 0, 42}};
};

uint64_t hamming_distance(const std::string & seq1, const std::string & seq2) {
  if (seq1.size() != seq2.size())
    throw Error("Sequence length mismatch in hamming")
        << seq1.size() << seq2.size();
  uint64_t result{0};
  for (uint64_t b{0}; b != seq1.size(); ++b) result += seq1[b] != seq2[b];
  return result;
}

constexpr uint64_t default_tag_length{18};
class SampleBarcodes {
  using BarcodeInfo = std::pair<std::string, uint64_t>;
  using BarcodeIds = std::unordered_map<std::string, uint64_t>;
  using IdNames = std::vector<std::string>;
  using IdLengths = std::vector<uint64_t>;

 public:
  explicit SampleBarcodes(const std::string & barcode_string,
                          const uint64_t max_hamming_ = 1) :
      max_hamming{max_hamming_} {
    using BarcodeNameInfos = std::unordered_map<std::string, BarcodeInfo>;
    const BarcodeNameInfos barcodes{
      [&barcode_string]() {
        BarcodeNameInfos result;
        std::istringstream barcode_stream{barcode_string};
        std::string line;
        uint64_t barcode_length{0};
        while (getline(barcode_stream, line, ',')) {
          std::istringstream line_stream{line};
          std::string barcode;
          std::string name;
          uint64_t length{0};
          getline(line_stream, barcode, ':');
          getline(line_stream, name, ':');
          if (name.empty()) name = "unnamed";
          line_stream >> length;
          if (!length) length = default_tag_length;
          result.emplace(barcode, BarcodeInfo{name, length});
          if (barcode_length && barcode.size() != barcode_length)
            throw Error("Inconsistent barcode lengths");
          barcode_length = barcode.size();
        }
        if (result.empty()) throw Error("No sample barcodes found");
        return result;
      }()};
    using NameIds = std::map<std::string, uint64_t>;
    const NameIds name_ids{[&barcodes]() {
        NameIds result;
        std::vector<std::string> names;
        for (const auto & info : barcodes) names.push_back(info.second.first);
        sort(names.begin(), names.end());
        for (const std::string name : names)
          result.emplace(name, result.size());
        return result;
      }()};
    const std::string first_barcode{barcodes.begin()->first};
    const uint64_t barcode_size{first_barcode.size()};
    pattern_ = barcodes.size() > 1 ?
        std::string(barcode_size, 'N') : first_barcode;
    id_names.resize(barcodes.size());
    id_barcodes.resize(barcodes.size());
    id_template_patterns.resize(barcodes.size());
    for (const auto & info : barcodes) {
      const std::string & name{info.second.first};
      const uint64_t id{name_ids.at(name)};
      barcode_ids.emplace(info.first, id);
      id_names[id] = name;
      id_barcodes[id] = info.first;
      id_template_patterns[id] = std::string(info.second.second, 'D');
    }

    id_names.push_back("unknown");
    id_barcodes.push_back(std::string(pattern_.size(), 'N'));
    id_template_patterns.push_back(std::string(default_tag_length, 'D'));
    std::cout << "Loaded barcodes:";
    for (const auto & info : barcode_ids)
      std::cout << " " << info.first
           << ":" << id_names[info.second]
           << ":" << id_template_patterns[info.second].size();
    std::cout << std::endl;
  }
  uint64_t size() const { return id_names.size(); }
  uint64_t operator[](const std::string & barcode) const {
    if (barcode.size() < pattern_.size()) return barcode_ids.size();
    auto found = barcode_ids.find(barcode);
    if (found == barcode_ids.end()) {
      if (max_hamming) {
        uint64_t best_id{barcode_ids.size()};
        uint64_t best_hamming{1000};
        for (const auto & info : barcode_ids) {
          const uint64_t hamming{hamming_distance(barcode, info.first)};
          if (hamming == best_hamming) {
            best_hamming = 1000;
            best_id = barcode_ids.size();
          } else if (hamming < best_hamming) {
            best_hamming = hamming;
            best_id = info.second;
          }
        }
        if (best_hamming <= max_hamming) return best_id;
      }
      return barcode_ids.size();
    }
    return found->second;
  }
  const std::string & pattern() const { return pattern_; }
  const std::string & name(const uint64_t index) const {
    return id_names[index];
  }
  const std::string & barcode(const uint64_t index) const {
    return id_barcodes[index];
  }
  const std::string & template_pattern(const uint64_t index) const {
    return id_template_patterns[index];
  }

 private:
  const uint64_t max_hamming;
  std::string pattern_{};
  BarcodeIds barcode_ids{};
  IdNames id_names{};
  IdNames id_barcodes{};
  IdNames id_template_patterns{};
};

class SimpleBinning {
 public:
  SimpleBinning(const Reference & reference_, const uint64_t n_bins) :
      reference{reference_}, counts(n_bins),
      bin_size{reference.size() / (n_bins ? n_bins : 1) + 1} {}
  uint64_t bin(const unsigned int chr, const unsigned int pos) const {
    const uint64_t abspos{reference.abspos(chr, pos)};
    return abspos / bin_size;
  }
  void increment(const unsigned int chr, const unsigned int pos) {
    ++counts[bin(chr, pos)];
  }
  uint64_t size() const { return counts.size(); }
  unsigned int operator[](const uint64_t bin) const { return counts[bin]; }
  unsigned int & operator[](const uint64_t bin) { return counts[bin]; }
  double evenness() const {
    std::vector<unsigned int> sorted{counts};
    sort(sorted.begin(), sorted.end());
    const uint64_t q1{sorted[sorted.size() / 4]};
    const uint64_t q3{sorted[3 * sorted.size() / 4]};
    return 1.0 * q1 / q3;
  }
  SimpleBinning & operator+=(const SimpleBinning & rhs) {
    for (uint64_t b{0}; b != counts.size(); ++b) counts[b] += rhs.counts[b];
    return *this;
  }

 private:
  const Reference & reference;
  std::vector<unsigned int> counts;
  uint64_t bin_size;
};

struct SequencePart {
  SequencePart(const std::string & name_, const std::string & sequence_,
               const uint64_t threshold_ = 1, const uint64_t min_ = 0) :
      name{name_}, sequence{sequence_}, threshold{threshold_}, min{min_} {}
  std::string name;
  std::string sequence;
  uint64_t threshold;
  uint64_t min;
  // Note these may vary in reality; here they are fixed size
  void offset(const uint64_t offset__) { offset_ = offset__; }
  uint64_t offset() const { return offset_; }
  uint64_t offset(const BeginEnd & be) const {
    return be ? end_offset() : offset();
  }

 private:
  uint64_t end_offset() const { return offset_ + sequence.size(); }
  uint64_t offset_{0};
};

class ReadStructure {
  using SequenceParts = std::vector<SequencePart>;
  using CI = SequenceParts::const_iterator;

 public:
  ReadStructure(const std::initializer_list<SequencePart> & parts_) :
      parts{parts_} {
    uint64_t offset{0};
    std::ostringstream structure_stream;
    for (uint64_t p{0}; p != parts.size(); ++p) {
      const std::string & name{parts[p].name};
      const std::string & seq{parts[p].sequence};
      structure_stream << (p ? "," : "") << name << ':' << seq.size();
      if (name == "C-MS" || name == "G-MS" || name == "KEY" || name == "MS") {
        key_index_ = p;
        key_sequence_ = seq;
      }
      parts[p].offset(offset);
      offset += seq.size();
      sequence_ += seq;
      if (lookup.emplace(name, p).second == false)
        throw Error("Duplicate sequence part name") << name;
    }
    if (key_index_ == -1ul) key_index_ = parts.size();
    if (offset != max_read_length) {
      static std::mutex mutex{};
      static std::set<std::string> seen;
      std::unique_lock<std::mutex> seen_lock{mutex};
      if (seen.insert(structure_stream.str()).second)
        std::cerr << "Expected read size is " << offset
                  << " for read " << structure_stream.str() << std::endl;
    }
  }
  CI begin() const { return parts.begin(); }
  CI end() const { return parts.end(); }
  uint64_t size() const { return parts.size(); }
  const SequencePart & operator[](const std::string & name) const {
    try {
      return (*this)[lookup.at(name)];
    } catch (...) {
      std::cerr << "Problem looking up sequence part name " << name
                << std::endl;
      throw;
    }
  }
  const SequencePart & operator[](const uint64_t part) const {
    return parts[part];
  }
  const std::string & sequence() const { return sequence_; }
  uint64_t key_index() const { return key_index_; }
  const SequencePart & key() const { return parts[key_index_]; }

 private:
  SequenceParts parts;
  std::unordered_map<std::string, uint64_t> lookup{};
  std::string sequence_{};
  uint64_t key_index_{-1ul};
  std::string key_sequence_{""};
};

struct PairStructure {
  PairStructure(const ReadStructure & read1_, const ReadStructure & read2_) :
      read1{read1_}, read2{read2_} {}
  ReadStructure read1;
  ReadStructure read2;
  const ReadStructure & operator[](const Read & read) const {
    return read ? read2 : read1;
  }
};

class ReadPiece {
  using Piece = std::array<std::string, SQPTs>;

 public:
  explicit ReadPiece(const SequencePart & part_) :
      part{part_}, piece{{"", "", part.sequence}} {}
  std::string & operator[](const SeqQual & sqp) { return piece[sqp]; }
  const std::string & operator[](const SeqQual & sqp) const {
    return piece[sqp];
  }
  uint64_t hamming() const { return hamming_; }
  bool is_bad() {
    const std::string & model{part.sequence};
    hamming_ = n_mismatches(model, piece[Seq]);
    const bool small{piece[Seq].size() < part.min};
    const bool different{hamming() > part.threshold};
    const bool bad{small || different};
    n_bad += bad;
    n_small += small;
    n_different += different;
    n_short += piece[Seq].size() < model.size();
    n_mismatch += hamming();
    n_bases += piece[Seq].size();
    return bad;
  }
  void report(const uint64_t n_pairs) const {
    std::cout << "Piece " << name() << ": "
              << "bad " << 100.0 * n_bad / n_pairs << "%, "
              << "small " << 100.0 * n_small / n_pairs << "%, "
              << "diff " << 100.0 * n_different / n_pairs << "%, "
              << "short " << 100.0 * n_short / n_pairs << "%, "
              << "mis " << 100.0 * n_mismatch / n_bases << "%\n";
  }
  const std::string & name() const { return part.name; }
  ReadPiece & operator+=(const ReadPiece & other) {
    n_bad += other.n_bad;
    n_short += other.n_short;
    n_small += other.n_small;
    n_different += other.n_different;
    n_mismatch += other.n_mismatch;
    n_bases += other.n_bases;
    return *this;
  }

 private:
  const SequencePart & part;
  Piece piece;
  uint64_t n_bad{0};
  uint64_t n_small{0};
  uint64_t n_different{0};
  uint64_t n_short{0};
  uint64_t n_mismatch{0};
  uint64_t n_bases{0};
  uint64_t hamming_{0};
};

class ReadPieces {
  using Lookup = std::unordered_map<std::string, uint64_t>;

 public:
  using Pieces = std::vector<ReadPiece>;
  explicit ReadPieces(const Read & read, const ReadStructure & structure,
                      const bool show_structure = false) :
      key_part{structure.key_index()},
      read_name{"r" + std::to_string(read + 1)},
      pieces(structure.begin(), structure.end()),
      bad_structure{setup_structure(show_structure, read_name)} {
        for (uint64_t p{0}; p != structure.size(); ++p) {
          const auto & part{structure[p]};
          pieces[p][Pattern] = part.sequence;
          lookup.emplace(part.name, lookup.size());
        }
      }
  ReadPiece * find(const std::string & part) {
    auto found = lookup.find(part);
    if (found == lookup.end()) return nullptr;
    return &pieces[found->second];
  }
  uint64_t size() const { return pieces.size(); }
  const ReadPiece & operator[](const std::string & part) const {
    return pieces[lookup.at(part)];
  }
  const ReadPiece & operator[](const uint64_t part) const {
    return pieces[part];
  }
  ReadPiece & operator[](const uint64_t part) { return pieces[part]; }
  void clear() {
    for (ReadPiece & piece : pieces)
      for (const SeqQual & sq : sq_types)
        piece[sq].clear();
  }
  bool get(const std::string data[SQTs]) {
    const std::string & seq{data[Seq]};
    clear();

    // Get first half of read
    uint64_t offset{0};
    for (uint64_t part{0}; part != key_part; ++part) {
      if (offset >= seq.size()) break;
      const uint64_t size{pieces[part][Pattern].size()};
      for (const SeqQual & sq : sq_types)
        pieces[part][sq] = data[sq].substr(offset, size);
      offset += size;
    }

    const bool has_key{key_part != pieces.size()};
    if (has_key) {
      // Search for string in rest of read
      const std::string & search{pieces[key_part + 1][Pattern]};
      uint64_t next_offset{0};
      uint64_t best_length{0};
      for (uint64_t start{offset}; start < seq.size(); ++start) {
        bool mismatch{false};
        uint64_t b{0};
        for (; b != search.size(); ++b) {
          const uint64_t s{start + b};
          if (s == seq.size()) break;
          if (is_mismatch(search[b], seq[s])) {
            if (mismatch) break;
            mismatch = true;
          }
          if (b > best_length) {
            best_length = b;
            next_offset = start;
          }
        }
      }

      // Get second half of read
      if (best_length > 4) {
        for (uint64_t part{key_part}; part != pieces.size(); ++part) {
          if (offset >= seq.size()) break;
          if (part != key_part) next_offset += pieces[part][Pattern].size();
          const uint64_t size{next_offset - offset};
          for (const SeqQual & sq : sq_types)
            pieces[part][sq] = seq.substr(offset, size);
          offset = next_offset;
        }
      } else {
        ++n_parse;
      }
    }

    // Update bad counts
    bool read_bad{false};
    for (uint64_t part{0}; part != pieces.size(); ++part)
      if (pieces[part].is_bad()) read_bad = true;
    n_struct += read_bad;

    // Output bad structure
    if (read_bad && bad_structure) {
      const uint64_t max_ms_size{has_key ?
            std::max(pieces[key_part][Pattern].size(),
                     pieces[key_part][Seq].size()) : 0};
      for (uint64_t p{0}; p != pieces.size(); ++p) {
        const std::string & pattern{pieces[p][Pattern]};
        (*bad_structure) << pattern;
        if (p == key_part)
          (*bad_structure) << std::string(max_ms_size - pattern.size(), '-');
      }
      (*bad_structure) << '\n';
      for (uint64_t p{0}; p != pieces.size(); ++p) {
        const ReadPiece & piece{pieces[p]};
        show_diffs((*bad_structure), piece[Pattern], 0, piece[Seq]);
      }
      (*bad_structure) << '\n';
      for (uint64_t p{0}; p != pieces.size(); ++p) {
        const std::string & sequence{pieces[p][Seq]};
        (*bad_structure) << sequence;
        if (p == key_part)
          (*bad_structure) << std::string(max_ms_size - sequence.size(), '-');
      }
      (*bad_structure) << "\n\n";
    }
    return read_bad;
  }
  void report(const NamedNumber & n_pairs) const {
    if (key_part != pieces.size()) paa::report(n_parse, n_pairs);
    paa::report(n_struct, n_pairs);
    for (const ReadPiece & piece : pieces) piece.report(n_pairs);
  }
  ReadPieces & operator+=(const ReadPieces & other) {
    for (uint64_t piece{0}; piece != size(); ++piece)
      pieces[piece] += other.pieces[piece];
    n_parse += other.n_parse;
    n_struct+= other.n_struct;
    return *this;
  }

 private:
  static std::unique_ptr<Gzip> setup_structure(
      const bool show_structure, const std::string & read_name) {
    if (show_structure) {
      mkdir("structure");
      return std::make_unique<Gzip>(
          "structure/" + read_name + "." + std::to_string(++n) + ".gz");
    } else {
      return nullptr;
    }
  }

  uint64_t key_part;
  std::string read_name;
  Pieces pieces{};
  Lookup lookup{};
  std::unique_ptr<Gzip> bad_structure;
  static std::atomic<uint64_t> n;

 public:
  NamedNumber n_parse{read_name + " bad parse"};
  NamedNumber n_struct{read_name + " bad structure"};
};
std::atomic<uint64_t> ReadPieces::n{0};

class PairPieces {
 public:
  explicit PairPieces(const PairStructure & pair_structure,
                      const bool show_structure = false) :
      pieces{{ReadPieces{R1, pair_structure[R1], show_structure},
          ReadPieces{R2, pair_structure[R2], show_structure}}} {}
  const ReadPieces & operator[](const Read & read) const {
    return pieces[read];
  }
  ReadPieces & operator[](const Read & read) { return pieces[read]; }
  ReadPiece & operator[](const std::string & part) {
    for (const Read & read : r12_types) {
      ReadPiece * piece{pieces[read].find(part)};
      if (piece != nullptr) return *piece;
    }
    throw Error("Could not find piece") << part;
  }
  std::array<bool, RRP> get(const RawPair & pair) {
    std::array<bool, RRP> result;
    result[R1] = pieces[R1].get(pair.data(R1));
    if (do_between) do_between();
    result[R2] = pieces[R2].get(pair.data(R2));
    result[RP] = result[R1] || result[R2];
    n_struct += result[RP];
    return result;
  }
  void report(const NamedNumber & n_pairs) const {
    paa::report(n_struct, n_pairs);
    for (const Read & read : r12_types) pieces[read].report(n_pairs);
  }
  PairPieces & operator+=(const PairPieces & other) {
    n_struct += other.n_struct;
    for (const Read & read : r12_types) pieces[read] += other.pieces[read];
    return *this;
  }
  NamedNumber n_struct{"pairs bad structure"};

  std::function<void(void)> do_between{nullptr};

 private:
  std::array<ReadPieces, RP> pieces;
};

}  // namespace paa

#endif  // PAA_SEQUENCING_READS_H_
