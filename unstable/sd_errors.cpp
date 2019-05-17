//
// sd_errors.cpp
//
// investigate sensitive detection error rates
//
// Copyright 2018 Peter Andrews @ CSHL
//

#include <algorithm>
#include <array>
#include <exception>
#include <fstream>
#include <future>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "bridges.h"
#include "encode.h"
#include "error.h"
#include "longSA.h"
#include "mumdex.h"
#include "threads.h"
#include "utility.h"

using std::array;
using std::cout;
using std::cerr;
using std::cref;
using std::endl;
using std::exception;
using std::future;
using std::min;
using std::ofstream;
using std::ostream;
using std::ostringstream;
using std::pair;
using std::ref;
using std::setprecision;
using std::string;
using std::to_string;
using std::vector;

using paa::commas;
using paa::longSA;
using paa::read_optional_formats;
using paa::BridgeInfo;
using paa::Error;
using paa::MUM;
using paa::MUMdex;
using paa::OneBridgeInfo;
using paa::OptionalSavers;
using paa::Pair;
using paa::Reference;
using paa::SimpleHit;
using paa::ThreadPool;

// Matrix with two compile-time dimensions
template <class Val, uint64_t II = 0, uint64_t JJ = 0>
class Matrix {
 public:
  static constexpr uint64_t I{II};
  static constexpr uint64_t J{JJ};
  explicit Matrix(const bool initialize = true) {
    if (initialize) {
      for (uint64_t v{0}; v != I * J; ++v) data[v] = 0;
    }
  }
  Matrix(const std::initializer_list<std::initializer_list<Val> > & input) {
    auto rows = input.begin();
    for (uint64_t i{0}; i != I; ++i, ++rows) {
      auto val = rows->begin();
      for (uint64_t j{0}; j != J; ++j, ++val) (*this)[i][j] = *val;
    }
  }
  Matrix & operator+=(const Matrix & rhs) {
    for (uint64_t v{0}; v != I * J; ++v) data[v] += rhs.data[v];
    return *this;
  }
  uint64_t n_rows() const { return I; }
  uint64_t n_cols() const { return J; }
  const Val * operator[](const uint64_t i) const { return data + i * J; }
  Val * operator[](const uint64_t i) { return data + i * J; }

 private:
  Val data[I * J];
};

constexpr int64_t max_offset{100};
constexpr int64_t n_offsets{2 * max_offset + 1};
constexpr int64_t max_invariant{500};
constexpr int64_t n_invariants{2 * max_invariant + 1};
using IOMatrix = Matrix<uint64_t, n_invariants, n_offsets>;
using UIOMatrix = vector<IOMatrix>;

UIOMatrix & operator+=(UIOMatrix & lhs, const UIOMatrix & rhs) {
  for (uint64_t u{0}; u != lhs.size(); ++u)
    lhs[u] += rhs[u];
  return lhs;
}

struct FullBridgeInfo {
  FullBridgeInfo(const OneBridgeInfo & bridge_,
                 const unsigned int count_,
                 const unsigned int read_offset_,
                 const bool read_2_) :
      bridge{bridge_},
    count{count_},
    read_offset{read_offset_},
    read_2{read_2_} { }
  OneBridgeInfo bridge;  // count of 1 for this sequence
  unsigned int count;  // tag count for this sequence
  unsigned int read_offset;
  bool read_2;
  bool operator<(const FullBridgeInfo & rhs) const {
    if (read_2 == rhs.read_2) {
      if (read_offset == rhs.read_offset) {
        return bridge < rhs.bridge;
      } else {
        return read_offset < rhs.read_offset;
      }
    } else {
      return read_2 < rhs.read_2;
    }
  }
};
bool operator==(const FullBridgeInfo & lhs, const FullBridgeInfo rhs) {
  return !(lhs.bridge < rhs.bridge || rhs.bridge < lhs.bridge);
}
using SimpleBridges = vector<OneBridgeInfo>;
using FullBridges = vector<FullBridgeInfo>;

struct CombinedBridgeInfo {
  explicit CombinedBridgeInfo(const FullBridgeInfo & bridge_) :
      bridge{bridge_.bridge},
    count{bridge_.count} {
  }
  void combine(const FullBridgeInfo & bridge_) {
    bridge.combine(bridge_.bridge);
    count += bridge_.count;
  }
  bool operator<(const FullBridgeInfo & rhs) const {
    return bridge < rhs.bridge;
  }
  bool operator<(const CombinedBridgeInfo & rhs) const {
    if (bridge.bridge_count() == rhs.bridge.bridge_count()) {
      return bridge < rhs.bridge;
    } else {
      return bridge.bridge_count() > rhs.bridge.bridge_count();
    }
  }
  BridgeInfo bridge;  // bridge_count() == n_sequences seen in
  unsigned int count;  // total tag count
};
using CombinedBridges = vector<CombinedBridgeInfo>;

struct SeqInfo {
  SeqInfo(const MUMdex & mumdex, const uint64_t pair_index,
          const string & tag_, const unsigned int count_) :
      tag{tag_},
    count{count_},
    sequences(mumdex.sequences(pair_index)) {
      static thread_local SimpleBridges simple_bridges;
      simple_bridges.clear();
      pair_bridges(mumdex, pair_index, simple_bridges);
      bridges.reserve(simple_bridges.size());
      for (const OneBridgeInfo & bridge : simple_bridges) {
        bridges.emplace_back(bridge, count, 0, false);
      }
    }
  string tag;
  unsigned int count;  // tag count for mate-pair-tag combo
  array<string, 2> sequences;
  FullBridges bridges{};  // All bridges in read 1/2
};
using SeqInfos = vector<SeqInfo>;

using TwoCounts = pair<unsigned int, unsigned int>;
struct TagInfo {
  explicit TagInfo(const SeqInfo & seq_info) :
      tag{seq_info.tag},
    count{seq_info.count},
    counts{{seq_info.count, seq_info.bridges.size()}},
    bridges(seq_info.bridges) { }
  void add(const SeqInfo & seq_info) {
    if (tag != seq_info.tag) throw Error("tag mismatch");
    count += seq_info.count;
    counts.emplace_back(seq_info.count, seq_info.bridges.size());
    bridges.insert(bridges.begin(),
                   seq_info.bridges.begin(), seq_info.bridges.end());
  }
  void finalize() {
    sort(counts.begin(), counts.end(),
         [](const TwoCounts lhs, const TwoCounts rhs) {
           return rhs < lhs;
         });
    sort(bridges.begin(), bridges.end());
    for (const FullBridgeInfo & bridge : bridges) {
      if (combined_bridges.empty() || combined_bridges.back() < bridge) {
        combined_bridges.emplace_back(bridge);
      } else {
        combined_bridges.back().combine(bridge);
      }
    }
    sort(combined_bridges.begin(), combined_bridges.end());
  }
  unsigned int out(ostream & stream, const uint64_t seq_id,
                   const int64_t invariant, const int64_t offset,
                   const vector<unsigned int> & pos_counts,
                   vector<unsigned int> & loci_hit) const {
    unsigned int result{0};
    for (unsigned int b{0}; b != combined_bridges.size(); ++b) {
      const CombinedBridgeInfo & cbridge{combined_bridges[b]};
      const BridgeInfo & bridge{cbridge.bridge};
      if (bridge.invariant() == invariant &&
          bridge.offset() == offset) {
        ++loci_hit[bridge.pos1()];
        result += bridge.bridge_count();
        stream << seq_id
               << " " << bridge.pos1()
               << " " << invariant
               << " " << offset
               << " " << pos_counts[bridge.pos1()]
               << " " << tag
               << " " << count
               << " " << counts.size()
               << " " << combined_bridges.size()
               << " " << bridge.bridge_count()
               << " " << cbridge.count
               << "\n";
      }
    }
    return result;
  }
  ostream & out(ostream & stream) const {
    // overall tag info
    stream << tag << " " << count << " " << counts.size() << " ";

    // tag counts
    for (unsigned int c{0}; c != counts.size(); ++c) {
      if (c) stream << ',';
      stream << counts[c].first;
    }

    // events per tag
    stream << ' ';
    for (unsigned int c{0}; c != counts.size(); ++c) {
      if (c) stream << ',';
      stream << counts[c].second;
    }

    // combined bridge invariants
    stream << ' ' << combined_bridges.size() << ' ';
    for (unsigned int b{0}; b != combined_bridges.size(); ++b) {
      if (b) stream << ',';
      stream << combined_bridges[b].bridge.invariant();
    }

    // event present counts
    stream << ' ';
    for (unsigned int b{0}; b != combined_bridges.size(); ++b) {
      if (b) stream << ',';
      stream << combined_bridges[b].bridge.bridge_count();
      if (0) stream << setprecision(3)
                    << 1.0 * combined_bridges[b].bridge.bridge_count() / count;
    }

    // event total counts
    stream << ' ';
    for (unsigned int b{0}; b != combined_bridges.size(); ++b) {
      if (b) stream << ',';
      stream << combined_bridges[b].count;
    }
    stream << '\n';
    return stream;
  }

  string tag;
  unsigned int count;
  vector<TwoCounts> counts;
  FullBridges bridges{};
  CombinedBridges combined_bridges{};
};
using TagInfos = vector<TagInfo>;

int main(int argc, char * argv[]) try {
  paa::exit_on_pipe_close();

  if (--argc < 1) throw Error("usage: sd_errors mumdex_name ...");

  const unsigned int min_mum_length{10};

  cerr << "Progress: loading MUMdex files in parallel" << endl;
  ThreadPool pool{static_cast<unsigned int>(min(argc, 32))};
  using TagInfosResult = pair<uint64_t, TagInfos>;
  using TagInfosFuture = future<TagInfosResult>;
  vector<TagInfosFuture> tag_infos_futures;
  const vector<int64_t> invariants_of_interest{
    -5, -4, -3, -2, -1};
  auto process_mumdex = [invariants_of_interest]
      (const string & mumdex_name, vector<unsigned int> & pos_counts,
       vector<unsigned int> & n_possible_indels, uint64_t & n_pairs,
       uint64_t & average_read_length, uint64_t & ref_length) {
    const MUMdex mumdex{mumdex_name};
    const Reference & ref{mumdex.reference()};
    ref_length = ref.size(0);
    const vector<string> optional_formats{read_optional_formats(mumdex_name)};
    OptionalSavers saver{optional_formats};
    saver.load(mumdex_name, mumdex.n_pairs() * 2);
    const uint64_t tag_id{saver.id("XT")};
    const uint64_t count_id{saver.id("XC")};
    // read in info for each read pair with tag combo
    SeqInfos seq_infos;
    seq_infos.reserve(mumdex.n_pairs());
    uint64_t n_bases{0};
    pos_counts.resize(ref.size(0));
    n_pairs = mumdex.n_pairs();
    for (uint64_t p = 0; p != mumdex.n_pairs(); ++p) {
      const string tag{saver[tag_id][p * 2]};
      const unsigned int count{saver[count_id].to_u32(p*2)};
      const Pair & pair{mumdex.pair(p)};
      for (bool read2 : {0, 1}) {
        n_bases += count * (pair.length(read2) - 2 * min_mum_length);
        average_read_length += pair.length(read2);
      }
      seq_infos.emplace_back(mumdex, p, tag, count);
      if (true) {
        for (uint64_t m{mumdex.mums_start(p)}; m != mumdex.mums_stop(p); ++m) {
          const MUM mum{mumdex.mum(m)};
          if (mum.length() < 20) continue;
          for (unsigned int b{0}; b != mum.length(); ++b)
            pos_counts[mum.position0() + b] += count;
        }
      }
    }
    average_read_length /= 2.0 * mumdex.n_pairs();
    if (false && average_read_length > ref.size(0))
      throw Error("unexpected average read length") << average_read_length;

    // sort pairs by tag
    sort(seq_infos.begin(), seq_infos.end(),
         [](const SeqInfo & lhs, const SeqInfo & rhs) {
           if (lhs.tag == rhs.tag) {
             return lhs.count > rhs.count;
           } else {
             return lhs.tag < rhs.tag;
           }
         });

    // group by tag
    TagInfos tag_infos;
    for (const SeqInfo & seq_info : seq_infos) {
      if (tag_infos.empty() || tag_infos.back().tag != seq_info.tag) {
        tag_infos.emplace_back(seq_info);
      } else {
        tag_infos.back().add(seq_info);
      }
    }

    // sort grouped tags by total count
    sort(tag_infos.begin(), tag_infos.end(),
         [](const TagInfo & lhs, const TagInfo & rhs) {
           return lhs.count > rhs.count;
         });

    // combine duplicate bridges for each tag
    for (uint64_t t{0}; t != tag_infos.size(); ++t)
      tag_infos[t].finalize();

    // determine where it is possible to have particular indels
    const longSA sa{ref.fasta_file(), true, true, false};
    n_possible_indels.back() = ref.size(0);
    for (uint64_t i{0}; i != invariants_of_interest.size(); ++i) {
      const int64_t invariant{invariants_of_interest[i]};
      const int64_t deletion_length{-invariant};
      if (invariant >= 0) throw Error("expect negative invariant");
      for (unsigned int b{0}; b != ref.size(0); ++b) {
        unsigned int start_base{0};
        unsigned int stop_base{0};
        const unsigned int sequence_length{static_cast<unsigned int>(
            average_read_length + deletion_length)};
        if (b < ref.size(0) / 2) {
          start_base = b > sequence_length / 2 ?
              b - sequence_length / 2 : 0;
          stop_base = start_base + sequence_length;
        } else {
          stop_base = b + sequence_length / 2 <= ref.size(0) ?
              b + sequence_length / 2 : ref.size(0);
          start_base = stop_base - sequence_length;
        }
        if (b < min_mum_length ||
            b + min_mum_length + deletion_length > ref.size(0)) continue;
        string sequence;
        for (unsigned int s{0}; s != ref.size(0); ++s) {
          if (s >= start_base && s < stop_base &&
              !(s >= b && s < b + deletion_length))
            sequence += ref[0][s];
        }
        if (sequence.size() != average_read_length)
          throw Error("Sequence construction problem");
        const uint64_t high_anchor_pos{b - 1};
        const uint64_t low_anchor_pos{static_cast<uint64_t>(
            b + deletion_length)};
        const std::vector<SimpleHit> mams{
          sa.find_mams(sequence, min_mum_length)};
        if (false) {
          cerr << "map of " << sequence << " found " << mams.size() << " mums"
               << endl;
          cerr << "anchors " << high_anchor_pos << " " << low_anchor_pos
               << endl;
          for (const SimpleHit & mam : mams) {
            cerr << mam.dir << " " << mam.chr << " " << mam.pos
                 << " " << mam.len << " " << mam.off << endl;
          }
        }
        bool has_high{false};
        bool has_low{false};
        for (const SimpleHit & mam : mams) {
          if (mam.pos + mam.len - 1 == high_anchor_pos) has_high = true;
          if (mam.pos == low_anchor_pos) has_low = true;
        }
        if (has_high && has_low) ++n_possible_indels[i];
      }
    }

    return TagInfosResult{n_bases, move(tag_infos)};
  };

  // process mumdex inputs in parallel
  vector<vector<unsigned int>> all_pos_counts(argc);
  vector<vector<unsigned int>> all_n_possible_indels(
      argc, vector<unsigned int>(invariants_of_interest.size() + 1));
  vector<uint64_t> all_n_pairs(argc);
  vector<uint64_t> all_read_lengths(argc);
  vector<uint64_t> all_ref_lengths(argc);
  for (int arg{0}; arg != argc; ++arg)
    tag_infos_futures.push_back(pool.run(process_mumdex, argv[arg + 1],
                                         ref(all_pos_counts[arg]),
                                         ref(all_n_possible_indels[arg]),
                                         ref(all_n_pairs[arg]),
                                         ref(all_read_lengths[arg]),
                                         ref(all_ref_lengths[arg])));

  // get future results, group by tags
  uint64_t n_bases{0};
  vector<TagInfos> all_tag_infos;
  for (TagInfosFuture & future : tag_infos_futures) {
    const TagInfosResult result{future.get()};
    n_bases += result.first;
    all_tag_infos.push_back(move(result.second));
  }
  cerr << "events detectable in a total of " << commas(n_bases) << " bases\n\n";

  uint64_t average_read_length{0};
  uint64_t total_reference_length{0};
  for (int a{0}; a != argc; ++a) {
    average_read_length += all_read_lengths[a];
    total_reference_length += all_ref_lengths[a];
  }
  average_read_length /= argc;

  vector<unsigned int> total_possible_indels(invariants_of_interest.size());
  for (const vector<unsigned int> & n_possible_indels : all_n_possible_indels) {
    for (unsigned int i{0}; i != n_possible_indels.size(); ++i) {
      if (i + 1 != n_possible_indels.size())
        total_possible_indels[i] += n_possible_indels[i];
      cerr << (i ? " " : "") << n_possible_indels[i];
    }
    cerr << endl;
  }

  // output detailed tag info
  if (true) {
    cerr << "Progress: output tag_info files" << endl;
    for (uint64_t i{0}; i != all_tag_infos.size(); ++i) {
      ofstream tag_file{"tag_info." + to_string(i) + ".txt"};
      tag_file << "tag total_count num_seqs tag_counts event_counts n_bridges "
               << "invariants present_counts total_counts" << endl;
      for (uint64_t t{0}; t != all_tag_infos[i].size(); ++t)
        all_tag_infos[i][t].out(tag_file);
    }
  }

  if (false) {
    cout << "all_pos_counts " << all_pos_counts.size() << endl;
    for (unsigned int s{0}; s != all_pos_counts.size(); ++s) {
      const vector<unsigned int> & pc{all_pos_counts[s]};
      cout << s << " " << pc.size() << endl;
      for (unsigned int b{0}; b != pc.size(); ++b) {
        cout << s << " " << b << " " << pc[b] << endl;
      }
    }
  }

  // Output info for -N 1 events
  if (true) {
    uint64_t n_pairs{0};
    for (const uint64_t pairs : all_n_pairs) n_pairs += pairs;
    cerr << "Progress: output event_info files" << endl;
    for (unsigned int inv{0}; inv != invariants_of_interest.size(); ++inv) {
      const int64_t invariant{invariants_of_interest[inv]};
      const int64_t offset_of_interest{invariant ?
            (invariant > 0 ? invariant + 1 : 1) : 2};
      vector<unsigned int> loci_hit(500);
      for (uint64_t i{0}; i != all_tag_infos.size(); ++i) {
        ofstream tag_file{"event_info." + to_string(i) + "." +
              to_string(invariant) + "." + to_string(offset_of_interest) +
              ".txt"};
        tag_file << "seq pos inv off coverage tag tag_count tag_seqs "
                 << "tag_events event_seqs event_tag_count" << endl;
        for (uint64_t t{0}; t != all_tag_infos[i].size(); ++t)
          all_tag_infos[i][t].out(tag_file, i, invariant, offset_of_interest,
                                  all_pos_counts[i], loci_hit);
      }
      unsigned int n_hits{0};
      unsigned int n_loci{0};
      for (const unsigned int hits : loci_hit) {
        if (hits) {
          ++n_loci;
          n_hits += hits;
        }
      }

      const uint64_t average_coverage{
        n_pairs * 2 * average_read_length / total_reference_length};
      const double ppse{1.0 * average_coverage * total_possible_indels[inv] /
            (n_hits ? n_hits : 1)};

      cerr << "inv " << invariant
           << " off " << offset_of_interest
           << " hits " << commas(n_hits)
           << " hit_loci " << n_loci
           << " possible_loci " << commas(total_possible_indels[inv])
           << " coverage " << commas(average_coverage)
           << " pairs_per_event "
           << commas(ppse)
           << endl;
    }
  }

  // count errors by type and offset
  cerr << "Progress: counting errors by invariant and offset" << endl;
  auto count_tags = [](const TagInfos & tag_infos) {
    UIOMatrix counts(2);
    for (uint64_t t{0}; t != tag_infos.size(); ++t) {
      const TagInfo & tag_info{tag_infos[t]};
      const CombinedBridges & bridges{tag_info.combined_bridges};
      for (unsigned int b{0}; b != bridges.size(); ++b) {
        const BridgeInfo & bridge{bridges[b].bridge};
        int invariant{static_cast<int>(bridge.invariant())};
        int offset{bridge.offset()};
        if (invariant < -max_invariant) invariant = -max_invariant;
        if (invariant > max_invariant) invariant = max_invariant;
        if (offset < -max_offset) offset = -max_offset;
        if (offset > max_offset) offset = max_offset;
        const bool unoriented{bridge.high1() == bridge.high2()};
        ++counts[unoriented][invariant + max_invariant][offset + max_offset];
      }
    }
    return counts;
  };

  // Process counts in parallel
  using CountsFuture = future<UIOMatrix>;
  vector<CountsFuture> counts_futures;
  for (uint64_t i{0}; i != all_tag_infos.size(); ++i)
    counts_futures.push_back(pool.run(count_tags, cref(all_tag_infos[i])));
  UIOMatrix counts(2);
  for (CountsFuture & future : counts_futures) counts += future.get();

  // output rates for errors
  auto type = [](
      const bool unoriented, const int invariant, const int offset) {
    if (unoriented) {
      return "inversion";
    } else {
      if (invariant - max_invariant == 0) {
        if (offset - max_offset == 2) {
          return "SNP";
        } else {
          return "substitution";
        }
      } else {
        if (invariant - max_invariant > 0) {
          return "insertion";
        } else {
          return "deletion";
        }
      }
    }
  };
  cerr << "Progress: calculate overall rate information" << endl;
  using Item = pair<uint64_t, string>;
  vector<Item> sorted_counts;
  for (int u{0}; u != static_cast<int>(counts.size()); ++u) {
    for (int i{0}; i != static_cast<int>(counts[u].n_rows()); ++i) {
      uint64_t total{0};
      ostringstream out;
      for (int o{0}; o != static_cast<int>(counts[u].n_cols()); ++o) {
        const uint64_t count{counts[u][i][o]};
        if (!count) continue;
        total += count;
        out << type(u, i, o)
            << " " << i - max_invariant << " " << o - max_offset;
        sorted_counts.emplace_back(count, out.str());
        out.str("");
      }
      if (!total) continue;
      out << type(u, i, 1) << " " << i - max_invariant << " " << "all";
      sorted_counts.emplace_back(total, out.str());
    }
  }
  sort(sorted_counts.begin(), sorted_counts.end(),
       [](const Item lhs, const Item rhs) {
         return rhs < lhs;
       });
  cout << "Type Invariant Offset EventCount 1EventInNbases\n";
  for (uint64_t c{0}; c != sorted_counts.size(); ++c) {
    const Item item{sorted_counts[c]};
    const double rate{1.0 * n_bases / item.first};
    cout << item.second << " " << item.first
         << " " << commas(rate) << "\n";
  }

  return 0;
} catch (Error & e) {
  cerr << "paa::Error:" << endl;
  cerr << e.what() << endl;
  return 1;
} catch (exception & e) {
  cerr << "std::exception" << endl;
  cerr << e.what() << endl;
  return 1;
} catch (...) {
  cerr << "unknown exception was caught" << endl;
  return 1;
}
