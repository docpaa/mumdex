//
// process_bisulfite_rna
//
// makes bisulfite converted references
//
// Copyright 2016 Peter Andrews @ CSHL
//

#include <algorithm>
#include <exception>
#include <fstream>
#include <future>
#include <iostream>
#include <map>
#include <mutex>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "error.h"
#include "files.h"
#include "gencode.h"
#include "mapper.h"
#include "mumdex.h"
#include "strings.h"
#include "utility.h"

using std::async;
using std::cerr;
using std::cout;
using std::condition_variable;
using std::endl;
using std::exception;
using std::future;
using std::ifstream;
using std::istream;
using std::map;
using std::move;
using std::mutex;
using std::ofstream;
using std::ostringstream;
using std::sort;
using std::set;
using std::string;
using std::unique_lock;
using std::vector;

using paa::ChromosomeIndexLookup;
using paa::Error;
using paa::Exon;
using paa::Interval;
using paa::GenCodeInfo;
using paa::Mapper;
using paa::Reference;
using paa::SimpleHit;
using paa::longSA;
using paa::readable;
using paa::replace_substring;
using paa::serr;
using paa::sout;
using paa::tout;

mutex cout_mutex;
const bool map_read_ahead{false};
const std::chrono::seconds sleep_duration{1};

// Convert c to t, keep others the same
inline char convert_c2t(const char base) {
  switch (base) {
    case 'c':
      return 't';
    case 'C':
      return 'T';
    default:
      return base;
  }
}

// Convert g to a, keep others the same
inline char convert_g2a(const char base) {
  switch (base) {
    case 'g':
      return 'a';
    case 'G':
      return 'A';
    default:
      return base;
  }
}

// Make converted references if not already exists
class ConvertedReferenceCreator {
 public:
  ConvertedReferenceCreator(const Reference & ref,
                            const string & g2a_name,
                            const string & c2t_name) {
    if (!readable(c2t_name) || !readable(g2a_name)) {
      ofstream c2t{c2t_name.c_str()};
      if (!c2t) throw Error("Problem opening file") << c2t_name;
      ofstream g2a{g2a_name.c_str()};
      if (!g2a) throw Error("Problem opening file") << g2a_name;

      for (unsigned int c = 0; c != ref.n_chromosomes(); ++c) {
        c2t << '>' << ref.name(c);
        g2a << '>' << ref.name(c);
        for (unsigned int b = 0; b != ref.size(c); ++b) {
          if (b % 50 == 0) {
            c2t << endl;
            g2a << endl;
          }
          const char base{ref[c][b]};
          c2t << convert_c2t(base) << endl;
          g2a << convert_g2a(base) << endl;
        }
      }
      c2t.close();
      g2a.close();
      // Create suffix arrays
      {
        Mapper c2t_mapper{c2t_name};
      }
      Mapper g2a_mapper{g2a_name};
    }
  }
};

void reverse_complement(string & bases) {
  reverse(bases.begin(), bases.end());
  for (auto & base : bases) {
    base = base == 'A' ? 'T' :
        (base == 'T' ? 'A' :
         (base == 'C' ? 'G' :
          (base == 'G' ? 'C' :
           (base == 'N' ? 'N' : 'X'))));
    if (base == 'X') throw Error("strange base seen");
  }
}

// Convert sam to object
class ReadInfo {
 public:
  explicit ReadInfo(istream & sam) {
    sam >> name >> flag >> chr >> pos
        >> mapq >> cigar >> mchr >> mpos >> tlen >> bases >> errors;
    getline(sam, optional);
    if (!sam) return;

    // Die on bad / inappropriate reads
    if ((flag & 0x100) ||  // secondary alignment
        (flag & 0x200) ||  // bad vendor quality
        (flag & 0x800)) {  // supplementary alignment
      throw Error("Bad read flag");
    }

    // Reverse complement if necessary
    if (flag & 0x10) {
      reverse(errors.begin(), errors.end());
      reverse_complement(bases);
    }
    read2 = flag & 128;
    complete = true;
  }

  string name{};
  unsigned int flag{};
  string chr{};
  unsigned int pos{};
  unsigned int mapq{};
  string cigar{};
  string mchr{};
  unsigned int mpos{};
  unsigned int tlen{};
  string bases{};
  string errors{};
  string optional{};
  bool read2{};
  bool complete{false};
};

class PairInfo {
 public:
  explicit PairInfo(istream & sam) : read1{sam}, read2{sam} {
    if (!read1.complete || !read2.complete) return;
    if (read1.name != read2.name)
      throw Error("Sam input not namesorted");
    if (read1.read2 || !read2.read2) {
      throw Error("Bad pair mate numbers")
          << read1.name << read1.flag << read2.flag;
    }
    complete = true;
  }

  ReadInfo read1;
  ReadInfo read2;
  bool complete{false};
};

class ReadMUM : public SimpleHit {
 public:
  ReadMUM(const SimpleHit & hit, const bool read2_,
          const unsigned int read_length) :
      SimpleHit{hit}, read2{read2_} {
        if (read2) {
          if (dir == '+') {
            dir = '-';
          } else {
            dir = '+';
          }
          off = read_length - off - len;
        }
      }
  bool read2;
  bool operator==(const ReadMUM & other) const {
    return static_cast<SimpleHit>(*this) == static_cast<SimpleHit>(other) &&
        read2 == other.read2;
  }
};

namespace paa {
bool operator<(const SimpleHit & lhs, const SimpleHit & rhs) {
  if (lhs.chr == rhs.chr) {
    if (lhs.dir == rhs.dir) {
        return lhs.off < rhs.off;
    } else {
      return lhs.dir < rhs.dir;
    }
  } else {
    return lhs.chr < rhs.chr;
  }
}
}  // namespace paa

bool operator<(const ReadMUM & lhs, const ReadMUM & rhs) {
  if (lhs.chr == rhs.chr) {
    if (lhs.dir == rhs.dir) {
      if (lhs.read2 == rhs.read2) {
        return lhs.off < rhs.off;
      } else {
        return lhs.read2 < rhs.read2;
      }
    } else {
      return lhs.dir < rhs.dir;
    }
  } else {
    return lhs.chr < rhs.chr;
  }
}

const unsigned int default_overlap_allowed{20};
unsigned int overlap_allowed(const SimpleHit &, const SimpleHit &) {
  return default_overlap_allowed;
}

unsigned int overlap_allowed(const ReadMUM & lhs, const ReadMUM & rhs) {
  return lhs.read2 ==rhs.read2 ? default_overlap_allowed : 1000U;
}

void update_positions(set<unsigned int> & positions,
                      const SimpleHit & mum) {
  for (unsigned int pos{mum.pos}; pos != mum.pos + mum.len; ++pos) {
    positions.insert(pos);
  }
}

// Get the best set of reads from mums
template<class Mapping>
unsigned int find_best_mappings(vector<Mapping> & mums,
                                vector<Mapping> & used_mums) {
  // Sort by chr, dir, (read,) offset
  sort(mums.begin(), mums.end());

  // Find best score
  unsigned int best_score{0};
  vector<unsigned int> ms;
  vector<unsigned int> best_ms;
  for (unsigned int m_start{0}; m_start != mums.size(); ++m_start) {
    // const Mapping & mum_start{mums[m_start]};
    ms.push_back(m_start);
    unsigned int first_skipped{0};
    for (unsigned int m_stop{m_start + 1}; m_stop != mums.size(); ++m_stop) {
      const Mapping & mum_last{mums[ms.back()]};
      const Mapping & mum_stop{mums[m_stop]};

      // Is run done?
      const bool stop{[&mum_last, &mum_stop]() {
          if (mum_last.chr != mum_stop.chr) return true;
          if (mum_last.dir != mum_stop.dir) return true;
          return false;
        }()};
      if (stop) {
        break;
      }

      // Is run possibly interrupted
      const bool skip{[&mum_last, &mum_stop]() {
          constexpr unsigned int max_distance{100000};
          const unsigned int allowed_overlap{
            overlap_allowed(mum_last, mum_stop)};
          if (mum_last.dir == '+') {
            if (mum_last.pos + mum_last.len + max_distance < mum_stop.pos)
              return true;
            if (mum_last.pos + mum_last.len > mum_stop.pos + allowed_overlap)
              return true;
          } else {
            if (mum_stop.pos + mum_stop.len + max_distance < mum_last.pos)
              return true;
            if (mum_stop.pos + mum_stop.len > mum_last.pos + allowed_overlap)
              return true;
          }
          return false;
        }()};
      if (skip) {
        if (!first_skipped) {
          first_skipped = m_stop;
        }
        continue;
      }

      // Update mum list
      ms.push_back(m_stop);
    }

    set<unsigned int> positions;
    for (const unsigned int m : ms) {
      update_positions(positions, mums[m]);
    }
    const unsigned int score{static_cast<unsigned int>(positions.size())};

    // See if score is best
    if (score > best_score) {
      best_score = score;
      best_ms = ms;
    }

    // Where to continue
    if (first_skipped) {
      m_start = first_skipped - 1;
    } else {
      m_start = ms.back();
    }

    // Reset for next check
    first_skipped = 0;
    ms.clear();
  }

  // Keep only selected mums
  for (const unsigned int m : best_ms) {
    used_mums.push_back(mums[m]);
  }
  // Diagnostic info
  if (false && best_ms.size() > 1) {
    unique_lock<mutex> cout_lock(cout_mutex);
    sout << "read stats"
         << best_score << mums.size() << used_mums.size() << endl;
    for (unsigned int m{0}; m != mums.size(); ++m) {
      const Mapping & mum{mums[m]};
      sout << (find(best_ms.begin(), best_ms.end(), m) != best_ms.end());
      sout << m << mum.chr << mum.pos << mum.dir << mum.off << mum.len << endl;
    }
  }

  return best_score;
}

// Map pair, get mums and mapping scores
class MappedInfo {
 public:
  MappedInfo(const Mapper * mappers[], const ReadInfo & read) {
    for (const bool c2t_read :
      {(read.read2 ? true : false), (read.read2 ? false : true)}) {
      string transformed_bases;
      for (const char base : read.bases) {
        transformed_bases += c2t_read ?
            convert_c2t(base) : convert_g2a(base);
      }
      for (const bool c2t_map : {false, true}) {
        all_mums.push_back(mappers[c2t_map]->sa().find_mams(transformed_bases));
        vector<bool> read_covered(read.bases.size());
        unsigned int biggest_mum{0};
        for (const auto & mum : all_mums.back()) {
          for (unsigned int b{mum.off}; b != mum.off + mum.len; ++b) {
            read_covered[b] = true;
          }
          if (biggest_mum < mum.len) {
            biggest_mum = mum.len;
          }
        }
        unsigned int n{0};
        for (const bool c : read_covered) {
          if (c) ++n;
        }

        ostringstream out;
        out << all_mums.back().size() << ':' << n << ':'
            << biggest_mum << ':';

        // Find best mum combo
        mums.push_back(vector<SimpleHit>());
        scores.push_back(find_best_mappings(all_mums.back(), mums.back()));

        out << mums.back().size();
        text.push_back(out.str());
      }
    }
  }

  vector<unsigned int> scores{};
  vector<vector<SimpleHit>> mums{};
  vector<vector<SimpleHit>> all_mums{};
  vector<string> text{};
};

// Pick the best mapping based on scores
class MappingPicker {
 public:
  MappingPicker(const string & sam_name,
                const Reference & ref_,
                const string & gencode_name,
                const string & g2a_name,
                const string & c2t_name,
                const unsigned int n_threads_) :
      sam{sam_name.c_str()},
    gencode{gencode_name, ref_},
    ref{ref_},
    g2a_mapper{g2a_name, map_read_ahead},
    c2t_mapper{c2t_name, map_read_ahead},
    n_threads{n_threads_},
    start_time{time(nullptr)} {
      if (!sam) throw Error("Problem opening sam file") << sam_name;

      // Launch threads
      for (unsigned int t{0}; t != n_threads; ++t) {
        futures.push_back(async(std::launch::async, std::ref(*this)));
      }
    }

  MappingPicker(const MappingPicker &) = delete;
  MappingPicker & operator=(const MappingPicker &) = delete;

  ~MappingPicker() try {
    // Wait for threads to finish, then collect
    // Check all occasionally for exceptions
    while (futures.size()) {
      for (vector<future<void>>::iterator f{futures.begin()};
           f != futures.end(); ++f) {
        const std::future_status status{f->wait_for(sleep_duration)};
        if (status == std::future_status::ready) {
          f->get();
          futures.erase(f);
          break;
        }
      }
    }

    // Output stats
    const time_t stop_time{time(nullptr)};
    const uint64_t elapsed(stop_time - start_time);
    cerr << "stats" << endl;
    serr << "N pairs is" << n_pairs << endl;
    serr << "choice made" << choice_made << "of" << n_pairs
         << "is" << 1.0 * choice_made / n_pairs << endl;
    serr << "choice average score is" << 1.0 * choice_average / choice_made
         << endl;
    serr << "choice average score difference is"
         << 1.0 * choice_difference / choice_made << endl;
    serr << "gene found" << gene_found << "of" << n_pairs
         << "is" << 1.0 * gene_found / n_pairs << endl;
    serr << "r1 r2 agree on choice" << r1r2_agree << "of" << n_pairs
         << "is" << 1.0 * r1r2_agree / n_pairs << endl;
    serr << "r1 r2 disagree on choice" << r1r2_disagree << "of" << n_pairs
         << "is" << 1.0 * r1r2_disagree / n_pairs << endl;
    serr << "all agree on choice" << all_agree << "of" << n_pairs
         << "is" << 1.0 * all_agree / n_pairs << endl;
    serr << "chr agree for choice" << chr_agree << "of" << all_agree
         << "is" << 1.0 * chr_agree / all_agree << endl;
    serr << "timing" << elapsed << "seconds or"
         << static_cast<unsigned int>(60.0 * n_pairs / elapsed)
         << "pairs per minute using" << n_threads << "threads or"
         << static_cast<unsigned int>(60.0 * n_pairs / elapsed / n_threads)
         << "per thread" << endl;
  } catch (Error & e) {
    cerr << "paa::Error:" << endl;
    cerr << e.what() << endl;
    return;
  } catch (exception & e) {
    cerr << "std::exception" << endl;
    cerr << e.what() << endl;
    return;
  }

  void operator()() {
    unique_lock<mutex> sam_lock(sam_mutex);
    while (sam) {
      const PairInfo pair{sam};
      if (!pair.complete) {
        if (sam) throw Error("Incomplete pair, but sam looks good");
        break;
      }
      sam_lock.unlock();
      process_pair(pair);
      sam_lock.lock();
    }
  }

  void process_pair(const PairInfo & pair) {
    const MappedInfo info1{mappers, pair.read1};
    const MappedInfo info2{mappers, pair.read2};

    const ReadInfo * reads[2]{&pair.read1, &pair.read2};
    const MappedInfo * info[2]{&info1, &info2};

    // Get best combines scores
    vector<vector<ReadMUM>> used_mums(4);
    vector<vector<ReadMUM>> all_mums(4);
    vector<unsigned int> combined_scores{0, 0, 0, 0};
    for (const unsigned int map : {0, 1, 2, 3}) {
      vector<ReadMUM> read_mums;
      for (const unsigned int read2 : { false, true }) {
        vector<SimpleHit> simple_mums{info[read2]->all_mums[map]};
        for (const SimpleHit & mum : simple_mums) {
          read_mums.emplace_back(mum, read2, reads[read2]->bases.size());
        }
      }
      combined_scores[map] = find_best_mappings(read_mums, used_mums[map]);
      all_mums[map] = read_mums;
    }

    // Pick best mappings by score
    unsigned int best_map[3]{0, 0, 0};  // read 1, read 2, overall
    bool equal_best[2]{false, false};
    for (const unsigned int map : {1, 2, 3}) {
      if (combined_scores[map] >= combined_scores[best_map[2]]) {
        best_map[2] = map;
      }
      for (const unsigned int read2 : { false, true }) {
        if (info[read2]->scores[map] >= info[read2]->scores[best_map[read2]]) {
          best_map[read2] = map;
        }
      }
    }

    // Check for equally good mappings by score
    for (const unsigned int map : {0, 1, 2, 3}) {
      for (const unsigned int read2 : { false, true }) {
        if (map != best_map[read2] && info[read2]->scores[map] ==
            info[read2]->scores[best_map[read2]]) {
          equal_best[read2] = true;
        }
      }
    }

    // Find next best map
    unsigned int next_best_map{best_map[2] ? 0U : 1U};
    for (const unsigned int map : {0, 1, 2, 3}) {
      if (map != best_map[2] &&
          combined_scores[map] > combined_scores[next_best_map]) {
        next_best_map = map;
      }
    }

    // Best score
    const unsigned int best_score{combined_scores[best_map[2]]};
    const unsigned int next_best_score{combined_scores[next_best_map]};

    // Find gene exon overlaps
    map<string, unsigned int> genes_seen;
    if (best_score != next_best_score) {
      const vector<ReadMUM> & mums{used_mums[best_map[2]]};
      for (const ReadMUM & mum : mums) {
        const Interval interval{mum.chr, mum.pos, mum.pos + mum.len};
        for (vector<Exon>::const_iterator check{gencode.begin(interval)};
             check != gencode.exons.end(); ++check) {
          if (check->chromosome > interval.chromosome ||
              check->lowest_later_start >= interval.stop_position) {
            break;
          }
          if (check->start_position < interval.stop_position) {
            // Found exon interval overlap
            ++genes_seen[gencode.entries[check->gene].gene_name];
          }
        }
      }
    }

    // Update stats
    unique_lock<mutex> stats_lock(stats_mutex);
    if (genes_seen.size()) ++gene_found;

    // See if reads agree on mapping
    // bool agree{false};
    if (!equal_best[0] && !equal_best[1]) {
      if (best_map[0] == best_map[1]) {
        ++r1r2_agree;
        // agree = true;
      } else {
        ++r1r2_disagree;
      }
    }

    // Make a choice
    if (best_score && next_best_score != best_score) {
      ++choice_made;
      choice_average += best_score;
      choice_difference += best_score - next_best_score;

      if (!equal_best[0] && !equal_best[1] &&
          best_map[0] == best_map[1] && best_map[2] == best_map[0]) {
        ++all_agree;
        chr_agree += info[0]->mums[best_map[0]].front().chr ==
            info[1]->mums[best_map[1]].front().chr;
      }
    }
    stats_lock.unlock();

#if 1
    ostringstream sam_out;


    // const bool read1_flipped{used_mums[best_map[2]].front().dir == '-'};
    // const bool read_flipped[2]{read1_flipped, !read1_flipped};

    // used_mums[best_map[2]].front().pos;

    const bool chosen{best_score != next_best_score};
    if (chosen) {
      if (used_mums[best_map[2]].front().chr !=
          used_mums[best_map[2]].back().chr)
        throw Error("Chromosome mismatch");
      // chr = used_mums[best_map[2]].front().chr;
    }

#if 0
    const bool read_mapped[2]{
      chosen ? (used_mums[best_map[2]].front().read2 == false) : false,
          chosen ? (used_mums[best_map[2]].back().read2 == true) : false};
#endif

    for (const unsigned int read2 : { false, true }) {
      unsigned int flag{0};
      string chr{chosen ?
            ref.name(used_mums[best_map[2]].front().chr) : "*"};
      unsigned int pos{0};  // used_mums[best_map[2]].front().pos};
      unsigned int mapq{combined_scores[best_map[2]]};
      unsigned int cigar{0};
      unsigned int mpos{0};
      unsigned int tlen{0};
      string extra{"moo"};

      string bases{reads[read2]->bases};
      string errors{reads[read2]->errors};
      if (true) {
        reverse_complement(bases);
        reverse(errors.begin(), errors.end());
      }

      sam_out << reads[read2]->name << '\t'
              << flag << '\t'
              << chr << '\t'
              << pos << '\t'
              << mapq << '\t'
              << cigar << '\t'
              << chr << '\t'
              << mpos << '\t'
              << tlen << '\t'
              << bases << '\t'
              << errors << '\t'
              << extra
              << reads[read2]->optional << endl;
    }
    if (1 && chosen) {
      sam_out << "mapping info "
              << combined_scores[best_map[2]] << " "
              << used_mums[best_map[2]].size()
              << endl;
      for (const ReadMUM & mum : all_mums[best_map[2]]) {
        sam_out << (find(used_mums[best_map[2]].begin(),
                         used_mums[best_map[2]].end(), mum) !=
                    used_mums[best_map[2]].end()) << " ";
        sam_out << mum.read2 << " "
                << ref.name(mum.chr) << " "
                << mum.pos << " "
                << mum.dir << " "
                << mum.off << " "
                << mum.len << endl;
      }
    }
    unique_lock<mutex> cout_lock(cout_mutex);
    ++n_pairs;
    cout << sam_out.str();

#else
    // Prepare for output
    ostringstream out;
    if (best_score == next_best_score) {
      out << ". .";
    } else {
      switch (best_map[2]) {
        case 0:
          out << "read1_g2a map_g2a";
          break;
        case 1:
          out << "read1_g2a map_c2t";
          break;
        case 2:
          out << "read1_c2t map_g2a";
          break;
        case 3:
          out << "read1_c2t map_c2t";
          break;
        default:
          throw Error("No best map");
      }
    }

    const char dir{best_score != next_best_score ?
          used_mums[best_map[2]].front().dir : '.'};
    out << " " << dir;
    if (best_score == next_best_score) {
      out << " .";
    } else {
      out << " " << best_map[2];
    }
    out << " " << best_score << " " << next_best_score;
    out << " " << agree;

    for (const unsigned int read2 : { false, true }) {
      out << " " << reads[read2]->bases.size();
    }

    for (const unsigned int read2 : { false, true }) {
      for (const unsigned int score : info[read2]->scores) {
        out << " " << score;
      }
    }

    for (const unsigned int read2 : { false, true }) {
      for (const string text : info[read2]->text) {
        out << " " << text;
      }
    }

    using keyval = std::pair<string, unsigned int>;
    vector<keyval> ordered_genes(genes_seen.begin(), genes_seen.end());
    sort(ordered_genes.begin(), ordered_genes.end(),
         [](const keyval & lhs, const keyval & rhs) {
           return lhs.second > rhs.second;
         });
    for (const auto & name_val : ordered_genes) {
      out << " " << name_val.first << ":" << name_val.second;
    }
    unique_lock<mutex> cout_lock(cout_mutex);
    cout << ++n_pairs << " " << out.str() << endl;
#endif
  }

  // Input
  mutex sam_mutex{};
  ifstream sam{};

  // Gene info
  const GenCodeInfo gencode;

  // Mapper
  const Reference & ref;
  Mapper g2a_mapper;
  Mapper c2t_mapper;
  const Mapper * mappers[2]{&g2a_mapper, &c2t_mapper};

  unsigned int n_threads{};

  // Threads
  vector<future<void>> futures{};

  // Statistics
  time_t start_time{};
  mutex stats_mutex{};
  uint64_t n_pairs{0};
  uint64_t r1r2_agree{0};
  uint64_t r1r2_disagree{0};
  uint64_t choice_made{0};
  uint64_t gene_found{0};
  uint64_t all_agree{0};
  uint64_t choice_average{0};
  uint64_t choice_difference{0};
  uint64_t chr_agree{0};
};

string convert_name(const string & fasta_name, const string type) {
  string new_name{fasta_name + "." + type + ".fa"};
  replace_substring(new_name, string(".fa.") + type, string(".") + type);
  replace_substring(new_name, string(".fasta.") + type, string(".") + type);
  return new_name;
}

int main(int argc, char * argv[]) try {
  paa::exit_on_pipe_close();

  --argc;
  if (argc != 3 && argc != 4)
    throw Error("usage: process_bisulfite_rna ref_fasta "
                "gencode_file sam [n_threads]");

  // Names of regular and converted reference files
  const string fasta_name{argv[1]};
  const string g2a_name{convert_name(fasta_name, "g2a")};
  const string c2t_name{convert_name(fasta_name, "c2t")};

  // Create regular and converted references and suffix arrays (only once)
  if (!readable(fasta_name + ".bin")) {
    const Mapper mapper{fasta_name};
  }
  const Reference ref{fasta_name};
  const ConvertedReferenceCreator reference_creator{ref, g2a_name, c2t_name};

  // Convert sam and process pairs
  const string gencode_name{argv[2]};
  const string sam_name{argv[3]};
  const unsigned int n_threads{argc == 4 ?
        static_cast<unsigned int>(atoi(argv[4])) : 32};
  MappingPicker picker{sam_name, ref, gencode_name,
        g2a_name, c2t_name, n_threads};

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
