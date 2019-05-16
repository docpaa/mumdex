//
// phase_correlation
//
// Look for phase correlation between fragments across genome
//
// Copyright 2018 Peter Andrews @ CSHL
//

#include <algorithm>
#include <exception>
#include <iostream>
#include <fstream>
#include <functional>
#include <future>
#include <mutex>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "error.h"
#include "files.h"
#include "haha.h"
#include "mumdex.h"
#include "psplot.h"
#include "strings.h"
#include "threads.h"

using std::cout;
using std::cerr;
using std::cref;
using std::endl;
using std::exception;
using std::function;
using std::future;
using std::lock_guard;
using std::mt19937_64;
using std::move;
using std::mutex;
using std::ofstream;
using std::ostringstream;
using std::random_device;
using std::set;
using std::stoul;
using std::string;
using std::to_string;
using std::uniform_real_distribution;
using std::vector;

using paa::bound;
using paa::mkdir;
using paa::replace_substring;
using paa::Bounds;
using paa::ChromosomeIndexLookup;
using paa::CoverageHMM;
using paa::Error;
using paa::FragmentInfo;
using paa::FragmentTable;
using paa::HaHaHMM;
using paa::HetTable;
using paa::Marker;
using paa::Progress;
using paa::PSDoc;
using paa::PSGraph;
using paa::PSHeat;
using paa::PSHSeries;
using paa::PSPage;
using paa::PSXYSeries;
using paa::PSXYSSeries;
using paa::Reference;
using paa::ThreadPool;

int main(int argc, char * argv[]) try {
  const std::string usage{"usage: phase_correlation ref in_dir chr ..."};
  if (--argc < 3) throw Error(usage);

  const unsigned int n_loci_window{5000};
  const unsigned int cutoff{n_loci_window * 20 / 5000};

  const Reference ref{argv[1]};
  const ChromosomeIndexLookup chr_lookup{ref};
  const string data_dir{argv[2]};
  argv += 2; argc -= 2;

  // Chromosomes to process
  const vector<string> chromosome_names{
    [&ref, &chr_lookup, argc, argv]() {
      vector<string> result;
      for (int a{0}; a != argc; ++a) {
        const std::string name{argv[a + 1]};
        chr_lookup[name];  // throws if name not found
        result.push_back(name);
      }
      return result;
    }()};
  std::cerr << "Running on chromosome"
            << (chromosome_names.size() > 1 ? "s" : "");
  for (const string name : chromosome_names) std::cerr << " " << name;
  std::cerr << endl;

  // Het input data
  const vector<HetTable> het_tables{
    [&chromosome_names, &chr_lookup, &data_dir]() {
      vector<HetTable> result;
      result.reserve(chromosome_names.size());
      for (const string name : chromosome_names) {
        const string table_name{data_dir + "/hets_" + name + ".txt"};
        result.emplace_back(table_name, chr_lookup);
      }
      return result;
    }()};
  const uint64_t n_samples{het_tables[0].n_samples()};

  using uVec = vector<uint64_t>;
  using Result1 = vector<uVec>;

  auto bin_fun = [n_loci_window, n_samples]
      (const HetTable & table, const uint64_t start) {
    Result1 result(4, uVec(n_samples));
    for (uint64_t h{start}; h != start + n_loci_window; ++h) {
      for (uint64_t sample{0}; sample != table.n_samples(); ++sample) {
        const unsigned int obs{table.unflipped(sample, h)};
        ++result[obs][sample];
      }
    }
    return result;
  };

  ThreadPool pool{160};
  vector<future<Result1>> futures1;
  uint64_t left_bin{0};
  ostringstream chr_ps;
  chr_ps << "0.25 lw\n";
  vector<uint64_t> chr_starts;
  vector<uint64_t> chr_stops;
  vector<uint64_t> chrs;
  for (const HetTable & left_table : het_tables) {
    chr_starts.push_back(left_bin);
    chr_ps << "np 1 1 1 c " << left_bin << " xc 0 yfc m "
           << left_bin << " xc 1 yfc l sp\n";
    chr_ps << "np 1 1 1 c 0 xfc " << left_bin << " yc m 1 xfc "
           << left_bin << " yc l sp\n";
    for (uint64_t left_start{0};
         left_start + n_loci_window < left_table.n_hets();
         left_start += n_loci_window, ++left_bin) {
      futures1.push_back(pool.run(bin_fun, cref(left_table), left_start));
      chrs.push_back(chr_stops.size());
    }
    chr_stops.push_back(left_bin);
  }
  for (uint64_t c{0}; c != chromosome_names.size(); ++c) {
    const uint64_t start{chr_starts[c]};
    const uint64_t stop{chr_stops[c]};
    chr_ps << 300.0 * (stop - start) / chrs.size() << " sf\n";
    chr_ps << "(" << replace_substring(chromosome_names[c], "chr", "")
           << ")" << stop << " xc " << start << " yc m gs jr show\n";
  }

  // plot
  PSDoc plots{"correlation"};
  plots.pdf(true);

  using Hist = PSHSeries<uint64_t, uint64_t>;
  unsigned int cover_bound((n_samples + 1) / 5);
  Hist * coverage_hist = Hist::create(
      plots, "Coverage;Single Allele Count;N",
      Bounds{0, 1.0 * cover_bound}, cover_bound);

  using dHist = PSHSeries<double, uint64_t>;
  PSPage corr_page{plots, "Correlations", "1 2"};
  dHist * agree_corr_hist = dHist::create(
      corr_page, "frac AA + BB;Correlation;N", Bounds{0, 1.0}, 100U);
  dHist * disagree_corr_hist = dHist::create(
      corr_page, "frac AB + BA;Correlation;N", Bounds{0, 1.0}, 100U);

  PSPage corr_dist_page{plots, "Correlations", "1 2"};
  const Marker marker{paa::circle(), 0.1, "0 0 0", 0.1, true};
  PSXYSeries corr_dist_agree{corr_dist_page, "Agree;dist;corr",
        marker, Bounds{0, 1.0 * ref.size(0), 0, 1}};
  PSXYSeries corr_dist_disagree{corr_dist_page, "Agree;dist;corr",
        marker, Bounds{0, 1.0 * ref.size(0), 0, 1}};
  const uint64_t max_count{n_loci_window / 10};
  Hist * segment_count_hist = Hist::create(
      plots, "Sample Segment Counts (cutoff=" + to_string(cutoff) + ");Count;N",
      Bounds{0, max_count}, static_cast<unsigned int>(max_count));

  PSXYSeries dan_plot{plots, "Dan;agree;disagree"};

  vector<Result1> counts;
  vector<vector<unsigned char>> has_allele(n_samples);
  for (auto & fut : futures1) {
    Result1 result{fut.get()};
    for (uint64_t sample{0}; sample != n_samples; ++sample) {
      const uint64_t a_count{result[HetTable::ObsA][sample] +
            result[HetTable::ObsAB][sample]};
      const uint64_t b_count{result[HetTable::ObsB][sample] +
            result[HetTable::ObsAB][sample]};
      unsigned char has{0};
      if (a_count > cutoff) has = has | 1;
      if (b_count > cutoff) has = has | 2;
      has_allele[sample].push_back(has);
      segment_count_hist->add_point(a_count);
      segment_count_hist->add_point(b_count);
    }
    counts.push_back(move(result));
  }

  cerr << "Generated " << counts.size() << " bins" << endl;

  struct Result2 {
    uint64_t left_bin{0};
    uint64_t right_bin{0};
    uint64_t joint[2][2]{{0, 0}, {0, 0}};
    uint64_t solo[2][2]{{0, 0}, {0, 0}};
  };

  auto bin_fun_2 = [&het_tables, &counts, &has_allele, n_samples]
      (const uint64_t lb, const uint64_t rb) noexcept {
    Result2 result;
    result.left_bin = lb;
    result.right_bin = rb;
    for (uint64_t sample{0}; sample != n_samples; ++sample) {
      const bool lA{has_allele[sample][lb] & 1};
      const bool lB{has_allele[sample][lb] & 2};
      const bool rA{has_allele[sample][rb] & 1};
      const bool rB{has_allele[sample][rb] & 2};
      if (lA) {
        ++result.solo[0][0];
        ++result.solo[0][1];
        if (rA) ++result.joint[0][0];
        if (rB) ++result.joint[0][1];
      }
      if (lB) {
        ++result.solo[1][0];
        ++result.solo[1][1];
        if (rA) ++result.joint[1][0];
        if (rB) ++result.joint[1][1];
      }
    }

    return result;
  };

  // Run segment-segment calc in parallel
  vector<future<Result2>> futures2;
  for (uint64_t lb{0}; lb != counts.size(); ++lb)
    for (uint64_t rb{0}; rb != counts.size(); ++rb)
      futures2.push_back(pool.run(bin_fun_2, lb, rb));

  // All loci plots
  PSHeat * correlations[2][2];
  char alleles[]{"AB"};
  for (const bool left_B : {false, true}) {
    for (const bool right_B : {false, true}) {
      (correlations[left_B][right_B] = PSHeat::create(
          plots, string() + "Common Alleles " + alleles[left_B] + " " +
          alleles[right_B] + ";Segment;Segment"))->expanded(true);
      correlations[left_B][right_B]->ps(chr_ps.str());
    }
  }

  // plots chromosome by chromosome
  using AlleleSeries = vector<vector<PSHeat *>>;
  vector<AlleleSeries> chr_correlations(
      chr_stops.size(),
      AlleleSeries(2, vector<PSHeat *>(2)));
  for (uint64_t c{0}; c != chromosome_names.size(); ++c) {
    for (const bool left_B : {false, true}) {
      for (const bool right_B : {false, true}) {
        chr_correlations[c][left_B][right_B] = PSHeat::create(
            plots, chromosome_names[c] + " Common Alleles " + alleles[left_B] +
            " " + alleles[right_B] + ";Segment;Segment");
        chr_correlations[c][left_B][right_B]->expanded(true);
      }
    }
  }

  // Get segment-segment results and plot them
  for (auto & fut : futures2) {
    const Result2 result{fut.get()};
    const uint64_t left_chr{chrs[result.left_bin]};
    const uint64_t right_chr{chrs[result.right_bin]};
    for (const bool left_B : {false, true}) {
      for (const bool right_B : {false, true}) {
        if (result.solo[left_B][right_B]) {
          const double value{bound(
              1.0 * result.joint[left_B][right_B] /
              result.solo[left_B][right_B], 0, 1)};
          correlations[left_B][right_B]->set_value(
              result.left_bin, result.right_bin, value);
          if (left_chr == right_chr) {
            chr_correlations[left_chr][left_B][right_B]->set_value(
                result.left_bin - chr_starts[left_chr],
                result.right_bin - chr_starts[right_chr], value);
          }
        }
      }
    }
  }

  // Get well-covered A or B loci, ignoring AB counts
  vector<vector<unsigned char>> loci;
  vector<unsigned int> loci_chr;
  vector<int> loci_pos;
  uint64_t max_cover{0};
  for (uint64_t t{0}; t != het_tables.size(); ++t) {
    const HetTable & table{het_tables[t]};
    for (uint64_t het{0}; het != table.n_hets(); ++het) {
      uint64_t n_single_alleles{0};
      vector<unsigned char> values;
      for (uint64_t sample{0}; sample != table.n_samples(); ++sample) {
        const unsigned char value{table.unflipped(sample, het)};
        values.push_back(value);
        n_single_alleles += value == 1 || value == 2;
      }
      if (max_cover < n_single_alleles) max_cover = n_single_alleles;
      coverage_hist->add_point(n_single_alleles);
      if (n_single_alleles > 20) {
        loci.push_back(move(values));
        loci_chr.push_back(t);
        loci_pos.push_back(table[het].position());
      }
    }
  }

  // Determine correlation between loci
  cerr << "Correlations" << endl;
  struct Result3 {
    uint64_t opportunity{0};
    uint64_t same{0};
    uint64_t different{0};
    unsigned int distance{0};
  };
  vector<future<vector<Result3>>> futures3;
  mutex plot_mutex;
  for (uint64_t ll_{0}; ll_ != loci.size(); ++ll_) {
    futures3.push_back(pool.run(
        [&loci, &loci_chr, &loci_pos,
         n_samples, agree_corr_hist, disagree_corr_hist,
         &plot_mutex]
        (const uint64_t ll) {
          vector<Result3> results;
          for (uint64_t rl{0}; rl != loci.size(); ++rl) {
            Result3 result;
            if (loci_chr[ll] < loci_chr[rl]) break;
            if (loci_chr[ll] > loci_chr[rl]) continue;
            for (uint64_t sam{0}; sam != n_samples ; ++sam) {
              const unsigned char left{loci[ll][sam]};
              if (left < 1 || left > 2) continue;
              const unsigned char right{loci[rl][sam]};
              if (right < 1 || right > 2) continue;
              ++result.opportunity;
              if (left == right) {
                ++result.same;
              } else {
                ++result.different;
              }
            }
            result.distance = loci_pos[ll] > loci_pos[rl] ?
                loci_pos[ll] - loci_pos[rl] : loci_pos[rl] - loci_pos[ll];
            results.push_back(result);
          }
          return results;
        }, ll_));
  }

  random_device rd;
  mt19937_64 mersenne(rd());
  uniform_real_distribution<double> dist{0, 1};
  function<double()> gen{bind(dist, std::ref(mersenne))};

  Progress progress{futures3.size(), 0.01, "loci"};
  uint64_t nnn{0};
  vector<double> vals;
  for (auto & future : futures3) {
    progress();
    vector<Result3> results{future.get()};
    for (const Result3 & result : results) {
      if (result.opportunity) {
        const double agree{1.0 * (result.same) / result.opportunity};
        const double disagree{
          1.0 * (result.different) / result.opportunity};
        agree_corr_hist->add_point(agree);
        disagree_corr_hist->add_point(disagree);
        if (result.opportunity >= 10 && result.distance > 1000) {
          if (nnn++ % 1 == 0) {
            corr_dist_agree.add_point(result.distance, agree);
            corr_dist_disagree.add_point(result.distance, disagree);
            vals.push_back(agree);
          }
        }
      }
    }
  }
  cerr << "vec size " << vals.size() << endl;
  sort(vals.begin(), vals.end());
  for (uint64_t i{0}; i != vals.size(); ++i) {
    dan_plot.add_point(1.0 * i / vals.size(), vals[i]);
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
