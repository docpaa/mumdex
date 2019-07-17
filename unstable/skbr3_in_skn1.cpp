//
// skbr3_in_skn1
//
// Detect known profile when added to flat profile
//
// Copyright 2017 Peter Andrews @ CSHL
//

#include <algorithm>
#include <array>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <random>
#include <sstream>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include "cn.h"
#include "error.h"
#include "files.h"
#include "mumdex.h"
#include "numerical.h"
#include "psplot.h"
#include "stats.h"
#include "threads.h"
#include "utility.h"

using std::array;
using std::binomial_distribution;
using std::cerr;
using std::cout;
using std::endl;
using std::exception;
using std::istringstream;
using std::lock_guard;
using std::make_unique;
using std::max;
using std::min;
using std::mt19937_64;
using std::mutex;
using std::numeric_limits;
using std::ofstream;
using std::ostringstream;
using std::pair;
using std::poisson_distribution;
using std::random_device;
using std::setprecision;
using std::string;
using std::to_string;
using std::uniform_int_distribution;
using std::unique_lock;
using std::unique_ptr;
using std::vector;

using paa::dne;
using paa::Bin;
using paa::Bounds;
using paa::CN_Bin;
using paa::CN_Bins;
using paa::Error;
using paa::FinestBins;
using paa::GoldenMinimizer;
using paa::LinearRegression;
using paa::MappedVector;
using paa::Marker;
using paa::NormalParams;
using paa::Progress;
using paa::PSDoc;
using paa::PSGraph;
using paa::PSHSeries;
using paa::PSPage;
using paa::PSXYSeries;
using paa::Reference;
using paa::ThreadPool;

#if 0
class Level {
 public:
  Level(const double value_, const bool in_segment_) :
      data(value_, in_segment_) { }
  double value() const { return data.first; }
  bool in_segment() const { return data.second; }

 private:
  pair<float, bool> data;
};
#else
class Level {
 public:
  Level(const double value_, const bool in_segment_) :
      data((in_segment_ ? -1 : 1) * value_) { }
  double value() const { return fabs(data); }
  bool in_segment() const { return data < 0; }

 private:
  float data;
};
#endif

class LogL {
 public:
  LogL(const vector<double> & pn_, const vector<double> & pc_,
       const vector<unsigned int> & cs_) : pn{pn_}, pc{pc_}, cs{cs_} { }

  double operator()(const double alpha) const {
    double result{0};
    for (unsigned int r{0}; r != cs.size(); ++r) {
      result += cs[r] * log((1 - alpha) * pn[r] + alpha * pc[r]);
    }
    return -result;
  }

 private:
  const vector<double> & pn;
  const vector<double> & pc;
  const vector<unsigned int> & cs;
};

int main(int argc, char* argv[]) try {
  // Check usage
  --argc;
  if (argc != 6) {
    throw Error("usage: skbr3_in_skn1 ref bins.txt skn1_counts skbr3_counts "
                "counts_per_bin skbr3_level");
  }

  //
  // Process command line arguments
  //

  // Genome reference
  const Reference ref{argv[1]};

  // Bins
  const vector<Bin> all_bins{load_bins(argv[2], ref)};
  const vector<Bin> bins{[&all_bins]() {
      vector<Bin> result;
      for (const Bin & bin : all_bins)
        if (!bin.bad()) result.push_back(bin);
      return result;
    }()};

  // Load processed copy number
  const CN_Bins skn1_data{argv[3]};
  const CN_Bins skbr3_data{argv[4]};
  const double cpb{atof(argv[5])};

  // Amount of skbr3 to add
  const string level_str{argv[6]};
  const double level{atof(argv[6])};
  if (level >= 1) throw Error("Assumes level less than 1");

  // Autosome bins
  const unsigned int x_chr{ref.find_x_chromosome()};
  const unsigned int n_bins{[&bins, x_chr]() {
      unsigned int result{0};
      for (unsigned int b{0}; b != bins.size(); ++b) {
        if (bins[b].chromosome() == x_chr) break;
        ++result;
      }
      return result;
    }()};

  // Autosome totals
  const double skn1_total{[&skn1_data, n_bins]() {
      double result{0};
      for (unsigned int b{0}; b != n_bins; ++b) {
        result += skn1_data[b].norm_count();
      }
      return result;
    }()};
  const double skbr3_total{[&skbr3_data, n_bins]() {
      double result{0};
      for (unsigned int b{0}; b != n_bins; ++b) {
        result += skbr3_data[b].norm_count();
      }
      return result;
    }()};
  cout << "SKN1 input cpb " << skn1_total / n_bins
       << " SKBR3 input cpb " << skbr3_total / n_bins << endl;

  // Randomization
  random_device rd;
  auto mersenne = mt19937_64(rd());

  vector<double> zero_mixtures;
  for (const bool add_skbr3 : {false, true}) {
    // Downsample rate
    const double skn1_cpb{skn1_total / n_bins};
    const double normal_downsample{cpb / skn1_cpb};
    const double normal_mixed_downsample{normal_downsample *
          (add_skbr3 ? (1 - level) : 1)};
    const double skbr3_cpb{skbr3_total / n_bins};
    const double skbr3_downsample{cpb / skbr3_cpb *
           (add_skbr3 ? level : 0)};

    const unsigned n_iter{add_skbr3 ? 1U : 100U};

    for (unsigned int i{0}; i != n_iter; ++i) {
      // cout << i << endl;
      // Generate simulated data
      vector<double> skbr3_counts;
      vector<double> skn1_counts;
      vector<unsigned int> counts;
      const double normal_rate{normal_mixed_downsample};
      const double skbr3_rate{skbr3_downsample};
      const CN_Bins * all_data[2]{&skn1_data, &skbr3_data};
      // vector<double> measures[2];
      for (unsigned int b{0}; b != n_bins; ++b) {
        skbr3_counts.push_back(skbr3_data[b].norm_count() * cpb / skbr3_cpb);
        counts.push_back(0);
        skn1_counts.push_back(skn1_data[b].norm_count() * cpb / skn1_cpb);
        for (const bool is_skbr3 : {false, true}) {
          const CN_Bins * data{all_data[is_skbr3]};
          const double rate{is_skbr3 ? skbr3_rate : normal_rate};
          if (rate > 0) {
            binomial_distribution<unsigned int> dist{
              static_cast<unsigned int>((*data)[b].norm_count()), rate};
            counts[b] += dist(mersenne);
          }
        }
      }

      // Better determination of mixture
      double sxy{0};
      double syy{0};
      for (unsigned int b{0}; b != counts.size(); ++b) {
        const double x{counts[b] - skn1_counts[b]};
        const double y{skn1_counts[b] - skbr3_counts[b]};
        sxy += x * y;
        syy += y * y;
      }
      const double mix{-sxy / syy};
      if (!add_skbr3) {
        zero_mixtures.push_back(mix);
      }

      if (add_skbr3) {
        NormalParams normal{zero_mixtures};

        // const double mix_ci{1.96 / sqrt(syy)};
        const double mix_ci{1.96 * normal.stdev};
        cout << "mixture is " << mix << " +/- " << mix_ci << endl;

        // Plot simulated data counts vs skbr3 counts
        const std::string dark{"0 0 0"};
        const Marker dark_circle_marker{paa::circle(), 0.2, dark, 0.2, true};
        PSDoc ps{"skbr3_in_skn1", "skbr3_in_skn1"};
        ostringstream sim_details;
        sim_details << "SKBR3 in SKN1 simulation: bin count = " << cpb
                    << ", SKBR3 level = " << level_str;
        PSGraph counts_graph{ps, sim_details.str() +
              ";SKBR3 bin count;Simulated bin count"};
        PSXYSeries counts_series{counts_graph, dark_circle_marker};

        for (unsigned int b{0}; b != counts.size(); ++b) {
          counts_series.add_point(skbr3_counts[b], counts[b]);
        }

        if (0) {
          // Linear fit of counts vs skbr3 counts
          const LinearRegression linear{skbr3_counts, counts};
          cout << "Slope " << linear.slope << " +/- " << linear.slope_error
               << " intercept " << linear.intercept
               << " +/- " << linear.intercept_error << endl;
        }

        // Draw fit line on graph
        auto minmax = minmax_element(skbr3_counts.begin(), skbr3_counts.end());
        const double low_x{*minmax.first};
        const double high_x{*minmax.second};
        ostringstream line_ps;
        line_ps << "newpath 1 0 0 c 2 setlinewidth "
                << low_x << " "
                << cpb - mix * cpb << " gc m "
                << high_x << " "
                << (high_x - cpb) * mix + cpb << " gc l stroke "
                << "newpath 0 0 1 c 1 setlinewidth "
                << low_x << " "
                << cpb - (mix - mix_ci) * cpb << " gc m "
                << high_x << " "
                << (high_x - cpb) * (mix + mix_ci) + cpb << " gc l stroke "
                << "newpath 0 0 1 c 1 setlinewidth "
                << low_x << " "
                << cpb - (mix + mix_ci) * cpb << " gc m "
                << high_x << " "
                << (high_x - cpb) * (mix - mix_ci) + cpb << " gc l stroke ";

        counts_graph.ps(line_ps.str());

        // Alex's first principles method
        cout << "Alex method" << endl;
        vector<double> pn;
        vector<double> pc;
        vector<unsigned int> cs;
        for (unsigned int b{0}; b != counts.size(); ++b) {
          if (b == 0 ||
              dne(skn1_data[b].seg_ratio(), skn1_data[b - 1].seg_ratio()) ||
              dne(skbr3_data[b].seg_ratio(), skbr3_data[b - 1].seg_ratio())) {
            pn.push_back(0);
            pc.push_back(0);
            cs.push_back(0);
          }
          pn.back() += skn1_data[b].ratio();
          pc.back() += skbr3_data[b].ratio();
          cs.back() += counts[b];
        }
        const double tn{accumulate(pn.begin(), pn.end(), 0.0)};
        const double tc{accumulate(pc.begin(), pc.end(), 0.0)};
        for (unsigned int r{0}; r != pn.size(); ++r) {
          pn[r] /= tn;
          pc[r] /= tc;
        }

        const LogL likelihood{pn, pc, cs};
        const GoldenMinimizer minimizer{likelihood, 0, 1};
        const double alpha{minimizer.min()};
        double err{0};
        for (unsigned int r{0}; r != pn.size(); ++r) {
          err += cs[r] * sqr((pc[r] - pn[r]) /
                             ((1 - alpha) * pn[r] + alpha * pc[r]));
        }
        err = 1 / sqrt(err);
        cout << "Minimum is at " << alpha << " +/- " << err << endl;
      }
    }
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


