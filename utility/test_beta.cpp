//
// test_beta.cpp
//
// test beta distribution code
//
// Copyright 2020 Peter Andrews @ CSHL
//

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "beta.h"
#include "error.h"
#include "psplot.h"
#include "utility.h"

using std::cerr;
using std::cin;
using std::cout;
using std::endl;
using std::exception;
using std::istream;
using std::ostringstream;
using std::string;
using std::vector;

using paa::nunset;
using paa::Beta;
using paa::Beta4;
using paa::Bounds;
using paa::Error;
using paa::GSLRNG;
using paa::PSDoc;
using paa::PSPage;
using paa::PSHSeries;
using paa::PSXYSeries;

int main(int argc, char ** argv) try {
  const string usage{"usage: test_beta"};
  if (0) cerr << argc << argv[0];

  GSLRNG rng{0};
  const uint64_t n_bins{100};
  const double step{1.0 / n_bins / 10};
  const vector<double> alphas{0.5, 1.0, 2.0, 10.0, 100};
  const vector<double> betas{0.5, 1.0, 2.0, 10.0, 100};

  PSDoc plots{"beta"};
  plots.pdf(true);
  for (const double alpha : alphas) {
    for (const double beta : betas) {
      if (beta > alpha) continue;
      Beta beta_dist_simple{alpha, beta};
      Beta beta_dist{beta_dist_simple.mean(), beta_dist_simple.stdev(), true};
      ostringstream title;
      title << "Beta Distribution with alpha " << beta_dist.alpha()
            << ", beta " << beta_dist.beta() << ", mean " << beta_dist.mean()
            << ", stdev " << beta_dist.stdev() << ";x;p";
      PSXYSeries & plot{*PSXYSeries::create(
          plots, title.str(),
          Bounds{beta_dist.low(), beta_dist.high(), 0, nunset(0.0)})};
      for (double val{step}; val < 1; val += step)
        plot.add_point(val, beta_dist(val));
      const double low{2};
      const double high{4};
      const double scale{high - low};
      Beta4 scaled_beta{
        low + beta_dist.mean() * scale, beta_dist.stdev() * scale,
            low, high, true};
      ostringstream scaled_title;
      scaled_title << "Beta Distribution with alpha " << scaled_beta.alpha()
                   << ", beta " << scaled_beta.beta()
                   << ", mean " << scaled_beta.mean()
                   << ", stdev " << scaled_beta.stdev() << ";x;p";
      PSXYSeries & scaled_plot{*PSXYSeries::create(
          plots, scaled_title.str(),
          Bounds{scaled_beta.low(), scaled_beta.high(), 0, nunset(0.0)})};
      for (double val{low + step * scale};
           val + step < high; val += step * scale)
        scaled_plot.add_point(val, scaled_beta(val));
      using Hist = PSHSeries<double, uint64_t>;
      Hist & hist{*Hist::create(
          plots, scaled_title.str(), n_bins,
          Bounds{scaled_beta.low(), scaled_beta.high()})};
      const uint64_t n_samples{10000000};
      for (uint64_t s{0}; s != n_samples; ++s)
        hist.add_point(scaled_beta(rng));
    }
  }

  return 0;
} catch (Error & e) {
  cerr << "paa::Error:" << endl;
  cerr << e.what() << endl;
  return 1;
} catch(exception & e) {
  cerr << e.what() << endl;
  return 1;
} catch(...) {
  cerr << "Some exception was caught." << endl;
  return 1;
}
