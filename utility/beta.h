//
// beta.h
//
// code using beta distribution
//
// copyright 2020 Peter Andrews
//

#ifndef PAA_BETA_H_
#define PAA_BETA_H_

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <cmath>
#include <mutex>
#include <utility>

#include "error.h"
#include "utility.h"

namespace paa {

class GSLRNG {
 public:
  explicit GSLRNG(const uint64_t seed) : gen{gsl_rng_alloc(gsl_rng_mt19937)} {
    if (!gen) throw Error("Problem creating gsl generator");
    gsl_rng_set(gen, seed);
  }
  ~GSLRNG() { gsl_rng_free(gen); }
  operator gsl_rng * () const { return gen; }
  GSLRNG(const GSLRNG &) = delete;
  GSLRNG & operator=(const GSLRNG &) = delete;
  void lock() { mutex.lock(); }
  void unlock() { mutex.unlock(); }

 private:
  gsl_rng * gen;
  std::mutex mutex{};
};

template<class Val>
class BetaT {
 public:
  using Params = std::pair<Val, Val>;
  BetaT(const Val alpha__, const Val beta__) :
      alpha_{alpha__}, beta_{beta__} {
        if (alpha_ <= 0) throw Error("Nonpositive alpha in Beta distribution");
        if (beta_ <= 0) throw Error("Nonpositive beta in Beta distribution");
      }
  explicit BetaT(const Params params) :
      BetaT{params.first, params.second} {}
  BetaT(const Val mean_, const Val stdev_, const bool) :
      alpha_{get_alpha(mean_, stdev_)},
      beta_{get_beta(mean_, stdev_)} {
        if (stdev_ >= max_stdev(mean_))
          throw Error("Maximum stdev for basic Beta distribution is")
              << max_stdev(mean_);
      }
  Val alpha() const { return alpha_; }
  Val beta() const { return beta_; }
  Val low() const { return 0.0; }
  Val high() const { return 1.0; }
  Val range() const { return 1.0; }
  Val middle() const { return 0.5; }
  static Val get_alpha(const Val mean_, const Val stdev_) {
    return mean_ * (mean_ * (1 - mean_) / sqr(stdev_) - 1);
  }
  static Val get_beta(const Val mean_, const Val stdev_) {
    return (1 - mean_) * (mean_ * (1 - mean_) / sqr(stdev_) - 1);
  }
  static Val max_stdev(const Val mean_) {
    return sqrt(mean_ * (1 - mean_));
  }
  Val mean() const { return alpha_ / (alpha_ + beta_); }
  Val stdev() const {
    return sqrt((alpha_ * beta_) /
                (sqr(alpha_ + beta_) * (alpha_ + beta_ + 1)));
  }
  Val operator()(const Val val) const {
    return gsl_ran_beta_pdf(val, alpha_, beta_);
  }
  Val operator()(gsl_rng * gen) const {
    return gsl_ran_beta(gen, alpha_, beta_);
  }

 private:
  Val alpha_;
  Val beta_;
};
using Beta = BetaT<double>;
using fBeta = BetaT<float>;

template<class Val>
class Beta4T {
 public:
  Beta4T(const Val alpha__, const Val beta__,
        const Val low__, const Val high__) :
      beta_dist{alpha__, beta__}, low_{low__}, high_{high__} {
        if (mean() < low())
          throw Error("Mean less than low in Beta distribution");
        if (mean() > high())
          throw Error("Mean greater than high in Beta distribution");
      }
  Beta4T(const Val mean_, const Val stdev_,
        const Val low__, const Val high__,
        const bool) :
      beta_dist{get_alpha_beta(mean_, stdev_, low__, high__)},
      low_{low__}, high_{high__} {}

  Val alpha() const { return beta_dist.alpha(); }
  Val beta() const { return beta_dist.beta(); }
  Val low() const { return low_; }
  Val high() const { return high_; }
  Val range() const { return (high_ - low_); }
  Val middle() const { return (high_ - low_) / 2; }
  Val mean() const { return low_ + scale * beta_dist.mean(); }
  Val stdev() const { return scale * beta_dist.stdev(); }
  static Val max_stdev(const Val mean_, const Val low__, const Val high__) {
    const Val sf{high__ - low__};
    const Val m{(mean_ - low__) / sf};
    return sf * BetaT<Val>::max_stdev(m);
  }
  Val operator()(const Val val) const {
    return static_cast<Val>(beta_dist((val - low_) / scale)) / scale;
  }
  Val operator()(gsl_rng * gen) const {
    return low_ + scale * static_cast<Val>(beta_dist(gen));
  }

 private:
  static typename BetaT<Val>::Params get_alpha_beta(
      const Val mean_, const Val stdev_,
      const Val low__, const Val high__) {
    const Val sf{high__ - low__};
    const Val m{(mean_ - low__) / sf};
    const Val sd{stdev_ / sf};
    return typename BetaT<Val>::Params{
      BetaT<Val>::get_alpha(m, sd), BetaT<Val>::get_beta(m, sd)};
  }

  BetaT<Val> beta_dist;
  Val low_;
  Val high_;
  Val scale{high_ - low_};
};
using Beta4 = Beta4T<double>;
using fBeta4 = Beta4T<float>;

}  // namespace paa



#endif  // PAA_BETA_H_
