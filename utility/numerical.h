//
// numerical.h
//
// function minimization routines
//
// copyright 2017 Peter Andrews
//

#ifndef PAA_NUMERICAL_H_
#define PAA_NUMERICAL_H_

#include <cmath>
#include <limits>

namespace paa {

class GoldenMinimizer {
 public:
  template <class F>
  GoldenMinimizer(const F & func, double a, double b,
                  const double tol =
                  std::sqrt(std::numeric_limits<double>::epsilon())) {
    const double gr{(sqrt(5) + 1) / 2};
    double c{b - (b - a) / gr};
    double d{a + (b - a) / gr};
    while (fabs(c - d) > tol * (fabs(a) + fabs(b))) {
      if (func(c) < func(d)) {
        b = d;
      } else {
        a = c;
      }
      c = b - (b - a) / gr;
      d = a + (b - a) / gr;
    }
    min_ = (b + a) / 2;
  }

  double min() const { return min_; }

 private:
  double min_{0};
};

}  // namespace paa

#endif  // PAA_NUMERICAL_H_
