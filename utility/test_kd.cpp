//
// test_kd.cpp
//
// test k-d tree implementation
//
// Copyright 2020 Peter Andrews @ CSHL
//

#include <iostream>
#include <random>
#include <string>
#include <vector>

#include "error.h"
#include "kdtree.h"
#include "utility.h"

using std::cerr;
using std::cout;
using std::endl;
using std::exception;
using std::istream;
using std::mt19937_64;
using std::random_device;
using std::string;
using std::uniform_real_distribution;
using std::vector;

using paa::commas;
using paa::Error;
using paa::KDTreeT;
using paa::PointD;
using paa::Timer;

// For random points between 0 and 1
uniform_real_distribution<double> pos_gen{0, 1};
mt19937_64 mersenne{random_device()()};

// Return a random point using pos_gen
template <class Type, uint64_t n_dimensions>
PointD<Type, n_dimensions> random_point() {
  PointD<Type, n_dimensions> point;
  for (uint64_t dimension{0}; dimension != n_dimensions; ++dimension)
    point[dimension] = pos_gen(mersenne);
  return point;
}

// Make sure the tree search structure is working
template <uint64_t n_dimensions>
void test_kd(const uint64_t n_points, const uint64_t n_trials) {
  // Give run info
  cerr << "Testing " << commas(n_points) << " points in "
       << n_dimensions << " dimensions for "
       << commas(n_trials) << " trials" << endl;

  // Setup points
  using Point = PointD<double, n_dimensions>;
  const vector<Point> points{[n_points]() {
      vector<Point> result;
      for (uint64_t p{0}; p != n_points; ++p)
        result.push_back(random_point<double, n_dimensions>());
      return result;
    }()};

  using KDTree = KDTreeT<uint32_t, double, n_dimensions>;
  constexpr double big{KDTree::big};

  // Find closest point, dumb method
  Timer timer;
  double average_distance{0};
  for (uint64_t trial{0}; trial != n_trials; ++trial) {
    const Point point{random_point<double, n_dimensions>()};
    double best_distance2{big};
    for (uint64_t p{0}; p != n_points; ++p) {
      const double distance2{point.distance2(points[p])};
      if (distance2 < best_distance2) best_distance2 = distance2;
    }
    average_distance += sqrt(best_distance2);
  }
  cerr << "Brute force method takes " << timer.seconds()
       << " seconds with average distance of " << average_distance / n_trials
       << endl;

  // KD tree method
  average_distance = 0;
  timer.reset();
  const KDTree tree{points};
  const double setup{timer.seconds()};
  timer.reset();
  for (uint64_t trial{0}; trial != n_trials; ++trial) {
    const Point point{random_point<double, n_dimensions>()};
    // const Point closest{tree.find_closest(point).point};
    const Point closest{points[tree.find_n_closest(point, 10).back()]};
    average_distance += sqrt(point.distance2(closest));
  }
  cerr << "KDTree method takes " << setup << " seconds for setup and "
       << timer.seconds() << " seconds for lookups with average distance of "
       << average_distance / n_trials << endl;

  // Check for different results
  cerr << "Checking for mismatches in results" << endl;
  for (uint64_t trial{0}; trial != n_trials; ++trial) {
    const Point point{random_point<double, n_dimensions>()};
    double best_distance2{big};
    uint64_t best_point{0};
    for (uint64_t p{0}; p != n_points; ++p) {
      const double distance2{point.distance2(points[p])};
      if (distance2 < best_distance2) {
        best_distance2 = distance2;
        best_point = p;
      }
    }
    const Point closest{points[tree.find_closest(point)]};
    if (closest != points[best_point]) throw Error("Closest mismatch");
  }
}

int main(int argc, char ** argv) try {
  // Process command line arguments
  const string usage{"usage: test_kd n_dimensions n_points n_trials"};
  if (--argc != 3) throw Error(usage);
  const uint64_t n_dimensions = atoi(argv[1]);
  const uint64_t n_points = atoi(argv[2]);
  const uint64_t n_trials = atoi(argv[3]);

  // Use compiled implementation that is specific to selected dimensionality
  switch (n_dimensions) {
    case 1:
      test_kd<1>(n_points, n_trials);
      break;
    case 2:
      test_kd<2>(n_points, n_trials);
      break;
    case 3:
      test_kd<3>(n_points, n_trials);
      break;
    case 4:
      test_kd<4>(n_points, n_trials);
      break;
    default:
      throw Error("Not (trivially) set up for more than 4 dimensions");
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
