//
// simulate_seating.cpp
//
// Determine 'seating chart' from distances between classmates
// actually has an application in spatial transcriptomics
//
// Copyright 2018 Peter Andrews @ CSHL
//

#include <algorithm>
#include <exception>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <functional>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "error.h"
#include "psplot.h"

using std::bind;
using std::cout;
using std::cerr;
using std::endl;
using std::exception;
using std::function;
using std::hex;
using std::ifstream;
using std::istringstream;
using std::min;
using std::mt19937_64;
using std::ofstream;
using std::ostream;
using std::ostringstream;
using std::random_device;
using std::setfill;
using std::setw;
using std::string;
using std::uniform_real_distribution;
using std::vector;

using paa::Bounds;
using paa::Error;
using paa::Marker;
using paa::PSDoc;
using paa::PSGraph;
using paa::PSXYSeries;

// Nasty global variables, set in main()
function<double()> gen;
uint64_t n_seats{0};
double error{0.0};
bool randomize_pos{true};

template <class Val> inline Val sqr(const Val val) { return val * val; }

// Encapsulates a grid layout
struct SeatGrid {
  SeatGrid() {}
  uint64_t row(const uint64_t s) const { return s / n_rows; }
  uint64_t col(const uint64_t s) const { return s % n_rows; }
  uint64_t seat(const uint64_t r, const uint64_t c) { return r * n_rows + c; }
  void add_seat(double & x, double & y) {
    y = 1.0 * row(n_points) / n_rows + off;
    x = 1.0 * col(n_points) / n_rows + off;
    ++n_points;
  }
  uint64_t n_rows{[]() {
      const uint64_t result{static_cast<uint64_t>(sqrt(n_seats))};
      return result + (sqr(result) == n_seats ? 0 : 1);
    }()};
  double off{0.5 / n_rows};
  uint64_t n_points{0};
} grid;

// Creates randomized colors, and formats as float or hex RGB
struct Color {
  string to_hex() const {
    ostringstream result;
    result << hex << setfill('0') << setw(2)
           << static_cast<unsigned int>(r * 256)
           << hex << setfill('0') << setw(2)
           << static_cast<unsigned int>(g * 256)
           << hex << setfill('0') << setw(2)
           << static_cast<unsigned int>(b * 256);
    return result.str();
  }
  double r{gen()};
  double g{gen()};
  double b{gen()};
};
ostream & operator<<(ostream & out, const Color & c) {
  return out << c.r << " " << c.g << " " << c.b;
}

// Seats placed in a unit square, in a grid or randomly
struct Seat {
  Seat() {
    if (randomize_pos) {
      x = gen();
      y = gen();
    } else {
      grid.add_seat(x, y);
    }
  }
  double distance_to(const Seat & seat) const {
    return sqrt(sqr(x - seat.x) + sqr(y - seat.y));
  }
  double noisy_distance_to(const Seat & seat) const {
    return distance_to(seat) * (1.0 + error * (0.5 - gen())) +
        error * 0.1 * gen();
  }
  double x{0};
  double y{0};
  Color color{};
};

// Matrix with two run-time dimensions
template <class Val>
struct Matrix {
  Matrix(const uint64_t I_, const uint64_t J_) : I{I_}, J{J_}, data(I * J) { }
  const Val * operator[](const uint64_t i) const { return &this->data[i * J]; }
  Val * operator[](const uint64_t i) { return &this->data[i * J]; }
  const uint64_t I;
  const uint64_t J;
  vector<Val> data;
};

int main(int argc, char* argv[]) try {
  // Read command line arguments
  if (--argc != 3)
    throw Error("usage: simulate_seating n_seats error randomize");
  n_seats = strtoul(argv[1], nullptr, 10);
  error = strtod(argv[2], nullptr);
  randomize_pos = strtoul(argv[3], nullptr, 10);
  if (!randomize_pos) grid = SeatGrid();

  // Set up uniform real random number generator
  random_device rd;
  auto mersenne = mt19937_64(rd());
  gen = bind(uniform_real_distribution<double>(0, 1), mersenne);

  // Place seats randomly
  const vector<Seat> seats{[]() {
      vector<Seat> result;
      result.reserve(n_seats);
      while (result.size() != n_seats) result.emplace_back();
      return result;
    }()};

  // Get full seat distance matrix with noise
  // Upper triangular = noisy (changes to reconstructed values later on)
  // lower triangular = actual
  Matrix<double> distances{[&seats]() {
      Matrix<double> result(seats.size(), seats.size());
      for (uint64_t i{0}; i != seats.size(); ++i) {
        for (uint64_t j{i + 1}; j != seats.size(); ++j) {
          result[i][j] = seats[i].noisy_distance_to(seats[j]);
          result[j][i] = seats[i].distance_to(seats[j]);
        }
      }
      return result;
    }()};

  // Output DOT for GraphViz neato rendering and make simple calculations
  ofstream dot{"seating.dot"};
  dot << "graph seating {\n"
      // << "epsilon=\"0.00001\"\n"
      // << "maxiter=\"100000\"\n"
      << "edge [style=invis];\n"
      << "node [shape=circle, style=filled, label=\"\", width=\""
      << 4 * sqrt(1.0 / n_seats) << "\"];\n";
  for (uint64_t i{0}; i != seats.size(); ++i) {
    const string hex{seats[i].color.to_hex()};
    dot << i << " [color=\"#" << hex << "\", fillcolor=\"#" << hex << "\"]\n";
  }
  const double theory_neighbor_distance{sqrt(1.0 / n_seats)};
  const double distance_limit{min(0.5, 10 * theory_neighbor_distance)};
  double average_distance{0};
  uint64_t n_distances{0};
  uint64_t n_edges{0};
  for (uint64_t i{0}; i != seats.size(); ++i) {
    for (uint64_t j{i + 1}; j != seats.size(); ++j) {
      const double distance{distances[i][j]};
      average_distance += static_cast<double>(distance);
      ++n_distances;
      if (distance < distance_limit) {
        dot << i << " -- " << j << " [len=\"" << 10 * distance << "\"];\n";
        ++n_edges;
      }
    }
  }
  dot << "}\n";
  dot.close();

  // Render reconstructed seating chart as fully specified DOT output
  if (system("neato -Txdot -o seating_reconstruction.dot seating.dot") == -1)
    cerr << "Problem rendering reconstructed seating chart" << endl;

  // Read in reconstructed layout - sloppy parsing, but works
  ifstream reconstruction{"seating_reconstruction.dot"};
  if (!reconstruction) throw Error("Problem opening seating reconstruction");
  string line;
  uint64_t n_read{0};
  vector<double> x(seats.size());
  vector<double> y(seats.size());
  while (getline(reconstruction, line)) {
    uint64_t id;
    string label_str;
    string width_str;
    string color_str;
    string fill_str;
    string pos_str;
    string height_str;
    istringstream line_stream{line.c_str()};
    if (!(line_stream >> id >> label_str >> width_str
          >> color_str >> fill_str >> pos_str >> height_str)) continue;
    if (label_str.substr(0, 6) != "[label") continue;
    if (fill_str.substr(0, 3) == "pos") pos_str = fill_str;
    if (pos_str.substr(0, 3) != "pos") continue;
    if (id != n_read)
      throw Error("Seat reordering in reconstruction") << id << n_read << line;
    istringstream pos_stream{pos_str.substr(5).c_str()};
    char comma;
    pos_stream >> x[n_read] >> comma >> y[n_read];
    const double tmp = x[n_read];
    x[n_read] = -y[n_read] + 800;
    y[n_read] = -tmp + 800;
    if (!pos_stream) throw Error("Problem parsing pos string") << pos_str;
    if (false) cout << id << '\t' << x[n_read] << '\t' << y[n_read] << '\n';
    ++n_read;
  }
  if (n_read != seats.size())
    throw Error("Wrong number of reconstructed seats read in")
        << n_read << seats.size();

  // Plots
  PSDoc plots{"seating", "seating", 600, 600};
  plots.pdf(true);
  const Marker marker{paa::circle(), 0.3, "0 0 0", 1, true, "0 0 0"};

  // Plot true vs perturbed distances
  const uint64_t points_sampling{n_distances /
        min(static_cast<uint64_t>(10000), n_distances)};
  PSGraph truth_vs_perturbed_graph{plots,
        "Seating Chart Added Noise;Actual Distance;Perturbed Distance"};
  PSXYSeries truth_vs_perturbed{truth_vs_perturbed_graph, marker};
  uint64_t point_no{0};
  for (uint64_t i{0}; i != seats.size(); ++i)
    for (uint64_t j{i + 1}; j != seats.size(); ++j)
      if ((++point_no % points_sampling) == 0)
        truth_vs_perturbed.add_point(distances[j][i], distances[i][j]);

  ostringstream limit_stream;
  limit_stream << "newpath 1 0 0 c 3 lw 0 xfc " << distance_limit << " yc m "
               << "1 xfc " << distance_limit << " yc l stroke "
               << "20 sf 0.3 xfc " << distance_limit / 2 << " yc 5 sub m ";
  // << "(<-- these distances are used as input) s";
  truth_vs_perturbed_graph.ps(limit_stream.str());

  // Get distances from reconstructed positions
  double scaled_average_distance{0};
  for (uint64_t i{0}; i != seats.size(); ++i) {
    for (uint64_t j{i + 1}; j != seats.size(); ++j) {
      distances[i][j] = sqrt(sqr(x[i] - x[j]) + sqr(y[i] - y[j]));
      scaled_average_distance += distances[i][j];
    }
  }
  const double distance_scale(scaled_average_distance / average_distance);
  double total_distance_deviation{0};
  for (uint64_t i{0}; i != seats.size(); ++i) {
    x[i] /= distance_scale;
    y[i] /= distance_scale;
    for (uint64_t j{i + 1}; j != seats.size(); ++j)
      total_distance_deviation +=
          fabs((distances[i][j] /= distance_scale) - distances[j][i]);
  }

  // Output some info to the terminal
  cout << "average distance " << average_distance / n_distances << endl
       << "theory neighbor distance " << theory_neighbor_distance << endl
       << "distance limit " << distance_limit << endl
       << "edges used per seat " << 1.0 * n_edges / n_seats << endl
       << "distance scale " << distance_scale << endl
       << "Average distance deviation "
       << total_distance_deviation / n_distances << endl;

  // Plot true vs reconstructed distances
  PSXYSeries truth_vs_reconstructed{plots,
        "Seating Chart Reconstruction;Actual Distance;Reconstructed Distance",
        marker};
  for (uint64_t i{0}; i != seats.size(); ++i)
    for (uint64_t j{i + 1}; j != seats.size(); ++j)
      if ((++point_no % points_sampling) == 0)
        truth_vs_reconstructed.add_point(distances[j][i], distances[i][j]);

  // Plot the ground truth seating chart
  PSGraph truth{plots, "Seating Chart Ground Truth;X;Y",
        Bounds{-0.05, 1.05, -0.05, 1.05}};
  const double seat_size{100 * sqrt(1.0 / n_seats)};
  ostringstream ps;
  for (const Seat & seat : seats)
    ps << seat.color << " c " << "np " << seat.x << " " << seat.y << " gc "
       << seat_size << " 0 360 arc fill" << "\n";
  truth.ps(ps.str());

  // Plot reconstructed seating chart
  auto x_mm = minmax_element(x.begin(), x.end());
  auto y_mm = minmax_element(y.begin(), y.end());
  PSGraph reconstruct{plots, "Seating Chart Reconstruction;X;Y",
        Bounds{*x_mm.first - 0.05, *x_mm.second + 0.05,
          *y_mm.first - 0.05, *y_mm.second + 0.05}};
  ps.str("");
  for (uint64_t s{0}; s != x.size(); ++s)
    ps << seats[s].color << " c " << "np " << x[s] << " " << y[s] << " gc "
       << seat_size << " 0 360 arc fill" << "\n";
  reconstruct.ps(ps.str());

  if (randomize_pos) return 0;

  // Plot average between-row distances to gauge warp
  const Marker bigger_marker{paa::circle(), 1, "0 0 0", 1, true, "0 0 0"};
  PSXYSeries row_distances{plots,
        "Seating Chart Reconstruction Warp;Row;Average Adjacent Row Distance",
        bigger_marker};
  for (uint64_t r{0}; r + 1 != grid.n_rows; ++r) {
    double total_distance{0};
    uint64_t n_added{0};
    for (uint64_t c{0}; c != grid.n_rows; ++c) {
      const uint64_t i{grid.seat(r, c)};
      if (i >= n_seats) continue;
      const uint64_t j{grid.seat(r + 1, c)};
      if (j >= n_seats) continue;
      total_distance += distances[i < j ? i : j][i < j ? j : i];
      ++n_added;
    }
    if (n_added) row_distances.add_point(r, total_distance / n_added);
  }

  return 0;
} catch (Error & e) {
  cerr << "Error:" << endl;
  cerr << e.what() << endl;
  return 1;
} catch (exception & e) {
  cerr << "exception" << endl;
  cerr << e.what() << endl;
  return 1;
} catch (...) {
  cerr << "unknown exception was caught" << endl;
  return 1;
}
