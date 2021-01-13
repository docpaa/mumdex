//
// grid.h
//
// Classes that use a grid for lookup
//
// Copyright 2020 Peter Andrews @ CSHL
//

#ifndef PAA_UTILITY_GRID_H_
#define PAA_UTILITY_GRID_H_

#include <algorithm>
#include <iostream>
#include <utility>
#include <vector>

#include "image.h"

namespace paa {

// A base class for grid-based access
template <class Cell>
class GridBase {
 public:
  using Data = std::vector<Cell>;
  using Point = PointT<float>;
  using Circle = CircleT<float>;
  using Bounds = std::pair<Point, Point>;
  template <class Val>
  GridBase(const PointT<Val> & limits_, const uint32_t divisions__) :
      limits{limits_}, divisions_{divisions__}, data_(sqr(divisions_)) { }
  uint32_t divisions() const { return divisions_; }
  void kill() { Data temp; data_.swap(temp); }

 protected:
  static uint32_t get_divisions(const uint64_t n_items) {
    return max(1.0, sqrt(n_items / 1000));
  }
  const Cell & data(const uint64_t row, const uint64_t col) const {
    return data_[index(row, col)];
  }
  Cell & data(const uint64_t row, const uint64_t col) {
    return data_[index(row, col)];
  }
  uPoint cell(const Point & location) const {
    Point point(location * scale);
    if (point.x() < 0) point.x(0);
    if (point.y() < 0) point.y(0);
    if (point.x() >= divisions_) point.x(divisions_ - 1);
    if (point.y() >= divisions_) point.y(divisions_ - 1);
    return uPoint{point};
  }
  uPoint cell(const uint32_t index_) const {
    return uPoint(index_ / divisions_, index_ % divisions_);
  }
  Point center(const uPoint & cell_) const {
    return Point{(cell_.x() + 0.5f) / scale.x(),
          (cell_.y() + 0.5f) / scale.y()};
  }
  uint64_t index(const uint64_t row, const uint64_t col) const {
    return row * divisions_ + col;
  }
  Bounds bounds(const uPoint & cell_) const {
    return Bounds{Point{cell_.x() / scale.x(), cell_.y() / scale.y()},
          Point{(cell_.x() + 1) / scale.x(), (cell_.y() + 1) / scale.y()}};
  }
  Point limits;
  uint32_t divisions_;
  Point scale{divisions_ / limits.x(), divisions_ / limits.y()};
  Data data_;
};

// A class for grid-based densities
class GridDensity : public GridBase<double> {
 public:
  using Cell = double;
  using Cols = std::vector<Cell>;
  using Rows = std::vector<Cols>;
  template <class Val>
  GridDensity(const PointT<Val> & limits_, const uint64_t n_items,
              const double initial = 0.0) :
      GridBase{limits_, get_divisions(n_items)} {
    for (Cell & c : data_) c = initial;
  }
  // static uint32_t get_divisions(const uint64_t n_items) {
  //   return max(1.0, sqrt(n_items / 10000));
  // }
  void clear() { for (Cell & c : data_) c = 0; }
  void add(const Point & point) {
    const uPoint cell_{cell(point)};
    data(cell_.x(), cell_.y()) += 1;
  }
  void add(const Circle & circle) {
    const Point low_point{circle.center().x() - circle.radius(),
          circle.center().y() - circle.radius()};
    const Point high_point{circle.center().x() + circle.radius(),
          circle.center().y() + circle.radius()};
    const uPoint center_cell{cell(circle.center())};
    const uPoint low_cell{cell(low_point)};
    const uPoint high_cell{cell(high_point)};
    const float distance2{sqr(circle.radius())};
    for (uint32_t row{low_cell.x()}; row <= high_cell.x(); ++row) {
      for (uint32_t col{low_cell.y()}; col <= high_cell.y(); ++col) {
        const Point cell_center{center(uPoint{row, col})};
        if ((center_cell.x() == row && center_cell.y() == col) ||
          cell_center.distance2(circle.center()) <= distance2)
          data(row, col) += 1;
      }
    }
  }
  void finalize() {
    double total{0};
    for (Cell & c : data_) if (c <= 0) std::cerr << "Zero density" << std::endl;
    for (Cell & c : data_) total += c;
    for (Cell & c : data_) c /= total;
    for (uint64_t c{1}; c != data_.size(); ++c)
      data_[c] += data_[c - 1];
    for (uint64_t c{0}; false && c != data_.size(); ++c)
      std::cerr << c << " " << cell(c) << " " << data_[c] << std::endl;
    }
  template <class RNG>
  Bounds select(RNG & rng_) {
    const double value{d_unit_dist(rng_)};
    const uint32_t index_{static_cast<uint32_t>(
        upper_bound(data_.begin(), data_.end(), value) - data_.begin())};
    const uPoint cell_{cell(index_)};
    return bounds(cell_);
  }
  std::uniform_real_distribution<double> d_unit_dist{0, 1};
};

#if 0
// A class for location-based person lookup
class GridFinder : public GridBase<vector<PersonId>> {
 public:
  using ID = PersonId;
  using Cell = vector<ID>;
  template <class Val>
  GridFinder(const PointT<Val> & limits_, const ID n_people__) :
      GridBase{limits_, get_divisions(n_people__)}, n_people_{n_people__} { }
  PersonId n_people() const { return n_people_; }
  void clear() { for (Cell & c : data_) c.clear(); }
  void add(const Point & location, const ID id) {
    const uPoint cell_{cell(location)};
    data(cell_.x(), cell_.y()).push_back(id);
  }
  // TODO(paa) make more efficient versions using sampling
  // And pass function for getting pos
  Cell & items(const Circle & circle, const People & people) const {
    const Point low_point{circle.center().x() - circle.radius(),
          circle.center().y() - circle.radius()};
    const Point high_point{circle.center().x() + circle.radius(),
          circle.center().y() + circle.radius()};
    const uPoint low_cell{cell(low_point)};
    const uPoint high_cell{cell(high_point)};
    result.clear();
    for (ID row{low_cell.x()}; row <= high_cell.x(); ++row)
      for (ID col{low_cell.y()}; col <= high_cell.y(); ++col)
        for (const ID id : data(row, col))
          if (circle.contains(people[id].pos())) result.push_back(id);
    return result;
  }
  Cell & items(const Rectangle & rectangle, const People & people) const {
    const uPoint low_cell{cell(rectangle.ll())};
    const uPoint high_cell{cell(rectangle.ur())};
    result.clear();
    for (ID row{low_cell.x()}; row <= high_cell.x(); ++row)
      for (ID col{low_cell.y()}; col <= high_cell.y(); ++col)
        for (const ID id : data(row, col))
          if (rectangle.contains(people[id].pos())) result.push_back(id);
    return result;
  }
  static Bytes bytes(const PersonId n_people) {
    return sizeof(GridFinder) + sqr(get_divisions(n_people)) * sizeof(Cell) +
        3 * n_people * sizeof(PersonId) / 2;
  }

 private:
  static thread_local Cell result;
  ID n_people_;
};
thread_local GridFinder::Cell GridFinder::result{};
#endif

}  // namespace paa

#endif  // PAA_UTILITY_GRID_H_
