//
// kdtree.h
//
// class that implements a kd tree
//
// Copyright 2020 Peter Andrews @ CSHL
//

#ifndef PAA_UTILITY_KDTREE_H_
#define PAA_UTILITY_KDTREE_H_

#include <algorithm>
#include <limits>
#include <queue>
#include <utility>
#include <vector>

#include "utility.h"

namespace paa {

// An n-dimensional position
template <class Type, uint64_t n_dimensions>
class PointD {
 public:
  PointD() = default;
  PointD(const PointD &) = default;
  template <class POINT>
  PointD(const POINT & point, bool, bool) :
      data{point.x(), point.y()} {}
  Type & operator[](const uint64_t dimension) { return data[dimension]; }
  Type operator[](const uint64_t dimension) const { return data[dimension]; }
  template <class POINT>
  double distance2(const POINT & other) const {
    double result{0};
    for (uint64_t dimension{0}; dimension != n_dimensions; ++dimension)
      result += sqr(static_cast<double>(data[dimension]) -
                    static_cast<double>(other[dimension]));
    return result;
  }
  template <class POINT>
  double distance(const POINT & other) const {
    return sqrt(distance2(other));
  }
  bool operator!=(const PointD & other) const {
    for (uint64_t dimension{0}; dimension != n_dimensions; ++dimension)
      if (data[dimension] < other[dimension] ||
          data[dimension] > other[dimension]) return true;
    return false;
  }

 private:
  Type data[n_dimensions];  // NOLINT
};

template <class Id, class Type, uint64_t n_dimensions>
class IdPointD : public PointD<Type, n_dimensions> {
 public:
  using Point = PointD<Type, n_dimensions>;
  IdPointD() : Point{}, id{bad} {}
  IdPointD(const Id id_, const Point & point) : Point{point}, id{id_} {}
  Id id;
  static constexpr Id bad{std::numeric_limits<Id>::max()};
};

// A simple k-d tree implementation
template <class Id, class Type, uint64_t n_dimensions>
class KDTreeT {
  using Point = IdPointD<Id, Type, n_dimensions>;
  using Points = std::vector<Point>;
  struct Node {
    explicit Node(const Point & point_,
                  const Node * left_, const Node * right_) :
        point{point_}, left{left_}, right{right_} {}
    Node(const Node &) = delete;
    Node & operator=(const Node &) = delete;
    ~Node() { for (const Node * node : {left, right}) if (node) delete node; }
    const Point point;
    const Node * const left;
    const Node * const right;
  };
  struct IdDist {
    IdDist() = default;
    IdDist(const Id & id_, const double distance2_) :
      id{id_}, distance2{distance2_} {}
    bool operator<(const IdDist & rhs) const {
      return distance2 < rhs.distance2;
    }
    Id id{std::numeric_limits<Id>::max()};
    double distance2{big};
  };

 public:
  static constexpr Id bad{Point::bad};
  static constexpr double big{std::numeric_limits<double>::max()};
  template <class POINTS>
  explicit KDTreeT(const POINTS & points) : root{build(points)} {}
  KDTreeT(const KDTreeT &) = delete;
  KDTreeT & operator=(const KDTreeT &) = delete;
  ~KDTreeT() { if (root) delete root; }
  template <class POINT>
  Id find_closest(const POINT & point) const {
    return find_closest(point, [](Id) { return true; }, root, 0).id;
  }
  template <class POINT, class FUN>
  Id find_closest(const POINT & point, FUN && fun) const {
    return find_closest(point, std::forward<FUN>(fun), root, 0).id;
  }
  using Ids = std::vector<Id>;
  using Queue = std::priority_queue<IdDist>;
  template <class POINT>
  Ids find_n_closest(const POINT & point, const uint64_t n) const {
    return find_n_closest(point, n, [](Id) { return true; });
  }
  template <class POINT, class FUN>
  Ids find_n_closest(const POINT & point, const uint64_t n, FUN && fun) const {
    Queue closest;
    find_n_closest(closest, point, n, std::forward<FUN>(fun), root, 0);
    Ids result;
    if (closest.size()) {
      result.reserve(closest.size());
      while (closest.size()) {
        result.push_back(closest.top().id);
        closest.pop();
      }
    }
    return result;
  }

 private:
  template <class POINT, class FUN>
  static IdDist find_closest(const POINT & point,
                             FUN && fun,
                             const Node * const node,
                             const uint64_t depth) {
    if (!node) return IdDist{};
    const Point & current{node->point};
    IdDist closest{};
    if (fun(current.id)) closest = IdDist{current.id, current.distance2(point)};
    const uint64_t dimension{depth % n_dimensions};
    const bool to_left{point[dimension] < current[dimension]};
    const Node * const inside_half{to_left ? node->left : node->right};
    const IdDist inside_best{find_closest(point, std::forward<FUN>(fun),
                                          inside_half, depth + 1)};
    if (inside_best.distance2 < closest.distance2) closest = inside_best;
    if (closest.distance2 >=
        sqr(static_cast<double>(point[dimension]) -
            static_cast<double>(current[dimension]))) {
      const Node * const overlap_half{to_left ? node->right : node->left};
      const IdDist overlap_best{
        find_closest(point, std::forward<FUN>(fun), overlap_half, depth + 1)};
      if (overlap_best.distance2 < closest.distance2) closest = overlap_best;
    }
    return closest;
  }
  template <class POINT, class FUN>
  static void find_n_closest(Queue & queue,
                             const POINT & point,
                             const uint64_t n,
                             FUN && fun,
                             const Node * const node,
                             const uint64_t depth) {
    if (!node) return;
    const Point & current{node->point};
    IdDist closest{};
    if (fun(current.id)) closest = IdDist{current.id, current.distance2(point)};
    if (queue.size() < n) {
      if (closest.distance2 < big) queue.push(closest);
    } else if (queue.top().distance2 > closest.distance2) {
      queue.pop();
      queue.push(closest);
    }
    const uint64_t dimension{depth % n_dimensions};
    const bool to_left{point[dimension] < current[dimension]};
    const Node * const inside_half{to_left ? node->left : node->right};
    find_n_closest(queue, point, n, std::forward<FUN>(fun),
                   inside_half, depth + 1);
    if (queue.size() < n || queue.top().distance2 >=
        sqr(static_cast<double>(point[dimension]) -
            static_cast<double>(current[dimension]))) {
      const Node * const overlap_half{to_left ? node->right : node->left};
      find_n_closest(queue, point, n, std::forward<FUN>(fun),
                     overlap_half, depth + 1);
    }
  }
  const Node * const root;
  // Build the tree from initial points
  template <class POINTS>
  static Node * build(POINTS points_) {
    Points points;
    points.reserve(points_.size());
    for (Id i{0}; i != points_.size(); ++i)
      points.emplace_back(i, points_[i]);
    return build(points, 0);
  }
  static Node * build(Points points, uint64_t depth) {
    if (points.empty()) return nullptr;

    // Get dimension to act upon
    const uint64_t dimension{depth % n_dimensions};

    // Sort to find median
    sort(points.begin(), points.end(),
         [dimension](const Point & lhs, const Point & rhs) {
           return lhs[dimension] < rhs[dimension];
         });
    const uint64_t median{points.size() / 2};
    const Points right_points{points.begin() + median + 1, points.end()};
    points.resize(median);
    ++depth;
    return new Node{
      points[median], build(points, depth), build(right_points, depth)};
  }
};

}  // namespace paa

#endif  // PAA_UTILITY_KDTREE_H_
