//
// union_find.h
//
// class that implements union find data structure
//
// Copyright 2020 Peter Andrews @ CSHL
//

#ifndef PAA_UTILITY_UNION_FIND_H_
#define PAA_UTILITY_UNION_FIND_H_

#include <algorithm>
#include <unordered_map>
#include <vector>

namespace paa {

template <class UInt>
class UnionFind {
 public:
  explicit UnionFind(const UInt n) : data(n), sizes(n, 1) {
    iota(data.begin(), data.end(), 0ul);
  }
  UInt size() const { return data.size(); }
  UInt size(const UInt i) { return sizes[find(i)]; }
  UInt find(UInt i) {
    UInt j{i};
    while (data[i] != i)
      i = data[i];
    while (data[j] != j) {
      const UInt k{data[j]};
      data[j] = i;
      j = k;
    }
    return i;
  }
  void unify(UInt i, UInt j) {
    i = find(i);
    j = find(j);
    if (i == j) return;
    if (sizes[i] < sizes[j]) {
      data[i] = j;
      sizes[j] += sizes[i];
    } else {
      data[j] = i;
      sizes[i] += sizes[j];
    }
  }

 private:
  std::vector<UInt> data;
  std::vector<UInt> sizes;
};

template <class UInt, class Val>
class UnionFindVal {
 public:
  template <class Iter>
  UnionFindVal(const Iter begin, const Iter end) :
      map{make_map(begin, end)}, uf(map.size()) {}
  UInt size() const { return uf.size(); }
  UInt size(const Val & val) { return uf.size(map.at(val)); }
  UInt find(const Val & val) { return uf.find(map.at(val)); }
  void unify(const Val & v1, const Val & v2) {
    return uf.unify(map.at(v1), map.at(v2));
  }

 private:
  using Map = std::unordered_map<Val, UInt>;
  template <class Iter>
  static Map make_map(const Iter begin, const Iter end) {
    Map map;
    for (Iter i{begin}; i != end; ++i) map.emplace(*i, map.size());
    return map;
  }
  std::unordered_map<Val, UInt> map;
  UnionFind<UInt> uf;
};

}  // namespace paa

#endif  // PAA_UTILITY_UNION_FIND_H_
