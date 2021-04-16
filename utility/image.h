//
// image.h
//
// Class for bitmap images and movies
//
// Copyright 2020 Peter Andrews @ CSHL
//

#ifndef PAA_IMAGE_H
#define PAA_IMAGE_H

#include <string>
#include <utility>
#include <vector>

#include "error.h"
#include "files.h"
#include "utility.h"

namespace paa {

using Coordinate = uint32_t;
using PixelIndex = uint32_t;
using ColorInt = uint8_t;
using ColorInts = std::vector<ColorInt>;
using ImgPixel = std::pair<PixelIndex, PixelIndex>;
using ImgPixels = std::vector<ImgPixel>;

// Point class
template <class Type>
class PointT {
 public:
  constexpr PointT() = default;
  constexpr PointT(const Type x__, const Type y__) : x_{x__}, y_{y__} { }
  template <class Type2>
  explicit constexpr PointT(const PointT<Type2> & other) :
  x_(other.x()), y_(other.y()) { }

  constexpr Type operator[](const bool y__) const { return y__ ? y_ : x_; }
  Type & operator[](const bool y__) { return y__ ? y_ : x_; }
  Type distance2(const PointT other) const {
    return Type(sqr(x_ - other.x_) + sqr(y_ - other.y_));
  }
  Type distance(const PointT other) const {
    return Type(sqrt(sqr(x_ - other.x_) + sqr(y_ - other.y_)));
  }
  template <class Type2>
  constexpr PointT operator*(const Type2 val) const {
    return PointT(x_ * val, y_ * val);
  }
  constexpr PointT operator/(const Type val) const {
    return PointT{x_ / val, y_ / val};
  }
  constexpr PointT operator*(const PointT rhs) const {
    return PointT{x_ * rhs.x_, y_ * rhs.y_};
  }
  constexpr PointT operator/(const PointT rhs) const {
    return PointT{x_ / rhs.x_, y_ / rhs.y_};
  }
  void x(const Type & val) { x_ = val; }
  void y(const Type & val) { y_ = val; }
  constexpr Type x() const { return x_; }
  constexpr Type y() const { return y_; }

 private:
  Type x_{0};
  Type y_{0};
};
using iPoint = PointT<int>;
using uPoint = PointT<unsigned int>;
using dPoint = PointT<double>;
using Point = dPoint;

template <class Type>
std::ostream & operator<<(std::ostream & out, const PointT<Type> point) {
  out << point.x() << " x " << point.y();
  return out;
}

// A rectangle
template <class Val>
class RectangleT {
 public:
  RectangleT(const PointT<Val> & ll__, const PointT<Val> & ur__) :
      ll_{ll__}, ur_{ur__} { }
  bool contains(const PointT<Val> & point) const {
    return (point.x() >= ll_.x() && point.x() < ur_.x() &&
            point.y() >= ll_.y() && point.y() < ur_.y());
  }
  PointT<Val> ll() const { return ll_; }
  PointT<Val> lr() const { return PointT<Val>{ur_.x(), ll_.y()}; }
  PointT<Val> ur() const { return ur_; }
  PointT<Val> ul() const { return PointT<Val>{ll_.x(), ur_.y()}; }
  Val width() const { return ur_.x() - ll_.x(); }
  Val height() const { return ur_.y() - ll_.y(); }
  double area() const {
    return 1.0 * static_cast<double>(width()) * static_cast<double>(height());
  }
  Val min_x() const { return ll_.x(); }
  Val max_x() const { return ur_.x(); }
  Val min_y() const { return ll_.y(); }
  Val max_y() const { return ur_.y(); }

 private:
  PointT<Val> ll_;
  PointT<Val> ur_;
};
using Rectangle = RectangleT<double>;

template <class Type>
class CircleT {
 public:
  using POINT = PointT<Type>;
  CircleT(const POINT & center__, const Type radius__) :
      center_{center__}, radius_{radius__} {}

  POINT center() const { return center_; }
  Type radius() const { return radius_; }
  double area() const { return pi * sqr(radius_); }
  Type x() const { return center_.x(); }
  Type y() const { return center_.y(); }
  bool contains(const POINT & point) const {
    return center_.distance(point) < radius_;
  }

 private:
  static constexpr Type pi = PI;
  POINT center_;
  Type radius_;
};
template <class Type>
std::ostream & operator<<(std::ostream & out, const CircleT<Type> circle) {
  out << circle.center() << " * " << circle.radius();
  return out;
}

// Movie image
class Image {
 public:
  using HexColors = std::vector<std::string>;
  constexpr static ColorInt reserved_colors{2};
  Image(const Coordinate n_x_, const Coordinate n_y_,
        const HexColors & hex_colors_) :
      hex_colors{hex_colors_}, size{n_x_, n_y_}, data(size.x() * size.y()) { }

  Coordinate n_x() const { return size.x(); }
  Coordinate n_y() const { return size.y(); }
  ColorInt operator()(const Coordinate x, const Coordinate y) const {
    return data[x * n_y() + y];
  }
  void set(const Coordinate x, const Coordinate y, ColorInt value) {
    data[x * n_y() + y] = value;
  }
  void set(const PointT<Coordinate> pos, ColorInt value) {
    data[pos.x() * n_y() + pos.y()] = value;
  }
  void clear(const ImgPixels & pixels) {
    for (const ImgPixel & pixel : pixels) set(pixel.first, pixel.second, 0);
  }
  void clear() { data.assign(data.size(), 0); }

  void save(const std::string & file_name) const {
    std::ofstream out{(file_name + ".xpm").c_str()};
    if (!out) throw Error{"Problem opening file in Image::save"} << file_name;
    out << "/* XPM */\nstatic char * XFACE[] = {\n/* <Values> */\n";
    out << "/* <width/cols> <height/rows> <colors> <char on pixel> */\n";
    out << "\"" << n_x() << " " << n_y() << " "
        << hex_colors.size() + reserved_colors
        << " 1\",\n/* <ColorInts> */\n";
    out << "\"" << static_cast<char>('#' + 0) << " c "
        << "#ffffff" << "\",\n";
    out << "\"" << static_cast<char>('#' + 1) << " c "
        << "#000000" << "\",\n";
    for (ColorInt c{0}; c != hex_colors.size(); ++c)
      out << "\"" << static_cast<char>('#' + c + reserved_colors) << " c "
          << hex_colors[c] << "\",\n";
    out << "/* <Pixels> */\n";
    for (Coordinate yy{0}; yy != n_y(); ++yy) {
      const Coordinate y{n_y() - yy - 1};
      out << "\"";
      for (Coordinate x{0}; x != n_x(); ++x)
        out << static_cast<char>('#' + data[x * n_y() + y]);
      if (y) {
        out << "\",\n";
      } else {
        out << "\"\n";
      }
    }
    out << "};\n";
  }
  static uint64_t bytes(const uint64_t n_pixels) {
    return sizeof(Image) + n_pixels * sizeof(ColorInt);
  }

 private:
  HexColors hex_colors{};
  uPoint size;
  ColorInts data;
};

}  // namespace paa

#endif
