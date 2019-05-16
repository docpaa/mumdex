//
// animate.cpp
//
// Copyright 2019 Peter Andrews
//

#include <algorithm>
#include <iostream>
#include <fstream>
#include <functional>
#include <random>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#include <climits>

#include "colors.h"
#include "error.h"
#include "utility.h"
// #include "x11plot.h"

using std::cerr;
using std::cin;
using std::cout;
using std::endl;
using std::exception;
using std::function;
using std::mt19937_64;
using std::ostringstream;
using std::random_device;
using std::normal_distribution;
using std::ofstream;
using std::pair;
using std::string;
using std::to_string;
using std::vector;

using paa::Colors;
using paa::Error;
using paa::Progress;

#if 0
using paa::Click;
using paa::Color;
using paa::Geometry;
using paa::X11App;
using paa::X11Win;
class X11Colors : public X11Win {
 public:
  X11Colors(const X11Colors &) = delete;
  X11Colors & operator=(const X11Colors &) = delete;

  static constexpr int side{600};

  explicit X11Colors(
      X11App & app__,
      const size_t n_colors_ = 10,
      const bool order = false,
      const Geometry & geometry__ = Geometry{{side, side}, {0, 0}},
      const std::string title = "") :
      X11Win{app__,
        Geometry{{geometry__.width(), geometry__.height()},
      {geometry__.x_offset() +
            (order ? geometry__.width() + geometry__.width() / 20 : 0),
            geometry__.y_offset()}},
        true, title},
      color_names{},
      n_colors{n_colors_} {
        XSelectInput(display(), window,
                     StructureNotifyMask | ExposureMask |
                     ButtonPressMask | ButtonReleaseMask);

      // Shrink initial color name list if too long
      if (color_names.size() > n_colors) {
        color_names.resize(n_colors);
      }

      // Make initial colors
      for (const std::string & color_name : color_names) {
        colors.emplace_back(color_name);
      }

      // Expand initial color list if necessary
      const bool progress{false};
      const size_t initial_size{color_names.size()};
      if (color_names.size() != n_colors) {
        if (progress) std::cerr << "Size";
        while (color_names.size() != n_colors) {
          const unsigned int step{static_cast<unsigned int>(
              256 / pow(colors.size(), 1.0 / 3) / 2 + 1)};
          if (progress)
            std::cerr << " " << color_names.size() << " " << step << std::flush;
          colors.emplace_back(colors);
          color_names.push_back(colors.back().to_string());
        }
        if (progress) std::cerr << std::endl;
        // Move first made color to end
        if (color_names.size() != initial_size) {
          colors.push_back(colors[initial_size]);
          color_names.push_back(color_names[initial_size]);
          colors.erase(colors.begin() + initial_size);
          color_names.erase(color_names.begin() + initial_size);
        }
      }

      // X Colors and GCs
      Xcolors.resize(color_names.size());
      gcs.resize(color_names.size());
      for (unsigned int c{0}; c != color_names.size(); ++c) {
        std::string & color_name{color_names[c]};
        XColor & color{Xcolors[c]};
        if (!XAllocNamedColor(display(), app.colormap, color_name.c_str(),
                              &color, &color))
          throw Error("Could not get color") << color_name;

        // For colored arcs for points
        gcs[c] = create_gc(color.pixel, app.white, 2,
                           LineSolid, CapButt, JoinMiter);
      }
    }

  virtual void button_press(const XButtonEvent & event) {
    const Click click{event};
    if (click == 0) return;
    if (click > 1) return;
  }

  virtual void button_release(const XButtonEvent & event) {
    const Click click{event};
    if (click == 0) return;
  }

  virtual void draw() {
    clear_window();
    XFlush(display());
  }

  virtual ~X11Colors() {
    for (GC & gc_ : gcs) XFreeGC(display(), gc_);
  }

  std::vector<std::string> color_names{};

 private:
  std::vector<Color> colors{};
  std::vector<XColor> Xcolors{};
  std::vector<GC> gcs{};
  size_t n_colors;
};
#endif

template <class TYPE>
double sqr(const TYPE value) {
  return value * value;
}

using Pixel = pair<uint64_t, uint64_t>;
using Pixels = vector<Pixel>;

class Image {
 public:
  Image(
    const uint64_t & n_x__,
    const uint64_t & n_y__) :
  n_x_{n_x__},
  n_y_{n_y__},
  data(n_x_ * n_y_) { }

  uint64_t n_x() const { return n_x_; }
  Image & n_x(const uint64_t & n_x__) {
    n_x_ = n_x__;
    return *this;
  }
  uint64_t n_y() const { return n_y_; }
  Image & n_y(const uint64_t & n_y__) {
    n_y_ = n_y__;
    return *this;
  }

  unsigned char & operator()(const uint64_t x, const uint64_t y) {
    return data[x * n_y_ + y];
  }

  void clear(const Pixels & pixels) {
    for (const Pixel pixel : pixels) {
      (*this)(pixel.first, pixel.second) = 0;
    }
  }

  void save(const std::string & file_name) const {
    ofstream out{(file_name + ".xpm").c_str()};
    if (!out) throw Error{"Problem opening file in Image::save"} << file_name;
    out << "/* XPM */\nstatic char * XFACE[] = {\n/* <Values> */\n";
    out << "/* <width/cols> <height/rows> <colors> <char on pixel>*/\n";
    out << "\"" << n_y_ << " " << n_x_
        << " " << colors.size() + 1 << " 1\",\n/* <Colors> */\n";
    out << "\"" << static_cast<char>('#' + 0) << " c "
        << "#ffffff" << "\",\n";
    for (uint64_t c{0}; c != colors.size(); ++c)
      out << "\"" << static_cast<char>('#' + c + 1) << " c "
          << colors[c].hex() << "\",\n";
    out << "/* <Pixels> */\n";
    for (uint64_t x{0}; x != n_x_; ++x) {
      out << "\"";
      for (uint64_t y{0}; y != n_y_; ++y)
        out << static_cast<char>('#' + data[x * n_y_ + y]);
      if (x + 1 == n_x_) {
        out << "\"\n";
      } else {
        out << "\",\n";
      }
    }
    out << "};\n";
  }

 private:
  uint64_t n_x_;
  uint64_t n_y_;
  vector<unsigned char> data;
  Colors colors{50};
};

class MovingObject {
 public:
  MovingObject(
    const double & x__,
    const double & y__,
    const double & vx__,
    const double & vy__) :
  x_{x__},
  y_{y__},
  vx_{vx__},
  vy_{vy__} { }

  double x() const { return x_; }
  MovingObject & x(const double & x__) {
    x_ = x__;
    return *this;
  }
  double y() const { return y_; }
  MovingObject & y(const double & y__) {
    y_ = y__;
    return *this;
  }
  double vx() const { return vx_; }
  MovingObject & vx(const double & vx__) {
    vx_ = vx__;
    return *this;
  }
  double vy() const { return vy_; }
  MovingObject & vy(const double & vy__) {
    vy_ = vy__;
    return *this;
  }

  Pixels draw(Image & image) const {
    uint64_t rad{10};
    Pixels changed;
    for (uint64_t x{x_ > rad ? static_cast<uint64_t>(x_ - rad) : 0UL};
         x < image.n_x() && x <= x_ + rad; ++x) {
      for (uint64_t y{static_cast<uint64_t>(y_ > rad ? y_ - rad : 0)};
           y < image.n_y() && y <= y_ + rad; ++y) {
        const double distance{sqr(x - x_) + sqr(y - y_)};
        if (distance <= sqr(rad)) {
          image(x, y) = id;
          changed.emplace_back(x, y);
        }
      }
    }
    return changed;
  }

  void update_position(const vector<MovingObject> & objects,
                       const double max_x, const double max_y,
                       function<double()> & gen) {
    const double decay{0.9995};
    double a_x{0};
    double a_y{0};
    for (const MovingObject & object : objects) {
      if (id == object.id) continue;
      const double distance{
        std::max(0.00001, sqrt(sqr(object.x_ - x_) + sqr(object.y_ - y_)))};
      const double acceleration{std::min(1.0, 2000.0 / sqr(distance))};
      a_x += acceleration * (object.x_ - x_) / distance;
      a_y += acceleration * (object.y_ - y_) / distance;
    }
    // vx_ += 0.01 * (y_ - 500);
    vx_ += a_x;
    vy_ += a_y;
    x_ += vx_ + 0.01 * gen();
    y_ += vy_;
    if (x_ < 0) {
      x_ = 0;
      vx_ = -vx_;
    }
    if (y_ < 0) {
      y_ = 0;
      vy_ = -vy_;
    }
    if (x_ > max_x) {
      x_ = max_x;
      vx_ = -vx_;
    }
    if (y_ > max_y) {
      y_ = max_y;
      vy_ = -vy_;
    }
    vx_ *= decay;
    vy_ *= decay;
    // const double cdist{sqrt(sqr(x_ - 500) + sqr(y_ - 500))};
    // vx_ += cdist * 0.00001 * ((y_ - 500) > 0 ? 1 : -1);
    // vy_ += cdist * 0.00001 * ((x_ - 500) > 0 ? -1 : 1);
  }

 private:
  double x_;
  double y_;
  double vx_;
  double vy_;
  static uint64_t gid;
  uint64_t id{++gid};
};
uint64_t MovingObject::gid{0};

int main(int argc, char **) try {
  std::ios_base::sync_with_stdio(false);
  cin.tie(nullptr);

  if (--argc != 0) throw Error("usage: animate");

  random_device rd;
  mt19937_64 mersenne(rd());
  normal_distribution<double> dist{1, 0.5};
  function<double()> gen{bind(dist, std::ref(mersenne))};

  const uint64_t n_frames{10000};
  const uint64_t n_objects{20};
  const uint64_t n_x{1000};
  const uint64_t n_y{2000};

  vector<MovingObject> objects;
  for (uint64_t o{0}; o != n_objects; ++o) {
    const double x{500.0 + 500 * (gen() - 0.5)};
    const double y{500.0 + 500 * (gen() - 0.5)};
    const double center_dist{sqrt(sqr(x - 500) + sqr(y - 500))};
    const double vx{0.1 * (gen() - 0.5) +
          0.005 * center_dist * (y > 0 ? 1 : -1) * gen()};
    const double vy{0.1 * (gen() - 0.5) +
          0.005 * center_dist * (x > 0 ? -1 : 1) * gen()};
    objects.emplace_back(x, y, vx, vy);
  }
  Image image{n_x, n_y};

  Progress progress{n_frames, 1.0 / n_frames, "render"};
  for (uint64_t frame{0}; frame != n_frames; ++frame) {
    vector<Pixels> all_pixels;
    vector<MovingObject> objects_copy{objects};
    for (uint64_t o{0}; o != objects.size(); ++o) {
      MovingObject & object{objects[o]};
      object.update_position(objects_copy, 1.0 * n_x, 1.0 * n_y, ref(gen));
      if (0) std::cerr << "Frame " << frame
                       << " object " << o
                       << " x " << object.x()
                       << " y " << object.y()
                       << endl;
      Pixels pixels{object.draw(image)};
      all_pixels.push_back(move(pixels));
    }
    // progress bar
    for (uint64_t y{0}; y != image.n_y() * frame / n_frames; ++y)
      for (uint64_t x{0}; x != image.n_x() / 200; ++x)
        image(x, y) = 1;
    image.save("image" + to_string(frame));
    progress();
    for (Pixels & pixels : all_pixels) image.clear(pixels);
  }

  ostringstream convert;
  convert << "echo 'convert -delay 0.2 -loop 0 image{1.."
          << n_frames - 1 << "}.xpm animate.gif' | bash";
  if (system(convert.str().c_str()) != 0) cerr << "Bad sys call" << endl;

  return 0;
} catch (Error & e) {
  cerr << "paa::Error:" << endl;
  cerr << e.what() << endl;
  return 1;
} catch (exception & e) {
  cerr << e.what() << endl;
  return 1;
} catch (...) {
  cerr << "Some exception was caught." << endl;
  return 1;
}
