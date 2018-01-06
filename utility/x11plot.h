//
// x11plot.h
//
// Plot on X11
//
// Copyright 2016 Peter Andrews @ CSHL
//

#ifndef PAA_X11PLOT_H
#define PAA_X11PLOT_H

#include <X11/Xlib.h>
#include <X11/Xutil.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <deque>
#include <iomanip>
#include <list>
#include <functional>
#include <fstream>
#include <map>
#include <memory>
#include <numeric>
#include <set>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include "error.h"
#include "files.h"
#include "plot.h"
#include "threads.h"
#include "utility.h"

namespace paa {

// Point class
template <class Type>
struct PointT {
  constexpr PointT() = default;
  constexpr PointT(const Type x_, const Type y_) : x{x_}, y{y_} { }
  template <class Event>
  PointT(const Event & event) {  // NOLINT
    x = event.x;
    y = event.y;
  }

  template <class Event>
  PointT & operator=(const Event & event) {
    x = event.x;
    y = event.y;
    return *this;
  }
  Type operator[](const bool y_) const { return y_ ? y : x; }
  Type & operator[](const bool y_) { return y_ ? y : x; }
  bool operator==(const PointT rhs) const { return x == rhs.x && y == rhs.y; }
  bool operator!=(const PointT rhs) const { return x != rhs.x || y != rhs.y; }
  double distance(const PointT other) const {
    return sqrt((x - other.x) * (x - other.x) + (y - other.y) * (y - other.y));
  }
  double distance(const Type x_, const Type y_) const {
    return distance(PointT{x_, y_});
  }

  Type x{0};
  Type y{0};
};
using Point = PointT<int>;  // Integer point
using dPoint = PointT<double>;  // Double point

class X11Font {
 public:
  X11Font(Display * display_, const unsigned int point_size,
          const std::string font_name = "helvetica",
          const bool fallback = false) : display{display_} {
    const std::string font_spec{std::string("-*-") + font_name + "-*-r-*-*-" +
          std::to_string(point_size) + "-*-*-*-*-*-*-*"};
    font = XLoadQueryFont(display, font_spec.c_str());
    if (fallback && !font) XLoadQueryFont(display, "fixed");
  }
  X11Font(X11Font &) = delete;
  X11Font & operator=(const X11Font &) = delete;
  X11Font & operator=(X11Font &&) = delete;
  X11Font(X11Font && other) : font{other.font} { other.font = nullptr; }

  operator bool() const { return font != nullptr; }
  Font id() const { return font->fid; }
  int width() const {
    return font->max_bounds.rbearing - font->max_bounds.lbearing;
  }
  int height() const {
    return font->max_bounds.ascent + font->max_bounds.descent;
  }
  int origin_height() const { return -font->max_bounds.descent; }
  int string_width(const std::string & text) const {
    // Fix this?
    return XTextWidth(font, text.c_str(),
                      static_cast<unsigned int>(text.size()));
  }
  int centered_y(const int y) const {
    return y + (font->max_bounds.ascent - font->max_bounds.descent) / 2;
  }
  int below_y(const int y) const {
    return y + font->max_bounds.ascent;
  }
  int centered_x(const std::string & text, const int x) const {
    return x - (string_width(text) + 1) / 2 -
        font->per_char[static_cast<int>(text[0])].lbearing + 1;
  }

  ~X11Font() { }  // segfault on if (font) XUnloadFont(display, id()); }

  Display * display{nullptr};
  XFontStruct * font{nullptr};
};

class X11Fonts {
 public:
  explicit X11Fonts(Display * display, const std::string & name = "helvetica") {
    for (unsigned int s{3}; s <= max_font_size; ++s) {
      X11Font font{display, s, name, s == max_font_size};
      if (font) {
        sizes.push_back(s);
        lookup[s] = static_cast<unsigned int>(fonts.size());
        fonts.push_back(std::move(font));
      }
    }
    if (fonts.empty()) throw Error("No fonts loaded");
  }
  X11Font * size(const unsigned int points) const {
    return &fonts[lookup.at(points)];  // can throw exception
  }
  X11Font * at_least(const unsigned int points) const {
    auto found = lower_bound(sizes.begin(), sizes.end(), points);
    if (found == sizes.end()) --found;
    return &fonts[found - sizes.begin()];
  }
  X11Font * at_most(const unsigned int points) const {
    auto found = upper_bound(sizes.begin(), sizes.end(), points);
    if (found == sizes.begin()) ++found;
    return &fonts[found - 1 - sizes.begin()];
  }
  X11Font * fits(const std::string text,
                 const int width, const int height) const {
    for (unsigned int f{0}; f != fonts.size(); ++f) {
      X11Font * font{&fonts[fonts.size() - f - 1]};
      if (font->height() > height) continue;
      if (font->string_width(text) < width) return font;
    }
    return &fonts[0];
  }
  void clear() { fonts.clear(); }

  static constexpr unsigned int max_font_size{60};
  std::map<unsigned int, unsigned int> lookup{};
  std::vector<unsigned int> sizes{};
  mutable std::vector<X11Font> fonts{};
};

using iBounds = std::vector<std::vector<int>>;
inline bool operator!=(const iBounds & lhs, const iBounds & rhs) {
  if (lhs.size() != rhs.size()) return true;
  for (unsigned int y{0}; y != lhs.size(); ++y) {
    if (lhs[y].size() != rhs[y].size()) return true;
    for (unsigned int d{0}; d != lhs[y].size(); ++d) {
      if (lhs[y][d] != rhs[y][d]) return true;
    }
  }
  return false;
}

std::string hex(const XColor & color) {
  std::string result;
  static std::string chars{"0123456789abcdefxxx"};
  for (const unsigned int component : {color.red, color.green, color.blue}) {
    result += chars[((component / 256) % 256) / 16];
    result += chars[(component % 256) / 16];
  }
  return result;
}

// Due to X11 Macros using hidden casts
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"

// Window class in an app
template <class X11App>
class X11WindowT {
 public:
  using X11Win = X11WindowT<X11App>;
  static constexpr double window_scale{1.65};
  static constexpr int radio_width{1};

  // Factory function
  static X11WindowT & create(X11App & app) {
    return app.add(std::make_unique<X11Win>(app));
  }

  X11WindowT(X11App & app_,
             const unsigned int width__ = default_window_width * window_scale,
             const unsigned int height__ = default_window_height * window_scale,
             const int x__ = 100, const int y__ = 100, const bool map_ = true) :
      app{app_}, size_{width__, height__},
    window{XCreateSimpleWindow(display(), DefaultRootWindow(display()),
                               x__, y__, width(), height(),
                               0, app.white, app.white)}
  {
    // Window properties
    XSelectInput(display(), window, StructureNotifyMask | ExposureMask);
    XSetWMProtocols(display(), window, app.wmDeleteMessage(), 1);
    XSetWindowBackgroundPixmap(display(), window, None);
    if (true) {
      if (0) std::cerr << "Resizing to "
                       << width() << "x" << height() << " "
                       << x__ << " " << y__ << std::endl;
      XSizeHints hints;
      hints.flags = PPosition | PSize;
      hints.x = x__;
      hints.y = y__;
      hints.width = width();
      hints.height = height();
      XSetNormalHints(display(), window, &hints);
    }

    const Font font{app.fonts.at_most(30)->id()};

    // Colors
    XSync(display(), False);
    colormap = DefaultColormap(display(), app.screen);

    // Graphics contexts
    gc = create_gc(app.black, app.white);
    XSetFont(display(), gc, font);
    fill_gc = create_gc(app.white, app.black);
    XSetFont(display(), fill_gc, font);
    radio_gc = create_gc(app.black, app.white, radio_width,
                         LineSolid, CapButt, JoinMiter);

    // Maximum request size for draw events (for points; divide for arcs, etc)
    max_request = XMaxRequestSize(display()) - 3;

    if (map_) XMapWindow(display(), window);
  }

  X11WindowT(const X11WindowT &) = delete;
  X11WindowT & operator=(const X11WindowT &) = delete;

 public:
  virtual ~X11WindowT() {
    for (GC g : {gc, fill_gc}) XFreeGC(display(), g);
    XDestroyWindow(display(), window);
  }

  // Window information
  Display * display() const { return app.display; }
  int width() const { return size_[0]; }
  int height() const { return size_[1]; }
  int extent(const bool y) { return size_[y]; }
  virtual bool slow() const { return false; }

  GC create_gc(const uint64_t foreground, const uint64_t background,
               const unsigned int line_width = 1,
               const unsigned int line_type = LineSolid,
               const unsigned int butt_type = CapButt,
               const unsigned int join_type = JoinMiter) const {
    const GC result{XCreateGC(display(), window, 0, nullptr)};
    XSetForeground(display(), result, foreground);
    XSetBackground(display(), result, background);
    XSetLineAttributes(display(), result, line_width,
                       line_type, butt_type, join_type);
    return result;
  }

  // Event actions
  virtual void configure(const XConfigureEvent & event) {
    if (size_[0] != static_cast<unsigned int>(event.width) ||
        size_[1] != static_cast<unsigned int>(event.height)) {
      just_configured = true;
      size_[0] = event.width;
      size_[1] = event.height;
      prepare_draw();
    }
  }
  virtual void expose(const XExposeEvent & event) {
    if (event.count == 0) return prepare_draw();
  }
  virtual void enter(const XCrossingEvent &) { }
  virtual void key(const XKeyEvent &) { }
  virtual void button_press(const XButtonEvent &) { }
  virtual void motion(const XMotionEvent &) { }
  virtual void button_release(const XButtonEvent &) { }
  virtual void leave(const XCrossingEvent &) { }
  virtual void client_message(const XClientMessageEvent &) { }

  void set_bounds(const bool y, const int l, const int h) {
    bounds[y][2] = (bounds[y][1] = h) - (bounds[y][0] = l);
  }
  void set_bounds(const int xl, const int xh, const int yl, const int yh) {
    const bool valid_initial_bounds{bounds.size() > 0};
    if (!valid_initial_bounds) bounds.assign(2, std::vector<int>(3));
    const iBounds last_bounds = bounds;
    set_bounds(0, xl, xh);
    set_bounds(1, yl, yh);
    if (bounds != last_bounds) {
      if (valid_initial_bounds) XFreePixmap(display(), pixmap);
      pixmap = XCreatePixmap(display(), window, width(), height(), app.depth);
    }
  }
  template <class POINT>
  bool in_bounds(const POINT & point) const {
    return in_bounds(point.x, point.y);
  }
  bool in_bounds(const int x, const int y) const {
    return x > bounds[0][0] && x < bounds[0][1] &&
        y > bounds[1][0] && y < bounds[1][1];
  }

  virtual void save_image(const std::string & base_name,
                          std::function<void()> call_back = [] () {}) {
    static const std::string pdf_name{get_next_file(base_name, "pdf")};
    const std::string image_name{get_next_file(base_name, "xpm")};
    // Not sure why I need to subtract 1 from width and height on Mac OS
    save_image(image_name, window, 0, 0, width() - 1, height() - 1, call_back);
    image_names.push_back(image_name);
    std::ostringstream pdf_command;
    pdf_command << "convert -quality 100 -density "
                << app.pixels_per_inch(0) << "x" << app.pixels_per_inch(1);
    for (const std::string & name : image_names) {
      pdf_command << " " << name;
    }
    pdf_command << " " << pdf_name;
    // std::cerr << pdf_command.str() << std::endl;
    std::cerr << "Saved " << image_names.size() << " image"
              << (image_names.size() == 1 ? "" : "s so far") << " in pdf file "
              << pdf_name << std::endl;
    if (system(pdf_command.str().c_str()) == -1) {
      std::cerr << "Problem creating pdf file" << std::endl;
    }
  }
  void save_image(const std::string & file_name, Drawable d,
                  const int xp, const int yp,
                  const unsigned int w, const unsigned int h,
                  std::function<void()> call_back = [] () {}) {
    if (0) std::cerr << xp << " " << yp << " "
                     << w << " " << h << " " << std::endl;
    XImage * image{XGetImage(display(), d, xp, yp, w, h, -1UL, XYPixmap)};
    if (!image) throw Error("Could not get image");
    call_back();
    std::map<uint64_t, char> colors;
    XColor color;
    std::ostringstream color_string;
    std::ostringstream image_string;
    for (unsigned int y{0}; y != h; ++y) {
      image_string << '"';
      for (unsigned int x{0}; x != w; ++x) {
        color.pixel = XGetPixel(image, x, y);
        auto inserted = colors.emplace(color.pixel, 'a' + colors.size());
        if (inserted.second == true) {
          XQueryColor(display(), colormap, &color);
          color_string << '"' << inserted.first->second << " "
                       << "c #" << hex(color) << '"' << ",\n";
        }
        image_string << inserted.first->second;
      }
      image_string << '"' << (y + 1 == h ? "" : ",") << "\n";
    }
    std::ofstream file{file_name.c_str()};
    if (!file) throw Error("Problem opening file") << file_name;
    file << "/* XPM */\n"
         << "static char * XFACE[] = {\n"
         << "/* <Values> */\n"
         << "/* <width/cols> <height/rows> <colors> <char on pixel>*/\n"
         << '"' << w << " " << h << " " << colors.size()
         << " 1" << '"' << ",\n"
         << "/* <Colors> */\n"
         << color_string.str()
         << "/* <Pixels> */\n"
         << image_string.str()
         << "};\n";
    std::cerr << "Saved image to " << file_name << std::endl;
    XDestroyImage(image);
  }

  // Data preparation and drawing functions (just an empty window here)
  virtual void prepare() { }
  virtual void draw() {
    XFillRectangle(display(), window, fill_gc, 0, 0, width(), height());
  }
  void prepare_draw() { prepare(); draw(); }

  X11App & app;
  unsigned int size_[2]{0, 0};
  iBounds bounds{};
  Window window{};
  Colormap colormap{};
  Pixmap pixmap{};
  std::vector<std::string> image_names{};
  bool inside{true};

  GC gc{}, fill_gc{}, radio_gc{};
  uint64_t max_request{};
  bool destroyed{false};
  mutable bool just_configured{true};
  static constexpr unsigned int default_window_width{default_doc_width};
  static constexpr unsigned int default_window_height{default_doc_height};
};

inline Display * open_default_display() {
  Display * display = XOpenDisplay(nullptr);
  if (display == nullptr) {
    throw Error("Could not open X display - "
                "is X windowing enabled in your terminal?");
  }
  return display;
}

// Application to display several windows of different types
class X11App {
 public:
  using X11Win = X11WindowT<X11App>;
  using WinPtr = std::unique_ptr<X11Win>;
  using WinList = std::list<WinPtr>;
  using WinIter = WinList::iterator;

  X11App() : display{open_default_display()},
    screen{DefaultScreen(display)}, depth{DefaultDepth(display, screen)},
    fonts{display} {
    display_size[0] = DisplayWidth(display, screen);
    display_size[1] = DisplayHeight(display, screen);
    display_mm[0] = DisplayWidthMM(display, screen);
    display_mm[1] = DisplayHeightMM(display, screen);
    black = BlackPixel(display, screen);
    white = WhitePixel(display, screen);
    if (!display) throw Error("Problem opening X11 display");
    wmDeleteMessage_ = {XInternAtom(display, "WM_DELETE_WINDOW", False)};
  }

  X11App(const X11App &) = delete;
  X11App & operator=(const X11App &) = delete;

  ~X11App() {
    fonts.clear();
    windows.clear();
    XCloseDisplay(display);
  }

  // Add a window to app
  X11Win & add(WinPtr && ptr) {
    WinIter win = windows.emplace(windows.end(), std::move(ptr));
    window_lookups[(*win)->window] = win;
    return **win;
  }

  unsigned int pixels_per_inch(const bool y) const {
    return 25.4 * display_size[y] / display_mm[y];
  }
  unsigned int pixels_per_inch() const {
    return std::max(pixels_per_inch(0), pixels_per_inch(1));
  }

  // Run the application
  void run() {
    // Only keep latest configure requests to avoid too many redraws
    std::map<Window, XConfigureEvent> configures;

    // Event loop
    while (windows.size()) {
      // Configure only when no events are pending
      if (!XPending(display)) {
        for (auto pair : configures) {
          (*window_lookups[pair.first])->configure(pair.second);
        }
        configures.clear();
      }

      // Get event
      XNextEvent(display, &event);

      // Process event based on type
      switch (event.type) {
        case ConfigureNotify:
          {
            const Window win{event.xconfigure.window};
            WinIter win_iter{window_lookups[win]};
            if ((*win_iter)->slow()) {
              configures[win] = event.xconfigure;
            } else {
              (*win_iter)->configure(event.xconfigure);
            }
          }
          break;

        case MapNotify:
          break;

        case VisibilityNotify:
          break;

        case Expose:
          (*window_lookups[event.xcrossing.window])->expose(event.xexpose);
          break;

        case EnterNotify:
          (*window_lookups[event.xcrossing.window])->enter(event.xcrossing);
          break;

        case KeyPress:
          (*window_lookups[event.xkey.window])->key(event.xkey);
          break;

        case ButtonPress:
          (*window_lookups[event.xbutton.window])->button_press(event.xbutton);
          break;

        case MotionNotify:
          (*window_lookups[event.xmotion.window])->motion(event.xmotion);
          break;

        case ButtonRelease:
          (*window_lookups[event.xbutton.window])->
              button_release(event.xbutton);
          break;

        case LeaveNotify:
          (*window_lookups[event.xcrossing.window])->leave(event.xcrossing);
          break;

        case ClientMessage:
          if (static_cast<uint64_t>(event.xclient.data.l[0]) ==
              *wmDeleteMessage()) {
            const Window window{event.xclient.window};
            close_window(window);
          } else {
            (*window_lookups[event.xmotion.window])->
                client_message(event.xclient);
          }
          break;

        case DestroyNotify:
          break;

        default:
          break;
      }
    }
  }

  void close_window(const Window window) {
    XSelectInput(display, window, 0);
    while (XCheckWindowEvent(display, window, -1UL, &event)) { }
    windows.erase(window_lookups[window]);
    window_lookups.erase(window);
  }

  Atom * wmDeleteMessage() const {  // Allows proper pass of window kill
    return const_cast<Atom *>(&wmDeleteMessage_);
  }

  Display * display{};
  int screen{};
  int depth{};
  X11Fonts fonts;
  unsigned int display_size[2]{0, 0};
  unsigned int display_mm[2]{0, 0};
  uint64_t black{}, white{};
  XEvent event{};

  Atom wmDeleteMessage_{};
  WinList windows{};

 private:
  std::map<Window, WinIter> window_lookups{};
  Point last_position{-1, -1};
};
using X11Win = X11App::X11Win;

// X11 convenience functions
void draw_centered_oval(Display * display, Window window, GC gc_,
                        const int x, const int y,
                        const int x_rad, const int y_rad) {
  XDrawArc(display, window, gc_, x - x_rad, y - y_rad,
           2 * x_rad + 1, 2 * y_rad + 1, 0, 64 * 360);
}
void fill_centered_oval(Display * display, Window window, GC gc_,
                        const int x, const int y,
                        const int x_rad, const int y_rad) {
  XFillArc(display, window, gc_, x - x_rad, y - y_rad,
           2 * x_rad + 1, 2 * y_rad + 1, 0, 64 * 360);
}
void draw_centered_rectangle(Display * display, Window window, GC gc_,
                             const int x, const int y,
                             const int x_rad, const int y_rad) {
  XDrawRectangle(display, window, gc_, x - x_rad, y - y_rad,
                 2 * x_rad + 1, 2 * y_rad + 1);
}
void fill_centered_rectangle(Display * display, Window window, GC gc_,
                             const int x, const int y,
                             const int x_rad, const int y_rad) {
  XFillRectangle(display, window, gc_, x - x_rad, y - y_rad,
                 2 * x_rad + 2, 2 * y_rad + 2);
}
XRectangle rect(const unsigned int x, const unsigned int y,
                const unsigned int width_, const unsigned int height_) {
  XRectangle rect_;
  rect_.x = x;
  rect_.y = y;
  rect_.width = width_;
  rect_.height = height_;
  return rect_;
}
XRectangle rect(const iBounds & bounds) {
  return rect(bounds[0][0], bounds[1][0], bounds[0][2], bounds[1][2]);
}
bool operator!=(const XRectangle lhs, const XRectangle & rhs) {
  return lhs.x != rhs.x || lhs.y != rhs.y ||
      lhs.width != rhs.width || lhs.height != rhs.height;
}
bool operator!=(const XPoint lhs, const XPoint & rhs) {
  return lhs.x != rhs.x || lhs.y != rhs.y;
}

using void_fun = std::function<void ()>;
using bool_fun = std::function<bool ()>;

struct Actions {
  Actions(void_fun press_ = []() {}, bool_fun visible_ = []() { return true; },
           void_fun release_ = []() {}) :
      press{press_}, visible{visible_}, release{release_} { }

  void_fun press;  // Function to take action on press
  bool_fun visible;  // Function to decide whether to show radio
  void_fun release;  // Function to take action on release
};

class Event {
 public:
  enum EventType {
    Prepare,
    PreDraw,
    Draw,
    X
  };
  explicit Event(const EventType type_ = Draw,
                 const XEvent * x_ = nullptr) :
      type{type_}, x{x_} {}

  EventType type{Draw};
  const XEvent * x{nullptr};
};

// Radio button widget
struct Radio {
  // Constructor
  Radio(const std::string & description_, X11Win * win_,
        const dPoint specification_, const Actions & actions_ = Actions(),
        const bool togglable_ = false, const bool start_toggle = false,
        const GC * gc_ = nullptr) :
      description{description_}, win{win_},
    specification{specification_},
    actions{actions_}, togglable{togglable_}, toggled{start_toggle},
    gc{gc_ == nullptr ? win->radio_gc : *gc_} { }

  // These are not strictly needed - just to avoid warning due to GC pointer
  Radio(const Radio &) = default;
  Radio & operator=(const Radio &) = default;

  // State testing and assignment
  operator bool() const { return toggled; }
  Radio & operator=(const bool state) {
    toggled = state;
    draw();
    return *this;
  }

  // Location of corners of graph region
  Point corner(const bool high_x, const bool high_y) const {
    return Point{high_x ? win->bounds[0][1] : win->bounds[0][0],
          high_y ? win->bounds[1][1] : win->bounds[1][0]};
  }

  int min_border() const {
    return min(win->bounds[0][0], win->bounds[1][0],
               win->width() - win->bounds[0][1],
               win->height() - win->bounds[1][1]);
  }

  // Location in window of radio button
  Point location() const {
    const bool high[2]{specification.x < 0, specification.y < 0};
    const Point anchor{corner(high[0], high[1])};
    Point point;
    // const int border{min_border()};
    const double border{1.0 * min_border()};
    for (const bool y : {false, true}) {
      // Specification near but not equal to zero is at edges
      if (fabs(specification[y]) > 0 && fabs(specification[y]) < 50) {
        point[y] = anchor[y] + border *
            (specification[y] + 0.5 + (high[y] ? 1 : -2));
      } else {
        // Specification of zero or around 100 is centered
        const double centered{specification[y] -
              (fabs(specification[y]) > 0 ? 100 : 0)};
        point[y] = (win->bounds[y][0] + win->bounds[y][1]) / 2 +
            centered * border;
      }
    }
    return point;
  }

  double radius() const { return radius_scale * min_border() / 3.0; }
  bool contains(const Point point) const {
    return location().distance(point) < radius();
  }
  bool visible() const { return actions.visible(); }

  // What to do when radio is pressed
  template <class Point>
  bool press(const Point point) {
    if (!actions.visible()) {
      skip_release = true;
      return contains(point);
    } else if (contains(point)) {
      toggled = !toggled;
      draw();
      actions.press();
      return true;
    }
    return false;
  }

  // What to do when radio is released
  template <class Point>
  bool release(const Point point) {
    // if (skip_release) return skip_release = false;
    if (skip_release) {
      skip_release = false;
      return contains(point);
    }
    if (contains(point)) {
      if (!togglable) toggled = !toggled;
      draw();
      actions.release();
      return true;
    }
    return false;
  }

  void erase() const {
    const Point point{location()};
    fill_centered_oval(win->display(), win->window, win->fill_gc,
                       point.x, point.y, radius() + 1, radius() + 1);
  }

  void draw() const {
    const Point point{location()};
    static GC grey_gc{[this]() {
        XColor grey;
        if (!XAllocNamedColor(win->display(), win->colormap, "rgb:dd/dd/dd",
                              &grey, &grey)) throw Error("Could not get grey");
        return win->create_gc(grey.pixel, win->app.white);
      }()};
    if (win->inside) {
      fill_centered_oval(win->display(), win->window, win->fill_gc,
                         point.x, point.y, radius() + 1, radius() + 1);
      // if (!actions.visible()) return;
      draw_centered_oval(win->display(), win->window,
                         actions.visible() ? gc : grey_gc,
                         point.x, point.y, radius(), radius());
      if (toggled) {
        fill_centered_oval(win->display(), win->window, gc,
                           point.x, point.y, radius() / 2, radius() / 2);
      } else {
        fill_centered_oval(win->display(), win->window, win->fill_gc,
                           point.x, point.y, radius() / 2, radius() / 2);
      }
    } else {
      erase();
    }
  }

  // Data
  std::string description;  // Help text for radio
  X11Win * win;  // The window attached to
  dPoint specification;  // Where on page
  Actions actions;  // Actions to perform
  bool togglable;  // Can radio be toggled
  bool toggled;  // State of radio
  GC gc;  // Color for radio, line width, etc
  bool skip_release{false};
  double radius_scale{1.0};
  unsigned int id{0};
};

using Range = std::vector<std::vector<double>>;

class SavedConfig {
 public:
  static constexpr double default_arc_radius{4};
  static constexpr double default_arc_width{2};
  static constexpr int default_line_width{4};
  static constexpr int default_line_type{LineSolid};

  SavedConfig() {}
  virtual ~SavedConfig() = default;
  bool operator!=(const SavedConfig & rhs) const {
    return dne(arc_radius, rhs.arc_radius) ||
        dne(arc_width, rhs.arc_width) ||
        line_width != rhs.line_width ||
        line_type != rhs.line_type ||
        series_order != rhs.series_order ||
        range != rhs.range ||
        max_range != rhs.max_range ||
        zoomed != rhs.zoomed ||
        drawn != rhs.drawn ||
        radio_states != rhs.radio_states;
  }
  void restore_config(const SavedConfig & rhs) {
    *this = rhs;
    return;
  }

  double arc_radius{default_arc_radius};
  double arc_width{default_arc_width};
  int line_width{default_line_width};
  int line_type{default_line_type};
  std::vector<unsigned int> series_order{};
  Range range{{unset(1.0), nunset(1.0), 0}, {unset(1.0), nunset(1.0), 0}};
  Range max_range{range};
  std::vector<unsigned char> zoomed{false, false};
  mutable bool drawn{false};
  std::vector<unsigned char> radio_states{};
};

class X11Graph : public X11Win, public SavedConfig {
 public:
  static constexpr int border_width{3};

  // Graph data typedefs
  using Values = std::vector<double>;
  using XYSeries = std::vector<const Values *>;
  using Data = std::vector<XYSeries>;
  using CallBack = std::function<bool (X11Graph &, Event &)>;

  // Graph plotting window
  // Creation factory from data in exact format needed
  static X11Graph & create_whole(X11App & app, const Data & data__,
                                 const unsigned int width_ = 1200,
                                 const unsigned int height_ = 1000,
                                 const int x_off_ = 0,
                                 const int y_off_ = 0) {
    return reinterpret_cast<X11Graph &>(app.add(
        std::make_unique<X11Graph>(app, data__,
                                   width_, height_, x_off_, y_off_)));
  }

  // Creation factory from a bunch of vectors x1, y1, x2, y2, ...
  template <class ... Input>
  static X11Graph & create(X11App & app, Input && ... input) {
    return reinterpret_cast<X11Graph &>(app.add(
        std::make_unique<X11Graph>(app, std::forward<Input>(input)...)));
  }

  // Construct from data in exact format needed
  X11Graph(X11App & app__, const Data & data__,
           const unsigned int width_ = 1200,
           const unsigned int height_ = 1000,
           const int x_off_ = 0,
           const int y_off_ = 0) :
      X11Win{app__, width_, height_, x_off_, y_off_}, input_data{data__},
    data{&input_data} {
        initialize();
      }

  // Construct from a bunch of vectors x1, y1, x2, y2, ...
  template <class ... Input>
  X11Graph(X11App & app__, Input && ... input) : X11Win{app__} {
    add_input(std::forward<Input>(input)...);
    data = &input_data;
    initialize();
  }

  X11Graph(const X11Graph &) = delete;
  X11Graph & operator=(const X11Graph &) = delete;

  // Add a bunch of series (only during construction right now)
  template <class ... Input>
  void add_input(Values & x__, Values & y__, Input && ... input) {
    input_data.emplace_back(0);
    input_data.back().push_back(&x__);
    input_data.back().push_back(&y__);
    add_input(std::forward<Input>(input)...);
  }
  void add_input() { }

  void add_call_back(const std::string & help_text,
                     const CallBack & call_back,
                     const bool full_draw = false,
                     const bool initially_on = true) {
    call_back_radios.reserve(100);
    call_backs.reserve(100);
    call_backs.push_back(call_back);
    call_back_radios.push_back(Radio{help_text, this,
        {1, call_backs.size() + 3.0}, {[this, full_draw]() {
            return full_draw ? draw() : redraw(); }}, true, initially_on});
    radios.push_back(&call_back_radios.back());
    // call_back_radios.back().draw();
  }

  void initialize();

  void set_line_widths(std::vector<GC> gcs, const int width_) {
    for (unsigned int g{0}; g != data->size(); ++g)
      XSetLineAttributes(display(), gcs[g], width_, LineSolid,
                         CapButt, JoinRound);
  }

  virtual bool slow() const { return true; }

  // Get range for axes x 0, y 1, or both 2
  void set_range(const bool y, const double low, const double high) {
    // range[y][0] = restrict_range_radios[y] ? max(low, max_range[y][0]) : low;
    // range[y][1] =
    // restrict_range_radios[y] ? min(high, max_range[y][1]) : high;
    range[y][0] = low;
    range[y][1] = high;
    range[y][2] = range[y][1] - range[y][0];
    zoomed[y] = (dne(range[y][0], max_range[y][0]) ||
                 dne(range[y][1], max_range[y][1]));
  }
  void get_range(const unsigned int a = 2);
  void range_jump(const bool y, const double dist) {
    set_range(y, range[y][0] + dist, range[y][1] + dist);
  }
  bool in_range(const double x, const double y) const {
    return x >= range[0][0] && x <= range[0][1] &&
        y >= range[1][0] && y <= range[1][1];
  }
  bool in_range(const dPoint pos) const { return in_range(pos.x, pos.y); }

  // Data to window coordinate transformation
  int coord(const bool y, const double val) const {
    if (y) return bounds[1][1] - (val - range[1][0]) * scale[1];
    return bounds[0][0] + (val - range[0][0]) * scale[0];
  }
  Point coord(const dPoint point) const {
    return Point{coord(0, point.x), coord(1, point.y)};
  }
  XPoint xcoord(const dPoint point) const {
    XPoint result;
    result.x = coord(0, point.x);
    result.y = coord(1, point.y);
    return result;
  }
  XPoint xcoord(const Point point) const {
    XPoint result;
    result.x = point.x;
    result.y = point.y;
    return result;
  }

  // Window to data coordinate transformation
  double icoord(const bool y, const int val) const {
    if (y) return (bounds[1][1] - val) / scale[1] + range[1][0];
    return (val - bounds[0][0]) / scale[0] + range[0][0];
  }
  dPoint icoord(const Point point) const {
    return dPoint{icoord(0, point.x), icoord(1, point.y)};
  }

  int min_border() const { return 0.05 * std::min(width(), height()); }
  unsigned int get_quadrant(const Point point) const;

  //
  // Event loop callbacks
  //
  virtual void expose(const XExposeEvent & event) {
    if (bounds.size() == 0) {
      get_range();
      prepare();
    }
    if (drawn) {
      XCopyArea(display(), pixmap, window, gc, event.x, event.y,
                event.width, event.height, event.x, event.y);
      draw_controls();
    } else {
      draw();
    }
  }
  virtual void enter(const XCrossingEvent &) {
    inside = true;
    erase_border();
    draw_controls();
  }
  virtual void key(const XKeyEvent & evnt) { last_motion = last_press = evnt; }
  virtual void button_press(const XButtonEvent & event);
  virtual void motion(const XMotionEvent & event);
  virtual void button_release(const XButtonEvent & event);
  virtual void leave(const XCrossingEvent &);

  bool do_arcs(const unsigned int s) const {
    if (!series_radios[s]) return false;
    return (arcs_radio && !series_only_lines[s]) || series_only_arcs[s];
  }
  bool do_arcs() const {
    for (unsigned int s{0}; s != series_only_arcs.size(); ++s)
      if (do_arcs(s)) return true;
    return false;
  }
  bool can_do_arcs() const {
    for (unsigned int s{0}; s != series_only_arcs.size(); ++s)
      if (series_radios[s] && !series_only_lines[s]) return true;
    return false;
  }
  bool do_lines(const unsigned int s) const {
    if (!series_radios[s]) return false;
    return (lines_radio && !series_only_arcs[s]) || series_only_lines[s];
  }
  bool do_lines() const {
    for (unsigned int s{0}; s != series_only_lines.size(); ++s)
      if (do_lines(s)) return true;
    return false;
  }
  bool can_do_lines() const {
    for (unsigned int s{0}; s != series_only_lines.size(); ++s)
      if (series_radios[s] && !series_only_arcs[s]) return true;
    return false;
  }

  void prepare_log();
  virtual void prepare();
  void draw_ticks() {
    if (inside) return;
    if (!tick_radios[0] && !tick_radios[1]) return;
    static std::vector<std::string> tick_labels;
    tick_labels.clear();
    if (0)
      std::cerr << "Bounds "
                << bounds[0][0] << " " << bounds[0][1] << " "
                << bounds[1][0] << " " << bounds[1][1] << std::endl;
    static GC tick_gc{create_gc(app.black, app.white)};
    static X11Font * tick_font{nullptr};
    const double avail_height{bounds[1][0] * 0.6};
    X11Font * fits{app.fonts.fits("moo", bounds[0][2], avail_height)};
    if (fits != tick_font) {
      tick_font = fits;
      XSetFont(display(), tick_gc, fits->id());
    }
    const int t_height{tick_font->height()};
    for (const bool y : {false, true}) {
      if (!tick_radios[y]) continue;
      const Axis axis{range[y][0], range[y][1], 3, log_radios[y]};
      for (const std::pair<double, bool> tick : axis.ticks()) {
        if (!tick.second) continue;
        const int loc{coord(y, tick.first)};
        std::ostringstream label;
        label << std::setprecision(6)
              << (log_radios[y] ? pow(10, tick.first) : tick.first);
        tick_labels.push_back(label.str());
        std::string & text{tick_labels.back()};
        const int t_width{tick_font->string_width(text)};
        XDrawString(display(), window, tick_gc,
                    y ? std::max(0, bounds[0][0] - t_width - 3) :
                    loc - t_width / 2,
                    y ? tick_font->centered_y(loc) : bounds[1][1] + t_height,
                    text.c_str(),
                    static_cast<unsigned int>(text.size()));
      }
    }
  }
  void set_clip_rectangle(const unsigned int x, const unsigned int y,
                          const unsigned int width_,
                          const unsigned int height_);

  void draw_grid() const;
  void erase_border() {
    XFillRectangle(display(), window, fill_gc,
                   0, bounds[1][1], width(), height());
    XFillRectangle(display(), window, fill_gc,
                   0, 0, bounds[0][0], height());
  }
  void draw_controls() {
    XDrawRectangle(display(), window, border_gc,
                   bounds[0][0], bounds[1][0], bounds[0][2], bounds[1][2]);
    draw_status();
    draw_grid();
    for (const Radio * radio : radios) radio->draw();

    // Draw tick labels
    draw_ticks();

    // Handle special drawing commands
    Event nothing;
    for (unsigned int c{0}; c != call_backs.size(); ++c)
      if (call_back_radios[c]) call_backs[c](*this, nothing);
  }

  static std::string long_status(const bool in, const bool y) {
    return std::string("Pointer (1 - 2/shift - 3/control) clicks ") +
        "(center - zoom in - zoom out) at point " +
        "and drags (select - scroll - zoom) for " +
        (in ? "X and Y axes" : (y ? "Y axis" : "X axis"));
  }

  void draw_status(const bool force = false) const {
    XFillRectangle(display(), window, fill_gc, bounds[0][0], 0,
                   bounds[0][1] - bounds[0][0], bounds[1][0] - border_width);
    if (force || help_radio || coord_radio) {
      const double avail_height{bounds[1][0] * 0.65};
      static X11Font * status_font{nullptr};
      X11Font * fits{app.fonts.fits(status, bounds[0][2], avail_height)};
      if (fits != status_font) {
        status_font = fits;
        XSetFont(display(), gc, fits->id());
      }
      XDrawString(display(), window, gc, bounds[0][0],
                  fits->centered_y((bounds[1][0] - border_width) / 2),
                  const_cast<char *>(status.c_str()),
                  static_cast<unsigned int>(status.size()));
    }
  }

  double line_vertical_y(const dPoint low_x, const dPoint high_x,
                         const double x) const {
    const double slope((high_x.y - low_x.y) / (high_x.x - low_x.x));
    return low_x.y + (x - low_x.x) * slope;
  }
  double line_horizontal_x(const dPoint low_x, const dPoint high_x,
                           const double y) const {
    const double slope((high_x.y - low_x.y) / (high_x.x - low_x.x));
    return low_x.x + (y - low_x.y) / slope;
  }

  XPoint line_bounds_intersection(const dPoint in, const dPoint out) const {
    // need to do better to include lines between points outside range
    // but intersecting the range.  also need to reduce unnecessary
    // computation here and need to place points outside of range...

    const bool out_high[2]{out[0] > range[0][1], out[1] > range[1][1]};
    const bool out_low[2]{out[0] < range[0][0], out[1] < range[1][0]};
    const bool is_out[2]{out_high[0] || out_low[0], out_high[1] || out_low[1]};
    const dPoint limit{range[0][out_high[0]], range[1][out_high[1]]};
    if (!dne(in.x, out.x)) return xcoord(dPoint{in.x, limit.y});
    if (!dne(in.y, out.y)) return xcoord(dPoint{limit.x, in.y});
    const double slope((out.y - in.y) / (out.x - in.x));
    const dPoint solutions{
      in.x + (limit.y - in.y) / slope, in.y + (limit.x - in.x) * slope};
    const dPoint trials[2]{{solutions.x, limit.y}, {limit.x, solutions.y}};
    const dPoint distance{in.distance(trials[0]), in.distance(trials[1])};
    const bool best_is_x{(is_out[0] && is_out[1]) ? distance.x < distance.y :
          is_out[1]};
    const dPoint intersection{best_is_x ? trials[0] : trials[1]};
    return xcoord(intersection);
  }

  virtual void draw();
  void redraw() {
    XCopyArea(display(), pixmap, window, gc, bounds[0][0], bounds[1][0],
              bounds[0][2], bounds[1][2], bounds[0][0], bounds[1][0]);
    draw_controls();
  }

  virtual void save_image(const std::string & base_name,
                          std::function<void()> call_back =
                          std::function<void()>()) {
    if (!call_back) {
      call_back = [this]() {
        status = "Saving Image";
        draw_controls();
        draw_status(true);
        XFlush(display());
      };
    }
    inside = false;
    const bool help_state = help_radio;
    help_radio = false;
    draw_controls();
    X11Win::save_image(base_name, call_back);
    inside = true;
    help_radio = help_state;
    status = "Done saving image";
    draw_controls();
    draw_status(true);
  }

  virtual ~X11Graph() {
    for (GC gc_ : {border_gc, minor_gc, major_gc}) XFreeGC(display(), gc_);
    // FIX Should also free all other GCs...
  }

  // Bug in code? sometimes movie cannot be started until after zoom out
  bool movie(const bool right) {
    status = "Playing the movie - click movie radio button again to stop";
    using Time = std::chrono::time_point<std::chrono::system_clock>;
    Time last_time{std::chrono::system_clock::now()};
    const double rate{0.5};  // page per second
    XEvent event;
    XWindowEvent(display(), window, ButtonReleaseMask, &event);
    while (true) {
      Time time{std::chrono::system_clock::now()};
      const double seconds{
        std::chrono::duration_cast<std::chrono::milliseconds>(
            time - last_time).count() / 1000.0};
      last_time = time;
      const double movement{rate * seconds * range[0][2]};
      range_jump(0, (right ? 1 : -1) * movement);
      prepare_draw();
      XFlush(display());
      if (!(range[0][0] > max_range[0][0] && range[0][1] < max_range[0][1])) {
        movie_radios[right] = false;
        return true;
      }
      if (XCheckWindowEvent(display(), window, ButtonPressMask, &event)) {
        XWindowEvent(display(), window, ButtonReleaseMask, &event);
        if (movie_radios[right].release(event.xbutton)) {
          movie_radios[right] = false;
          return true;
        }
      }
    }
  }

  void show_order(const std::string prefix) const {
    std::cout << prefix;
    for (const unsigned int s : series_order) {
      std::cout << " " << s;
    }
    std::cout << std::endl;
  }
  void show_range(const std::string prefix) const {
    std::cout << prefix << " range";
    for (const bool y : {false, true}) {
      for (const double val : range[y]) {
        std::cout << " " << val;
      }
    }
    if (bounds.size()) {
      std::cout << " bounds";
      for (const bool y : {false, true}) {
        for (const double val : bounds[y]) {
          std::cout << " " << val;
        }
      }
    }
    std::cout << " scale";
    for (const double val : scale) {
      std::cout << " " << val;
    }
    std::cout << std::endl;
  }



  // Data is series * (x, y) -> point
  Data input_data{}, log_data{}, log_x_data{}, log_y_data{}, *data{};
  std::vector<std::unique_ptr<Values> > log_series{};

  // Drawing data for series
  std::vector<std::vector<XArc> > arcs{};
  std::vector<std::vector<XPoint> > points{};

  // Graphics contexts
  GC border_gc{}, border_fill_gc{}, minor_gc{}, major_gc{};

  // Colors for series
  std::vector<std::string> color_names{
    "rgb:e5/00/00", "rgb:25/00/9e", "rgb:00/b7/00", "rgb:e5/be/00",
        "rgb:06/56/93", "rgb:b7/dd/00", "rgb:e5/83/00", "rgb:95/00/95",
        "red", "blue", "orange", "green",
        "brown", "purple", "cyan", "yellow", "black", "gray",
        "rgb:e5/00/00", "rgb:25/00/9e", "rgb:00/b7/00", "rgb:e5/be/00",
        "rgb:06/56/93", "rgb:b7/dd/00", "rgb:e5/83/00", "rgb:95/00/95",
        "red", "blue", "orange", "green",
        "brown", "purple", "cyan", "yellow", "black", "gray",
        "rgb:e5/00/00", "rgb:25/00/9e", "rgb:00/b7/00", "rgb:e5/be/00",
        "rgb:06/56/93", "rgb:b7/dd/00", "rgb:e5/83/00", "rgb:95/00/95",
        "red", "blue", "orange", "green",
        "brown", "purple", "cyan", "yellow", "black", "gray",
        "rgb:e5/00/00", "rgb:25/00/9e", "rgb:00/b7/00", "rgb:e5/be/00",
        "rgb:06/56/93", "rgb:b7/dd/00", "rgb:e5/83/00", "rgb:95/00/95",
        "red", "blue", "orange", "green",
        "brown", "purple", "cyan", "yellow", "black", "gray",
    "rgb:e5/00/00", "rgb:25/00/9e", "rgb:00/b7/00", "rgb:e5/be/00",
        "rgb:06/56/93", "rgb:b7/dd/00", "rgb:e5/83/00", "rgb:95/00/95",
        "red", "blue", "orange", "green",
        "brown", "purple", "cyan", "yellow", "black", "gray",
        "rgb:e5/00/00", "rgb:25/00/9e", "rgb:00/b7/00", "rgb:e5/be/00",
        "rgb:06/56/93", "rgb:b7/dd/00", "rgb:e5/83/00", "rgb:95/00/95",
        "red", "blue", "orange", "green",
        "brown", "purple", "cyan", "yellow", "black", "gray",
        "rgb:e5/00/00", "rgb:25/00/9e", "rgb:00/b7/00", "rgb:e5/be/00",
        "rgb:06/56/93", "rgb:b7/dd/00", "rgb:e5/83/00", "rgb:95/00/95",
        "red", "blue", "orange", "green",
        "brown", "purple", "cyan", "yellow", "black", "gray",
        "rgb:e5/00/00", "rgb:25/00/9e", "rgb:00/b7/00", "rgb:e5/be/00",
        "rgb:06/56/93", "rgb:b7/dd/00", "rgb:e5/83/00", "rgb:95/00/95",
        "red", "blue", "orange", "green",
        "brown", "purple", "cyan", "yellow", "black", "gray",
    "rgb:e5/00/00", "rgb:25/00/9e", "rgb:00/b7/00", "rgb:e5/be/00",
        "rgb:06/56/93", "rgb:b7/dd/00", "rgb:e5/83/00", "rgb:95/00/95",
        "red", "blue", "orange", "green",
        "brown", "purple", "cyan", "yellow", "black", "gray",
        "rgb:e5/00/00", "rgb:25/00/9e", "rgb:00/b7/00", "rgb:e5/be/00",
        "rgb:06/56/93", "rgb:b7/dd/00", "rgb:e5/83/00", "rgb:95/00/95",
        "red", "blue", "orange", "green",
        "brown", "purple", "cyan", "yellow", "black", "gray",
        "rgb:e5/00/00", "rgb:25/00/9e", "rgb:00/b7/00", "rgb:e5/be/00",
        "rgb:06/56/93", "rgb:b7/dd/00", "rgb:e5/83/00", "rgb:95/00/95",
        "red", "blue", "orange", "green",
        "brown", "purple", "cyan", "yellow", "black", "gray",
        "rgb:e5/00/00", "rgb:25/00/9e", "rgb:00/b7/00", "rgb:e5/be/00",
        "rgb:06/56/93", "rgb:b7/dd/00", "rgb:e5/83/00", "rgb:95/00/95",
        "red", "blue", "orange", "green",
        "brown", "purple", "cyan", "yellow", "black", "gray"};
#if 0
  color_names.reserve(color_names.size() * 4);
  color_names.insert(color_names.end(), color_names.begin(), color_names.end());
  color_names.insert(color_names.end(), color_names.begin(), color_names.end());
#endif
  std::vector<XColor> series_colors{color_names.size()};
  std::vector<GC> series_arc_gcs{color_names.size()};
  std::vector<GC> series_line_gcs{color_names.size()};
  std::vector<GC> series_radio_gcs{color_names.size()};
  std::vector<Radio> series_radios{};
  std::vector<uint8_t> series_only_arcs{}, series_only_lines{};

  // Graph state information
  std::string status{""};
  std::vector<double> scale{};  // x, y, y / x
  Point last_press{}, last_motion{};
  bool moved{false};
  bool small_move{false};

  bool_fun radio_tester(const Radio & radio, const bool state = true) {
    return [&radio, state]() { return radio.toggled == state; };
  }
  bool_fun zoom_tester(const bool y) {
    return [this, y]() { return zoomed[y]; };
  }

  //
  // Radio controls
  //
  Radio help_radio{"Toggle showing help text for controls", this, {1, 2},
    {[this]() { coord_radio = false; draw_controls(); }}, true, true};
  Radio coord_radio{"Toggle showing coordinates of cursor", this, {1, 3},
    {[this]() { help_radio = false; status = ""; draw_controls(); }}, true};
  Radio arcs_radio{"Draw a marker at each graph point", this, {-1, -2},
    {[this]() { return arcs_radio ? prepare_draw() : draw(); },
          [this]() { return can_do_arcs(); }}, true, true};
  Radio outlines_radio{"Toggle between solid and outlined markers", this,
    {-1, -5.5}, {[this]() { draw(); }, [this]() { return do_arcs(); }}, true};
  Radio lines_radio{"Connect graph points by lines", this, {-2, -1},
    {[this]() { return lines_radio ? prepare_draw() : draw(); },
          [this]() { return can_do_lines(); }},  true};
  Radio tick_radios[2]{
    {"Toggle axis labels on X axis (shown when cursor leaves window)",
          this, {5.5, -1},
      {[this]() { }}, true},
    {"Toggle axis labels on Y axis (shown when cursor leaves window)",
          this, {1, -5.5},
      {[this]() { }}, true}};
  Radio log_radios[2]{
    {"Toggle logarithmic scale on X axis", this, {6.5, -1},
      {[this]() { prepare_log(); prepare_draw(); }}, true},
    {"Toggle logarithmic scale on Y axis", this, {1, -6.5},
      {[this]() { prepare_log(); prepare_draw(); }}, true}};
  Radio grid_radios[2][2]{
    {{"Toggle major grid lines on X axis", this, {4.25, -1},
        {[this]() { if (grid_radios[0][0]) { return draw_grid(); }
            grid_radios[1][0] = false; redraw(); }}, true, true},
      {"Toggle major grid lines on Y axis", this, {1, -4.25},
        {[this]() { if (grid_radios[0][1]) { return draw_grid(); }
            grid_radios[1][1] = false; redraw(); }}, true, true}},
    {{"Toggle minor grid lines on X axis", this, {3.25, -1},
        {[this]() { if (grid_radios[1][0]) { grid_radios[0][0] = true;
              draw_grid(); } redraw(); }}, true, true},
      {"Toggle minor grid lines on Y axis", this, {1, -3.25},
        {[this]() { if (grid_radios[1][1]) { grid_radios[0][1] = true;
              draw_grid(); } redraw(); }}, true, true}}};
  Radio movie_radios[2]{
    {"Play a movie traveling left", this, {97.5, -1},
      {[this]() { movie(false); }, [this]() { return zoomed[0]; }}, true},
    {"Play a movie traveling right", this, {102.5, -1},
      {[this]() { movie(true); }, [this]() { return zoomed[0]; }}, true}};
  Radio restrict_range_radios[2]{
    {"Toggle range restriction on X axis to actual data range", this, {3, -1},
      {[this]() { get_range(0); prepare_draw(); }}, true, true},
    {"Toggle range restriction on Y axis to actual data range", this, {1, -3},
      {[this]() { get_range(1); prepare_draw(); }}, true, true}};

  // Saved history
  std::deque<SavedConfig> saved_config{};
  std::vector<Radio *> saved_radios{
        &arcs_radio, &outlines_radio, &lines_radio};
  SavedConfig current_config() const {
    SavedConfig current{*this};
    current.radio_states.clear();
    for (Radio * radio : saved_radios) {
      current.radio_states.push_back(*radio);
    }
    return current;
  }

  void restore_config(const SavedConfig & config) {
    if (dne(config.line_width, line_width)) {
      set_line_widths(series_line_gcs,
                      (config.line_width == 1 ? 0 : config.line_width));
    }
    if (dne(config.arc_width, arc_width)) {
      set_line_widths(series_arc_gcs, config.arc_width);
    }
    for (unsigned int r{0}; r != saved_radios.size(); ++r) {
      saved_radios[r]->toggled = config.radio_states[r];
    }
    SavedConfig::restore_config(config);
  }
  void save_config(const SavedConfig & config) {
    saved_config.push_back(std::move(config));
  }

  // Cannot go into saved config
  Radio previous_views_radio{"Show previous view",
        this, {-1, 1}, {[this]() { },
          [this]() {return saved_config.size() > 1; },
          {[this]() {
              if (saved_config.size() > 1) {
                saved_config.pop_back();
                restore_config(saved_config.back());
                prepare_draw();
              }}}}};

  // All radios
  std::deque<Radio *> radios{&help_radio, &coord_radio,
        &arcs_radio, &outlines_radio, &lines_radio,
        &tick_radios[0], &tick_radios[1],
        &log_radios[0], &log_radios[1],
        &grid_radios[0][0], &grid_radios[0][1],
        &grid_radios[1][0], &grid_radios[1][1],
        &movie_radios[0], &movie_radios[1],
        &previous_views_radio};
  //        &restrict_range_radios[0], &restrict_range_radios[1]};
  std::vector<Radio> create_unnamed_radios();
  std::vector<Radio> extra_radios{};

  // Callback functions and radios for customization
  std::vector<CallBack> call_backs{};
  std::vector<Radio> call_back_radios{};
  std::vector<Radio> unnamed_radios{};

  // Number of threads to use
  unsigned int n_threads_{std::thread::hardware_concurrency()};
  unsigned int n_threads(const unsigned int n_threads__ = 0) {
    if (n_threads__ != 0) n_threads_ = n_threads__;
    return n_threads_;
  }
};

// Common initialization for constructors
inline void X11Graph::initialize() {
  if (data->size() > color_names.size())
    throw Error(std::string("Too many series to display (max is ") +
                std::to_string(color_names.size()) + ")");

  scale.resize(3);

  // Events to watch out for
  XSelectInput(display(), window, StructureNotifyMask | ExposureMask |
               KeyPressMask | EnterWindowMask | LeaveWindowMask |
               ButtonPressMask | PointerMotionMask | ButtonReleaseMask);

  // Graphics contexts
  border_gc = create_gc(app.black, app.white, border_width,
                        LineSolid, CapButt, JoinMiter);
  border_fill_gc = create_gc(app.white, app.black, border_width,
                             LineSolid, CapButt, JoinMiter);
  minor_gc = create_gc(app.black, app.white, 1, LineOnOffDash,
                       CapButt, JoinMiter);
  major_gc = create_gc(app.black, app.white, 2, LineOnOffDash,
                       CapButt, JoinMiter);

  // Series colors
  for (unsigned int c{0}; c != color_names.size(); ++c) {
    std::string & color_name{color_names[c]};
    XColor & color{series_colors[c]};
    if (!XAllocNamedColor(display(), colormap, color_name.c_str(),
                          &color, &color))
      throw Error("Could not get color") << color_name;

    // For colored arcs for points
    series_arc_gcs[c] = create_gc(color.pixel, app.white, arc_width,
                                  LineSolid, CapButt, JoinMiter);

    // For colored lines to connect points
    series_line_gcs[c] = create_gc(color.pixel, app.white, line_width,
                                   line_type, CapProjecting, JoinRound);

    // For colored thin lines to outline radios
    series_radio_gcs[c] = create_gc(color.pixel, app.white, radio_width,
                                   LineSolid, CapButt, JoinMiter);
  }

  // Create series radios, and add to master radio list, adjust properties
  series_radios.reserve(data->size());
  for (unsigned int c{0}; c != data->size(); ++c) {
    series_radios.push_back(Radio{"Toggle display of this series", this,
        {-1, data->size() + 1.0 - c},
        {[this]() { prepare_draw(); }, [this]() { return inside; }},
            true, true, &series_radio_gcs[c]});
    saved_radios.push_back(&series_radios.back());
    series_only_arcs.emplace_back(0);
    series_only_lines.emplace_back(0);
    series_order.push_back(c);
  }
  series_order.reserve(series_order.size() + 1);

  // Add radios to master radio list
  unnamed_radios = create_unnamed_radios();

  for (std::vector<Radio> * radio_vec : {&series_radios, &unnamed_radios})
    for (Radio & radio : (*radio_vec)) radios.push_back(&radio);
  extra_radios.reserve(1000);
}

// Get range for axes a: x 0, y 1, or both 2
inline void X11Graph::get_range(const unsigned int a) {
  constexpr double padding{0.01};
  for (const bool y : {0, 1}) {
    if (a != 2 && a != y) continue;
    range[y] = {unset(1.0), nunset(1.0), 0};
    for (unsigned int s{0}; s != data->size(); ++s) {
      if (!series_radios[s]) continue;
      for (const double val : *(*data)[s][y]) {
        if (!std::isfinite(val)) continue;
        if (range[y][0] > val) range[y][0] = val;
        if (range[y][1] < val) range[y][1] = val;
      }
    }
    range[y][2] = range[y][1] - range[y][0];
    range[y][0] -= padding * range[y][2];
    range[y][1] += padding * range[y][2];
    range[y][2] = range[y][1] - range[y][0];
    zoomed[y] = false;
  }
  if (a == 2) max_range = range;
}

// Which quadrant of graph is a point in
inline unsigned int X11Graph::get_quadrant(const Point point) const {
  const bool below_pos{point.y > bounds[1][1] + (point.x - bounds[0][0]) *
        (bounds[1][0] - bounds[1][1]) / (bounds[0][1] - bounds[0][0])};
  const bool below_neg{point.y > bounds[1][0] + (point.x - bounds[0][0]) *
        (bounds[1][1] - bounds[1][0]) / (bounds[0][1] - bounds[0][0])};
  return below_pos ? (below_neg ? 0 : 3) : (below_neg ? 1 : 2);
}

inline void X11Graph::button_press(const XButtonEvent & event) {
  last_motion = last_press = event;
  moved = false;
  for (Radio * radio : radios) { if (radio->press(event)) return; }

  // Button actions - no work done here other than record last press
  // Inside graph: both dimensions.  Outside graph: one dimension
  //
  // key   button : press only (on release)  / drag or move (during)
  // --------------------------------------------------------------------
  // none        1: center                   / select a new view
  // shift   or  2: center and zoom in       / scroll
  // control or  3: center and zoom out      / zoom
}

inline void X11Graph::motion(const XMotionEvent & event) {
  // std::cerr << "In motion" << std::endl;
  moved = true;
  if (XPending(display())) return;

  Event motion_event{Event::X, &app.event};
  bool call_back_acted{false};
  for (unsigned int c{0}; c != call_backs.size(); ++c) {
    if (call_back_radios[c]) {
      if (call_backs[c](*this, motion_event)) {
        call_back_acted = true;
        break;
      }
    }
  }

  if (!call_back_acted) {
    // Status text
    if (!help_radio && help_radio.contains(event)) {
      status = help_radio.description;
      return draw_status(true);
    }
    status = "";
    if (help_radio) {
      for (Radio * radio : radios) {
        // if (radio->visible()) 0;
        if (radio->contains(event)) {
          status = radio->description;
          if (!radio->visible()) status += " (inactive)";
          break;
        }
      }
      if (status.empty())
        status = long_status(in_bounds(event), get_quadrant(event) % 2);
    } else if (coord_radio && in_bounds(event)) {
      std::ostringstream coordinates;
      coordinates << std::setprecision(12) << "(";
      const Point point{event};
      for (const bool y : {false, true}) {
        const double val{icoord(y, point[y])};
        const double res{range[y][2] / bounds[y][2]};
        const double pres{pow(10, floor(log10(res)))};
        const double rval{round(val / pres) * pres};
        const double nval{log_radios[y] ? pow(10, rval) : rval};
        coordinates << (y ? " , " : " ") << nval;
      }
      coordinates << " )";
      status = coordinates.str();
    }
    draw_status();
  }
  if (event.state == 0) return;
  for (Radio * radio : radios) if (radio->contains(last_press)) return;

  if (!(event.state & (Button1Mask | Button2Mask | Button3Mask))) return;

  const Point point{event};
  const unsigned int quadrant{get_quadrant(last_press)};
  const bool y_press{(quadrant % 2) == 1};
  const Range old_range(range);

  const bool scroll{event.state & Button2Mask || event.state & ShiftMask};
  const bool zoom{event.state & Button3Mask || event.state & ControlMask};
  const bool select{event.state & Button1Mask && !scroll && !zoom};

  if (scroll) {
    for (const bool y : {false, true}) {
      if (!in_bounds(last_press) && y_press != y) continue;
      const int distance{point[y] - last_motion[y]};
      const double move{(y ? 1 : -1) * distance / scale[y]};
      range_jump(y, move);
    }
  } else if (select) {
    draw_controls();
    if (in_bounds(last_press)) {
      const Point min_point{min(last_motion.x, last_press.x, point.x),
            min(last_motion.y, last_press.y, point.y)};
      const Point max_point{max(last_motion.x, last_press.x, point.x),
            max(last_motion.y, last_press.y, point.y)};
      // Cover vertical lines
      const int y_start{min(last_motion.y, last_press.y)};
      const int y_height{abs(last_motion.y - last_press.y) + 1};
      XCopyArea(display(), pixmap, window, gc, last_motion.x, y_start,
                1, y_height, last_motion.x, y_start);
      XCopyArea(display(), pixmap, window, gc, last_press.x, y_start,
                1, y_height, last_press.x, y_start);

      // Cover horizontal lines
      const int x_start{min(last_motion.x, last_press.x)};
      const int x_width{abs(last_motion.x - last_press.x) + 1};
      XCopyArea(display(), pixmap, window, gc, x_start, last_motion.y,
                x_width, 1, x_start, last_motion.y);
      XCopyArea(display(), pixmap, window, gc, x_start, last_press.y,
                x_width, 1, x_start, last_press.y);
      XDrawRectangle(display(), window, gc,
                     min(last_press.x, point.x), min(last_press.y, point.y),
                     abs(last_press.x - point.x), abs(last_press.y - point.y));
    } else {
      const bool above(quadrant == 0 || quadrant == 3);
      const int loc{bounds[!y_press][above] + (above ? 2 : -2) * border_width};
      XDrawLine(display(), window, border_fill_gc,
                y_press ? loc : last_press.x, y_press ? last_press.y : loc,
                y_press ? loc : last_motion.x, y_press ? last_motion.y : loc);
      XDrawLine(display(), window, border_gc,
                y_press ? loc : last_press.x, y_press ? last_press.y : loc,
                y_press ? loc : point.x, y_press ? point.y : loc);
    }
  } else if (zoom) {
    for (const bool y : {false, true}) {
      if (!in_bounds(last_press) && y_press != y) continue;
      const int distance{point[y] - last_motion[y]};
      const double change{(y ? 1 : -1) * range[y][2] *
            distance / bounds[y][2]};
      set_range(y, range[y][0] - change, range[y][1] + change);
    }
  }
  last_motion = point;
  if (range != old_range) {
    small_move = true;
    return prepare_draw();
  }
}

inline void X11Graph::button_release(const XButtonEvent & event) {
  for (Radio * radio : radios) { if (radio->release(last_press)) return; }

  Event button_event{Event::X, &app.event};
  for (unsigned int c{0}; c != call_backs.size(); ++c) {
    if (call_back_radios[c]) {
      if (call_backs[c](*this, button_event)) return;
    }
  }

  const Point release{event};
  const unsigned int quadrant{get_quadrant(last_press)};
  const bool y_press{(quadrant % 2) == 1};
  const Range old_range(range);
  const bool shift(event.state & ShiftMask);
  const bool control(event.state & ControlMask);

  if (moved) {
    if (event.button == Button1 && !shift && !control) {
      // Drag event defines an x, y zoom
      for (const bool y : {false, true}) {
        if (!in_bounds(last_press) && y_press != y) continue;
        const double min_c{icoord(y, std::min(release[y], last_press[y]))};
        const double max_c{icoord(y, std::max(release[y], last_press[y]))};
        set_range(y, (y ? max_c : min_c), (y ? min_c : max_c));
      }
    }
    moved = false;
  } else {
    // Button 1, 2, 3 click only behavior
    const bool in{event.button == Button2 || shift};
    const bool out{event.button == Button3 || control};
    const bool center{event.button == Button1 && !out && !in};
    if (center || in || out) {
      for (const bool y : {false, true}) {
        // if ((center || out) && range[y] == max_range[y]) continue;
        if (!in_bounds(last_press) && y_press != y) continue;
        const double zoom{center ? 1.0 : (in ? 0.1 : 10.0)};
        const double half{0.5 * range[y][2] * zoom};
        const double mid{icoord(y, last_press[y])};
        set_range(y, std::max(max_range[y][0], mid - half),
                  std::min(max_range[y][1], mid + half));
      }
    }
  }
  if (range != old_range || small_move) {
    small_move = false;
    prepare_draw();
  }
}

inline void X11Graph::leave(const XCrossingEvent &) {
  inside = false;
  if (destroyed) return;
  //  for (Radio * radio : radios) radio->erase();
  status = "";
  draw_controls();
}

inline void X11Graph::prepare_log() {
  // A one time operation to set up
  if (log_data.size() != input_data.size()) {
    log_data = Data(input_data.size());
    log_x_data = log_y_data = input_data;
    for (unsigned int s{0}; s != input_data.size(); ++s) {
      for (const bool y : {false, true}) {
        log_series.emplace_back(
            std::make_unique<Values>(input_data[s][y]->size()));
        Values & log_values{*log_series.back()};
        for (unsigned int p{0}; p != input_data[s][y]->size(); ++p)
          log_values[p] = log10((*input_data[s][y])[p]);
        log_data[s].push_back(&log_values);
      }
      log_x_data[s][0] = log_data[s][0];
      log_y_data[s][1] = log_data[s][1];
    }
  }

  // Select data view each time (precomputed to avoid log() delay)
  if (log_radios[0] && log_radios[1]) {
    data = &log_data;
  } else if (log_radios[0]) {
    data = &log_x_data;
  } else if (log_radios[1]) {
    data = &log_y_data;
  } else {
    data = &input_data;
  }
  get_range();
}

inline void X11Graph::set_clip_rectangle(
    const unsigned int x, const unsigned int y,
    const unsigned int width_, const unsigned int height_) {
  XRectangle clip_rectangle(rect(x, y, width_, height_));
  static XRectangle last_arc_clip_rectangle(rect(0, 0, 0, 0));
  if (clip_rectangle != last_arc_clip_rectangle) {
    for (const GC gc_ : series_arc_gcs) {
      XSetClipRectangles(display(), gc_, 0, 0, &clip_rectangle, 1, YXBanded);
      last_arc_clip_rectangle = clip_rectangle;
    }
  }
  static XRectangle last_line_clip_rectangle(rect(0, 0, 0, 0));
  if (clip_rectangle != last_line_clip_rectangle) {
    for (const GC gc_ : series_line_gcs)
      XSetClipRectangles(display(), gc_, 0, 0, &clip_rectangle, 1, YXBanded);
    last_line_clip_rectangle = clip_rectangle;
  }
}

inline void X11Graph::prepare() {
  drawn = false;
  // Set graph area and clip rectangle
  const int border{min_border()};
  set_bounds(border, extent(0) - border, border, extent(1) - border);
  set_clip_rectangle(bounds[0][0], bounds[1][0], bounds[0][2], bounds[1][2]);

  // Set scales
  for (const bool y : {false, true}) { scale[y] = bounds[y][2] / range[y][2]; }
  scale[2] = scale[1] / scale[0];

  // Determine arcs and points to connect with lines to display
  arcs.resize(data->size());
  points.resize(data->size());
  const Range erange{{range[0][0] - line_width, range[0][1] + line_width},
    {range[1][0] - line_width, range[1][1] + line_width}};

  auto series_fun = [this, &erange](const unsigned int s) {
    arcs[s].clear();
    points[s].clear();
    if (!series_radios[s]) return;

    // Arc properties
    const double radius{arc_radius};
    const double diam{radius * 2};
    XArc arc;  // FN
    arc.width = diam;
    arc.height = diam;
    arc.angle1 = 0;
    arc.angle2 = 64 * 360;

    const XYSeries & series{(*data)[s]};
    for (unsigned int p{0}; p != series[0]->size(); ++p) {
      const dPoint vals{(*series[0])[p], (*series[1])[p]};
      if (std::isfinite(vals[0]) && std::isfinite(vals[1])) {
        if (in_range(vals) && do_arcs(s)) {
          for (const bool y : {false, true})
            (y ? arc.y : arc.x) = (coord(y, vals[y])) - radius;
          arcs[s].emplace_back(arc);
        }
        // Complicated to handle lines exiting the display area properly
        if (do_lines(s)) {  // Assumes points ordered by X!!!
          if (vals[0] < erange[0][0]) continue;
          if (p) {
            // left - right range transition
            const dPoint last{(*series[0])[p - 1], (*series[1])[p - 1]};
            if (last.x > erange[0][1]) break;
            if (last.x < erange[0][0] || vals.x > erange[0][1]) {
              for (const bool left : {true, false}) {
                const double x_vertical{left ? erange[0][0] : erange[0][1]};
                if (last.x >= x_vertical || vals.x <= x_vertical) continue;
                const double y_int{line_vertical_y(last, vals, x_vertical)};
                if (y_int > erange[1][0] && y_int < erange[1][1])
                  points[s].push_back(xcoord(dPoint{x_vertical, y_int}));
              }
            }
            // top-bottom transition points, properly ordered by x
            if ((last.y < erange[1][0]) != (vals.y < erange[1][0]) ||
                (last.y < erange[1][1]) != (vals.y < erange[1][1])) {
              const bool last_low{last.y < erange[1][0]};
              for (const bool high :
                {last_low ? false : true, last_low ? true : false}) {
                const double y_horizontal{high ? erange[1][1] : erange[1][0]};
                if ((last.y > y_horizontal) == (vals.y > y_horizontal))
                  continue;
                const double x_int{line_horizontal_x(last, vals, y_horizontal)};
                if (x_int >= erange[0][0] && x_int <= erange[0][1])
                  points[s].push_back(xcoord(dPoint{x_int, y_horizontal}));
              }
            }
          }
          if (vals.x >= erange[0][0] && vals.x <= erange[0][1] &&
              vals.y >= erange[1][0] && vals.y <= erange[1][1])
            points[s].push_back(xcoord(vals));
        }
      }
    }
  };

  ThreadPool pool{n_threads()};
  std::vector<std::future<void>> futures;
  for (unsigned int s{0}; s != data->size(); ++s)
    futures.emplace_back(pool.run(series_fun, s));

  for (std::future<void> & result : futures) result.get();
}

inline void X11Graph::draw_grid() const {
  for (const bool y : {false, true}) {
    const Axis axis{range[y][0], range[y][1], 3, log_radios[y]};
    for (const std::pair<double, bool> tick : axis.ticks()) {
      if (!grid_radios[!tick.second][y]) continue;
      const int loc{coord(y, tick.first)};
      XDrawLine(display(), window, tick.second ? major_gc : minor_gc,
                y ? bounds[0][0] : loc, y ? loc : bounds[1][0],
                y ? bounds[0][1] : loc, y ? loc : bounds[1][1]);
    }
  }
}

inline void X11Graph::draw() {
  if (just_configured) {
    XFillRectangle(display(), pixmap, fill_gc, 0, 0, width(), height());
  } else {
    XFillRectangle(display(), pixmap, fill_gc,
                   bounds[0][0], bounds[1][0], bounds[0][2], bounds[1][2]);
  }

  // Handle special drawing commands
  Event pre_draw(Event::PreDraw);
  for (unsigned int c{0}; c != call_backs.size(); ++c)
    if (call_back_radios[c]) call_backs[c](*this, pre_draw);

  // The graph arcs and points to connect with lines
  const uint64_t arc_block{max_request / 3};
  const uint64_t line_block{max_request / 2};
  for (const unsigned int s : series_order) {
    if (do_arcs(s)) {
      for (unsigned int bs{0}; bs < arcs[s].size(); bs += arc_block) {
        const uint64_t n{bs + arc_block < arcs[s].size() ?
              arc_block : arcs[s].size() - bs};
        if (n) (outlines_radio ? XDrawArcs : XFillArcs)(
                display(), pixmap, series_arc_gcs[s], &arcs[s][bs],
                static_cast<unsigned int>(n));
      }
    }
    if (do_lines(s)) {
      for (unsigned int bs{0}; bs < points[s].size(); bs += line_block) {
        const uint64_t n{bs + line_block < points[s].size() ?
              line_block : points[s].size() - bs};
        if (n) XDrawLines(display(), pixmap, series_line_gcs[s],
                          &points[s][bs], static_cast<unsigned int>(n),
                          CoordModeOrigin);
      }
    }
  }
  if (true || just_configured) {
    just_configured = false;
    XCopyArea(display(), pixmap, window, gc, 0, 0, width(), height(), 0, 0);
    draw_controls();
  } else {
    redraw();
  }
  drawn = true;
  if (!small_move) {
    const SavedConfig current{current_config()};
    if (saved_config.empty() || current != saved_config.back()) {
      save_config(std::move(current));
      draw_controls();
    }
  }
  // std::cout << "Config stack size is " << saved_config.size() << std::endl;
  // show_range("After draw");
}

inline std::vector<Radio> X11Graph::create_unnamed_radios() {
  return std::vector<Radio>{
    {"Save an image of graph, and add all images so far to a pdf",
          this, {1, 1}, {[this]() { save_image("cn"); }}},
    {"Zoom out both axes", this, {1, -1}, {[this]() { get_range(2);
          prepare_draw(); }, [this]() { return zoomed[0] || zoomed[1]; }}},
    {"Zoom out X axis", this, {2, -1}, {[this]() {
          get_range(0); prepare_draw(); }, zoom_tester(0)}},
    {"Zoom out Y axis", this, {1, -2}, {[this]() {
          get_range(1); prepare_draw(); }, zoom_tester(1)}},
    {"Jump left X axis by one screen", this, {98.5, -1}, {[this]() {
          range_jump(0, -range[0][2]); prepare_draw(); }, zoom_tester(0)}},
    {"Jump left X axis by half a screen", this, {99.5, -1}, {[this]() {
          range_jump(0, -range[0][2] / 2); prepare_draw(); }, zoom_tester(0)}},
    {"Jump right X axis by half a screen", this, {100.5, -1}, {[this]() {
          range_jump(0, range[0][2] / 2); prepare_draw(); }, zoom_tester(0)}},
    {"Jump right X axis by one screen", this, {101.5, -1}, {[this]() {
          range_jump(0, range[0][2]); prepare_draw(); }, zoom_tester(0)}},
    {"Jump up Y axis by one screen", this, {1, 98.5}, {[this]() {
          range_jump(1, range[1][2]); prepare_draw(); }, zoom_tester(1)}},
    {"Jump up Y axis by half a screen", this, {1, 99.5}, {[this]() {
          range_jump(1, range[1][2] / 2); prepare_draw(); }, zoom_tester(1)}},
    {"Jump down Y axis by half a screen", this, {1, 100.5}, {[this]() {
          range_jump(1, -range[1][2] / 2); prepare_draw(); }, zoom_tester(1)}},
    {"Jump down Y axis by one screen", this, {1, 101.5}, {[this]() {
          range_jump(1, -range[1][2]); prepare_draw(); }, zoom_tester(1)}},
    {"Make markers bigger", this, {-1, -4.25}, {[this]() {
          arc_radius += 1; prepare_draw(); }, [this]() { return do_arcs(); }}},
    {"Make markers smaller", this, {-1, -3.25}, {[this]() {
          arc_radius -= 1; prepare_draw(); }, [this]() {
          return do_arcs() && arc_radius >= 2; }}},
    {"Make marker outlines thicker", this, {-1, -7.75}, {[this]() {
          set_line_widths(series_arc_gcs, arc_width += 1); draw(); }, [this]() {
          return do_arcs() && outlines_radio; }}},
    {"Make marker outlines thinner", this, {-1, -6.75}, {[this]() {
          set_line_widths(series_arc_gcs, arc_width -= 1); draw(); }, [this]() {
          return do_arcs() && outlines_radio && arc_width > 0; }}},
    {"Make series lines thicker", this, {-3.25, -1}, {[this]() {
          set_line_widths(series_line_gcs, line_width += 1); draw(); },
            [this]() { return do_lines(); }}},
    {"Make series lines thinner", this, {-4.25, -1}, {[this]() {
          line_width -= 1;
          set_line_widths(series_line_gcs, (line_width == 1 ? 0 : line_width));
          draw(); }, [this]() {return do_lines() && line_width >= 2; }}},
    {"Set default values for line and marker properties", this, {-1, -1},
      {[this]() {
          arcs_radio = true; outlines_radio = false; lines_radio = false;
          arc_radius = default_arc_radius; arc_width = default_arc_width;
          line_width = default_line_width;
          set_line_widths(series_arc_gcs, arc_width);
          set_line_widths(series_line_gcs, line_width);
          prepare_draw(); }, [this]() {
          return ((do_lines() &&
                   (lines_radio || dne(line_width, default_line_width))) ||
                  (do_arcs() && (!arcs_radio || outlines_radio ||
                                 dne(arc_radius, default_arc_radius) ||
                                 dne(arc_width, default_arc_width)))); }}}};
}

class X11TextGrid : public X11Win {
 public:
  using X11Win::X11Win;
  using Column = std::vector<std::string>;
  using Data = std::vector<Column>;
  using CellStatus = std::vector<std::vector<unsigned char>>;
  using CallBack = std::function<bool (const CellStatus &)>;

  static void create(X11App & app, const Data & data_,
                     const std::vector<unsigned int> & inactive_cols_ = {},
                     const std::vector<unsigned int> & inactive_rows_ = {},
                     const std::vector<unsigned int> & exclusive_cols_ = {},
                     const std::vector<unsigned int> & exclusive_rows_ = {},
                     CallBack call_back_ = CallBack(),
                     CallBack cell_test_ = CallBack(),
                     const unsigned int width__ = 1000,
                     const unsigned int height__ = 800,
                     const int x_off__ = 0,
                     const int y_off__ = 0) {
    app.add(std::make_unique<X11TextGrid>(
        app, data_, inactive_cols_, inactive_rows_,
        exclusive_cols_, exclusive_rows_,
        call_back_, cell_test_,
        width__, height__, x_off__, y_off__));
  }

  explicit X11TextGrid(X11App & app__, const Data & data_,
                       const std::vector<unsigned int> & inactive_cols_ = {},
                       const std::vector<unsigned int> & inactive_rows_ = {},
                       const std::vector<unsigned int> & exclusive_cols_ = {},
                       const std::vector<unsigned int> & exclusive_rows_ = {},
                       CallBack call_back_ = CallBack(),
                       CallBack cell_test_ = CallBack(),
                       const unsigned int width__ = 1000,
                       const unsigned int height__ = 800,
                       const int x_off__ = 0,
                       const int y_off__ = 0) :
      X11Win(app__, width__, height__, x_off__, y_off__, false),
      data{data_.size() ? data_ : Data(1, Column{"Empty"})},
    inactive_cols{inactive_cols_}, inactive_rows{inactive_rows_},
    exclusive_cols{exclusive_cols_}, exclusive_rows{exclusive_rows_},
    cell_status_(n_cols(), std::vector<unsigned char>(n_rows(), 0)),
    max_widths(n_cols()), call_back{call_back_}, cell_test{cell_test_} {
    // Events to watch out for
    XSelectInput(display(), window,
                 StructureNotifyMask | ExposureMask |
                 ButtonPressMask | ButtonReleaseMask);

    // Fonts of various size
    fonts.reserve(max_font_size);
    for (unsigned int s{3}; s <= max_font_size; ++s) {
      fonts.emplace_back(display(), s);
      if (fonts.back()) {
        font_sizes.push_back(s);
      } else {
        fonts.pop_back();
      }
    }
    if (fonts.empty()) throw Error("No fonts loaded");
    font = &fonts[fonts.size() / 2];

    XSync(display(), False);
    XColor grey;
    if (!XAllocNamedColor(display(), colormap, "rgb:cc/cc/cc",
                          &grey, &grey)) throw Error("Could not get grey");
    grey_gc = create_gc(grey.pixel, app.white);

    prepare();
    shrink_window_to_fit();
    XMapWindow(display(), window);
  }

  X11TextGrid(const X11TextGrid &) = delete;
  X11TextGrid & operator=(const X11TextGrid &) = delete;

  void shrink_window_to_fit() {
    XWindowChanges window_params;
    XResizeWindow(display(), window, layout_width(), layout_height());
    return;
    window_params.width = layout_width();
    window_params.height = layout_height();
    XConfigureWindow(display(), window, CWWidth | CWHeight, &window_params);
  }

  virtual void configure(const XConfigureEvent & event) {
    this->X11Win::configure(event);
    if (event.width > layout_width() || event.height > layout_height())
      shrink_window_to_fit();
  }
  virtual void button_press(const XButtonEvent & event) {
    last_motion = last_press = event;
    if (!in_bounds(event)) {
      for (Radio * radio : radios) { if (radio->press(event)) return; }
      return;
    }

    const Point cell{inside_cell(event)};
    if (count(inactive_rows.begin(), inactive_rows.end(), cell.y) ||
        count(inactive_cols.begin(), inactive_cols.end(), cell.x)) return;
    if (!cell_status(cell) &&
        count(exclusive_cols.begin(), exclusive_cols.end(), cell.x))
      set_column_status(cell.x, false);
    if (!cell_status(cell) &&
        count(exclusive_rows.begin(), exclusive_rows.end(), cell.y))
      set_row_status(cell.y, false);
    toggle_cell_status(cell);
    draw();
  }
  virtual void motion(const XMotionEvent &) { draw(); }
  virtual void button_release(const XButtonEvent &) {
    for (Radio * radio : radios) if (radio->release(last_press)) return;
  }

  template <class POINT>
  Point inside_cell(const POINT & point) const {
    return Point(static_cast<unsigned int>(upper_bound(
        column_offsets.begin(), column_offsets.end(),
        point.x) - column_offsets.begin() - 1),
                 (point.y - border_padding()) / cell_height());
  }
  bool cell_status(const Point cell) const {
    return cell_status_[cell.x][cell.y];
  }
  bool cell_status(const unsigned int x, const unsigned int y) const {
    return cell_status_[x][y];
  }
  void toggle_cell_status(const Point cell) {
    cell_status_[cell.x][cell.y] = !cell_status(cell);
  }
  void set_column_status(const unsigned int col, bool status) {
    cell_status_[col].assign(cell_status_[col].size(), status);
  }
  void set_row_status(const unsigned int row, bool status) {
    for (unsigned int col{0}; col != cell_status_.size(); ++col) {
      cell_status_[col][row] = status;
    }
  }
  void clear_status() {
    for (unsigned int c{0}; c != data.size(); ++c) set_column_status(c, false);
  }
  bool cells_selected() const {
    for (unsigned int c{0}; c != data.size(); ++c)
      for (unsigned int r{0}; r != data[c].size(); ++r)
        if (cell_status(c, r)) return true;
    return false;
  }
  unsigned int n_cells_selected() const {
    unsigned int result{0};
    for (unsigned int c{0}; c != data.size(); ++c)
      for (unsigned int r{0}; r != data[c].size(); ++r)
        if (cell_status(c, r)) ++result;
    return result;
  }

  int border_padding() const { return 50; }
  int cell_padding() const { return std::max(0.3 * font->height(), 10.0); }
  unsigned int n_rows() const { return data.empty() ? 0U :
        static_cast<unsigned int>(data[0].size()); }
  unsigned int n_cols() const { return static_cast<unsigned int>(data.size()); }
  int cell_width(const unsigned int col) const {
    return 2 * cell_padding() + max_widths[col];
  }
  int cell_height() const { return 2 * cell_padding() + font->height(); }
  int grid_width() const { return grid_width_; }
  int grid_height() const { return n_rows() * cell_height(); }
  int layout_width() const { return 2 * border_padding() + grid_width(); }
  int layout_height() const { return 2 * border_padding() + grid_height(); }
  unsigned int font_index() const {
    return static_cast<unsigned int>(font - &fonts[0]); }
  unsigned int font_size() const { return font_sizes[font_index()]; }
  int column_offset(const unsigned int col) const {
    return column_offsets[col];
  }
  int row_offset(const unsigned int row) const {
    return border_padding() + row * cell_height();
  }

  bool layout() {
    // Determines paddings, column and row sizes and sees if fits in window
    for (unsigned int c{0}; c != data.size(); ++c) {
      max_widths[c] = 0;
      for (unsigned int r{0}; r != data[c].size(); ++r)
        max_widths[c] = std::max(max_widths[c], font->string_width(data[c][r]));
    }

    const unsigned int total_text_max_width{
      std::accumulate(max_widths.begin(), max_widths.end(), 0U)};

    grid_width_ = total_text_max_width + 2 * n_cols() * cell_padding();

    if (layout_width() > width()) return false;
    if (layout_height() > height()) return false;
    if (0)
      std::cerr << "Layout "
                << n_cols() << " " << n_rows() << " "
                << width() << " " << height() << " "
                << grid_width() << " " << grid_height() << " "
                << total_text_max_width << " " << font_size() << std::endl;
    return true;
  }

  virtual void prepare() {
    // Get optimal layout to just fit
    if (0) {
      while (layout()) {
        if (font == &fonts.back()) break;
        ++font;
      }
    }
    font = &fonts.back();
    while (!layout()) {
      if (font == &fonts.front()) break;
      --font;
    }

    int column_offset_{border_padding()};
    column_offsets.clear();
    for (unsigned int c{0}; c != data.size(); ++c) {
      column_offsets.push_back(column_offset_);
      column_offset_ += cell_width(c);
    }
    column_offsets.push_back(column_offset_);

    XSetFont(display(), gc, font->id());
    XSetFont(display(), fill_gc, font->id());

    set_bounds(border_padding(), border_padding() + grid_width(),
               border_padding(), border_padding() + grid_height());
  }

  virtual void draw() {
    // White page
    XFillRectangle(display(), window, fill_gc, 0, 0, width(), height());
    just_configured = false;

    std::vector<XRectangle> rectangles{
      rect(border_padding(), border_padding(), grid_width(), grid_height())};
    std::vector<XRectangle> fill_rectangles;
    for (unsigned int c{0}; c != data.size(); ++c) {
      for (unsigned int r{0}; r != data[c].size(); ++r) {
        XRectangle cell_rectangle(rect(column_offset(c), row_offset(r),
                                       cell_width(c), cell_height()));
        rectangles.push_back(cell_rectangle);
        cell_rectangle.x += 1;  // FN
        cell_rectangle.y += 1;
        cell_rectangle.width -= 1;
        cell_rectangle.height -= 1;
        if (cell_status(c, r)) fill_rectangles.push_back(cell_rectangle);
      }
    }
    XDrawRectangles(display(), window, gc, &rectangles[0],
                    static_cast<unsigned int>(rectangles.size()));
    if (fill_rectangles.size())
      XFillRectangles(display(), window, grey_gc, &fill_rectangles[0],
                      static_cast<unsigned int>(fill_rectangles.size()));
    for (unsigned int c{0}; c != data.size(); ++c)
      for (unsigned int r{0}; r != data[c].size(); ++r)
        XDrawString(display(), window, gc, column_offset(c) + cell_padding(),
                    row_offset(r) + cell_padding() + font->height(),
                    data[c][r].c_str(),
                    static_cast<unsigned int>(data[c][r].size()));
    for (const Radio * radio : radios) { radio->draw(); }
  }

  virtual ~X11TextGrid() { }

  Data data{};
  std::vector<unsigned int> inactive_cols{};
  std::vector<unsigned int> inactive_rows{};
  std::vector<unsigned int> exclusive_cols{};
  std::vector<unsigned int> exclusive_rows{};
  CellStatus cell_status_{};

  GC grey_gc{};

  static constexpr unsigned int max_font_size{60};
  std::vector<unsigned int> font_sizes{};
  std::vector<X11Font> fonts{};
  X11Font * font{};

  double border_padding_factor{0.0};
  int border_width{3};
  int cell_border_width{2};
  int grid_width_{0};
  std::vector<int> max_widths, column_offsets{};
  CallBack call_back, cell_test{};

  Point last_press{}, last_motion{};

  Radio bigger_radio{"Bigger_text", this, {1, 98.5},
    {[this]() { }, [this]() {
        return font_index() + 1 != fonts.size(); }, [this]() {
        ++font; layout(); shrink_window_to_fit(); prepare_draw(); }}};
  Radio smaller_radio{"Bigger_text", this, {1, 99.5},
    {[this]() { }, [this]() {
        return font_index() != 0; }, [this]() {
        --font; layout(); shrink_window_to_fit(); prepare_draw(); }}};
  Radio clear_radio{"Clear all selections", this, {1, 1},
    {[this]() { clear_status(); draw(); },
          [this]() { return cells_selected(); }}};
  Radio plot_radio{"Plot selected data", this, {-1, 1},
    {[this]() { call_back(cell_status_); clear_status(); draw(); },
          [this]() { return call_back && cell_test &&
                cell_test(cell_status_); }}};
  std::vector<Radio *> radios{&clear_radio, &plot_radio};
};

class X11Plotter {
  using Names = std::vector<std::string>;
  using Col = std::vector<double>;
  using Data = std::vector<Col>;

 public:
  template <class TSV>
  explicit X11Plotter(const TSV & tsv) {
    std::vector<std::vector<std::string>> text{
      {"Data Field"}, {"Plot X"}, {"Plot Y"}};
    for (unsigned int c{0}; c != tsv.n_cols(); ++c) {
      const auto & col = tsv(c);
      if (col.is_real()) {
        names.push_back(col.name());
        text[0].push_back(col.name());
        text[1].push_back("");
        text[2].push_back("");
        data.emplace_back();
        for (unsigned int r{0}; r != tsv.n_rows(); ++r)
          if (tsv[c].is_integral()) {
            data.back().push_back(tsv.as_jitter(c, r));
          } else {
            data.back().push_back(tsv.as_real(c, r));
          }
      }
    }
    X11TextGrid::CallBack call_back{std::bind(&X11Plotter::launch_graph, this,
                                              std::placeholders::_1)};
    X11TextGrid::CallBack cell_test{std::bind(&X11Plotter::launch_ready, this,
                                              std::placeholders::_1)};
    X11TextGrid::create(app, text, {0}, {0}, {1}, {}, call_back, cell_test);
    app.run();
  }

  bool launch_graph(const X11TextGrid::CellStatus & status) {
    X11Graph::Data gd;
    X11Graph::Values * xs{&data[std::find(
        status[1].begin(), status[1].end(), 1) - status[1].begin() - 1]};
    for (unsigned int n{0}; n != names.size(); ++n)
      if (status[2][n + 1]) gd.emplace_back(X11Graph::XYSeries{xs, &data[n]});
    X11Graph & graph{X11Graph::create_whole(app, gd, 1300, 850)};
    graph.arc_radius = 1;
    return true;
  }
  bool launch_ready(const X11TextGrid::CellStatus & status) const {
    return std::find(status[1].begin(), status[1].end(), 1) != status[1].end()
        && std::find(status[2].begin(), status[2].end(), 1) != status[2].end();
  }

  X11App app{};
  Data data{};
  Names names{};
};



#pragma GCC diagnostic pop

}  // namespace paa

#endif
