//
// psplot.h
//
// postscript plots
//
// Copyright 2016 Peter Andrews @ CSHL
//

#ifndef PAA_PSPLOT_H
#define PAA_PSPLOT_H

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <functional>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
#include "error.h"
#include "files.h"
#include "plot.h"
#include "stats.h"
#include "utility.h"

#define VERBOSE 0

namespace paa {

class DocSettings {
 public:
  DocSettings() { }
  explicit DocSettings(const double padding__) : padding_{padding__} { }
  DocSettings(const unsigned int width__,
              const unsigned int height__,
              const double padding__ = default_doc_padding) :
      width_{width__},
    height_{height__},
    padding_{padding__} { }

  DocSettings & reset() { return *this = DocSettings(); }
  bool eps() const { return eps_; }
  bool & eps() { return eps_; }
  DocSettings & eps(const bool eps__) {
    eps_ = eps__;
    return *this;
  }
  bool ps() const { return ps_; }
  bool & ps() { return ps_; }
  DocSettings & ps(const bool ps__) {
    ps_ = ps__;
    return *this;
  }
  bool pdf() const { return pdf_; }
  bool & pdf() { return pdf_; }
  DocSettings & pdf(const bool pdf__) {
    pdf_ = pdf__;
    return *this;
  }
  bool png() const { return png_; }
  DocSettings & png(const bool png__) {
    png_ = png__;
    return *this;
  }
  void pdf_only() {
    pdf_ = true;
    ps_ = false;
    eps_ = false;
    png_ = false;
  }
  bool ticks() const { return ticks_; }
  bool & ticks() { return ticks_; }
  DocSettings & ticks(const bool ticks__) {
    ticks_ = ticks__;
    return *this;
  }
  unsigned int width() const { return width_; }
  unsigned int & width() { return width_; }
  DocSettings & width(const unsigned int width__) {
    width_ = width__;
    return *this;
  }
  unsigned int height() const { return height_; }
  unsigned int & height() { return height_; }
  DocSettings & height(const unsigned int height__) {
    height_ = height__;
    return *this;
  }
  double padding() const { return padding_; }
  double & padding() { return padding_; }
  DocSettings & padding(const double padding__) {
    padding_ = padding__;
    return *this;
  }

 private:
  bool eps_{false};
  bool ps_{true};
  bool pdf_{false};
  bool png_{false};
  bool ticks_{true};
  unsigned int width_{default_doc_width};
  unsigned int height_{default_doc_height};
  double padding_{default_doc_padding};
};
// Move to cpp file if linking problems ensue
DocSettings doc_defaults;

inline std::string date_time(std::time_t t = std::time(nullptr)) {
  char output[100];
  //  if (std::strftime(output, sizeof(output), "%c", std::localtime(&t))) {
  if (std::strftime(output, sizeof(output), "%a %b %d %T %Y",
                    std::localtime(&t))) {
    return std::string{output};
  } else {
    return "bad time";
  }
}

// Layout of a page
class Layout {
 public:
  // spaces must be exact. ^ means any character, repeated
  // Number X [optional x fractions]
  // Number Y [optional y fractions]
  // Sub-formats: space X index Y index CHAR format CHAR
  // X [(xf1 xf2 ... xfN-1X)] Y [(yf1 yf2 ... yfN-1Y)][ ^x y layout^]...
  explicit Layout(Bounds bounds__, const std::string format = "",
                  double padding = 0.01) :
      bounds_{bounds__} {
    if (format != "") {
      unsigned int nxy[2];
      std::vector<double> fractions[2];
      double paddings[2];
      std::vector<std::vector<std::string>> formats;
      std::istringstream format_stream{format.c_str()};
      for (const bool y : {false, true}) {
        char c;
        format_stream >> nxy[y];
        format_stream.get();
        paddings[y] = padding * bounds_.w(y);
        double total{(nxy[y] - 1) * padding};
        const double fcorr{1 - total};
        const bool explicit_fractions{format_stream.peek() == '('};
        if (explicit_fractions) {
          format_stream.get();
          for (unsigned int x{0}; x + 1 != nxy[y]; ++x) {
            double fraction;
            format_stream >> fraction;
            fraction *= fcorr;
            total += fraction;
            fractions[y].push_back(fraction);
          }
          if (total >= 1)
            throw Error("Total layout fractions too big") << format;
          fractions[y].push_back(1 - total);
          format_stream >> c;
          if (c != ')') {
            std::string rest;
            getline(format_stream, rest);
            throw Error("Layout expected close parentheses")
                << "in" << format << " at " << rest;
          }
        } else {
          for (unsigned int x{0}; x != nxy[y]; ++x) {
            fractions[y].push_back(1.0 / nxy[y] * fcorr);
          }
        }
      }

      // Read in sub-formats
      formats.resize(nxy[0], std::vector<std::string>(nxy[1]));
      char c;
      unsigned int sx;
      unsigned int sy;
      while (format_stream >> c >> sx >> sy) {
        format_stream.get();
        getline(format_stream, formats[sx][sy], c);
      }

      // Set up sub components
      double x_start{bounds_.xl()};
      for (unsigned int x{0}; x != nxy[0]; ++x) {
        components.push_back(std::vector<Layout>());
        double x_stop{x_start + bounds_.xw() * fractions[0][x]};
        double y_stop{bounds_.yh()};
        for (unsigned int y{0}; y != nxy[1]; ++y) {
          double y_start{y_stop - bounds_.yw() * fractions[1][y]};
          components.back().emplace_back(
              Bounds{x_start, x_stop, y_start, y_stop},
              formats[x][y]);
          y_stop = y_start - paddings[1];
        }
        x_start = x_stop + paddings[0];
      }
    }
  }

  void describe(const unsigned int level = 0) const {
    std::cout << std::string(level * 2, ' ');
    sout << nx() << ny() << n() << bounds_ << std::endl;
    if (nx() != 1 || ny() != 1) {
      for (unsigned int x{0}; x != nx(); ++x) {
        for (unsigned int y{0}; y != ny(); ++y) {
          std::cout << std::string(level * 2 + 1, ' ');
          sout << x << y << std::endl;
          components[x][y].describe(level + 1);
        }
      }
    }
  }

  unsigned int nx() const {
    if (components.empty()) {
      return 1;
    } else {
      return static_cast<unsigned int>(components.size());
    }
  }
  unsigned int ny() const {
    if (components.empty()) {
      return 1;
    } else {
      return static_cast<unsigned int>(components.front().size());
    }
  }
  unsigned int n() const {
    if (components.empty()) return 1;
    unsigned int result{0};
    for (unsigned int x{0}; x != nx(); ++x) {
      for (unsigned int y{0}; y != ny(); ++y) {
        result += components[x][y].n();
      }
    }
    return result;
  }
  const Layout & operator[](const unsigned int i) const {
    unsigned int total{0};
    if (components.empty()) {
      if (i) throw Error("Unexpected i in Layout::operator[]") << i;
      return *this;
    }
    for (unsigned int x{0}; x != nx(); ++x) {
      for (unsigned int y{0}; y != ny(); ++y) {
        const unsigned int this_n{components[x][y].n()};
        if (i < total + this_n) {
          return components[x][y][i - total];
        }
        total += this_n;
      }
    }
    throw Error("Unexpected i in Layout::operator[] 2") << i;
  }

  const Bounds & bounds() const { return bounds_; }
  template <class ... ARGS>
  const Bounds & bounds(const unsigned int x, const unsigned int y,
                        ARGS ... rest) const {
    if (components.empty()) {
      if (x || y) throw Error("Unexpected index in Layout::bounds") << x << y;
      return bounds();
    }
    return components[x][y].bounds(rest ...);
  }

 private:
  Bounds bounds_;
  std::vector<std::vector<Layout>> components{};
};

// Connects parent, self and child classes to each other
template <class Parent, class Self, class Child>
class Multiplexer {
 public:
  // Construct
  explicit Multiplexer(Self * self__) : self_{self__} { }
  template <class ... Parents>
  Multiplexer(Self * self__, Parents ... parents__) : self_{self__} {
    add(parents__ ...);
  }
  Multiplexer(Multiplexer && other) = default;

  Multiplexer(const Multiplexer &) = delete;
  Multiplexer & operator=(const Multiplexer &) = delete;

  void erase_parent(void * parent_) {
    typename std::vector<Parent *>::iterator found{find(
        parents_.begin(), parents_.end(), reinterpret_cast<Parent *>(parent_))};
    if (found != parents_.end()) {
      parents_.erase(found);
    }
  }
  void erase_child(void * child) {
    typename std::vector<Child *>::iterator found{find(
        children_.begin(), children_.end(), reinterpret_cast<Child *>(child))};
    if (found != children_.end()) {
      children_.erase(found);
    }
  }

  // Add nodes
  Multiplexer & add(Parent * parent_) {
    if (parent_) {
      parents_.push_back(parent_);
      parent_->children().push_back(self_);
    }
    return *this;
  }
  Multiplexer & add(Child * child) {
    if (child) {
      children_.push_back(child);
      child->parents().push_back(self_);
    }
    return *this;
  }
  Multiplexer & add(std::initializer_list<Parent *> ps) {
    for (Parent * p : ps) {
      add(p);
    }
    return *this;
  }
  Multiplexer & add(std::initializer_list<Child *> cs) {
    for (Child * c : cs) {
      add(c);
    }
    return *this;
  }
  Multiplexer & manage(std::unique_ptr<Parent> && m) {
    add(m.get());
    managed_parents_.push_back(std::move(m));
    return *this;
  }
  Multiplexer & manage(std::unique_ptr<Child> && m) {
    add(m.get());
    managed_children_.push_back(std::move(m));
    return *this;
  }
  template <class PTR>
  Multiplexer & ownp(PTR * m) {
    parents().front()->managed_children_.emplace_back(m);
    return *this;
  }

  // Access
  const std::vector<Parent *> & parents() const { return parents_; }
  const std::vector<Child *> & children() const { return children_; }
  std::vector<Parent *> & parents() { return parents_; }
  std::vector<Child *> & children() { return children_; }
  Parent & parent() { return *parents_.front(); }

  // Finalize parents
  void finalize_doc() {
    for (Parent * parent_ : parents_) {
      parent_->finalize_doc();
    }
  }

 protected:
  ~Multiplexer() {
    for (Parent * parent_ : parents_) {
      parent_->erase_child(self_);
    }
    for (Child * child : children_) {
      child->erase_parent(self_);
    }
  }

 private:
  Self * self_{};
  std::vector<Parent *> parents_{};
  std::vector<Child *> children_{};

 public:
  std::vector<std::unique_ptr<Parent>> managed_parents_{};
  std::vector<std::unique_ptr<Child>> managed_children_{};
};

template <class PTR>
void ownp(const PTR ptr) {
  if (ptr) ptr->ownp(ptr);
}

// Postscript document
template<class PSPage>
class PSDocT : public DocSettings,
               public Multiplexer<PSDocT<PSPage>, PSDocT<PSPage>, PSPage> {
 public:
  using Multi = Multiplexer<PSDocT, PSDocT, PSPage>;
  using Multi::children;

  // Construct
  PSDocT(const std::string & filename,
         const std::string & title_ = "",
         const unsigned int width__ = doc_defaults.width(),
         const unsigned int height__ = doc_defaults.height(),
         const double padding__ = doc_defaults.padding()) :
      // DocSettings{width__, height__, padding__},
      DocSettings{doc_defaults},
    Multi{this},
    file_name{filename},
    title{title_.size() ? title_ : file_name} {
      width(width__);
      height(height__);
      padding(padding__);
    }

  // Destroy
  ~PSDocT() {
    if (VERBOSE) std::cerr << "Destroy doc " << title << std::endl;
    if (!children().empty()) finalize_doc();
  }

  std::string ext() const { return eps() ? ".eps" : ".ps"; }

  // Finalize
  void output() { finalize_doc(); }
  void finalize_doc() try {
    if (finalized) return;
    finalized = true;
    if (eps() && children().size() != 1) {
      std::cerr << "Must have one page in an eps file: " << children().size()
                << std::endl;
      return;
    }

    // Open output file
    psout.open((file_name + ext() + ".tmp").c_str());
    if (!psout) std::cerr << "Problem opening ps file " << file_name
                          << std::endl;

    // Write postscript header
    psout << std::setprecision(default_precision);
    if (eps()) {
      psout << "%!PS-Adobe-3.0 EPSF-3.0\n";
    } else {
      psout << "%!PS-Adobe-3.0\n";
    }
    psout << "%%Pages: " << children().size() << "\n";
    psout << "%%BoundingBox: 0 0 " << width() << " " << height()
          << R"foo(
%%Title: )foo"
          << title
          << R"foo(
%%Creator: Peter Andrews
%%CreationDate: )foo" << date_time() << R"foo(
%%EndComments
%%BeginProlog
% stack
/e {exch} bind def
% graphics
/gs {gsave} bind def
/gr {grestore} bind def
/tr {translate} bind def
/ra {rotate} bind def
/sc {dup scale} bind def
/cxy {currentpoint} bind def
/ps {showpage} bind def
% movement
/m {moveto} bind def
/rm {rmoveto} bind def
% paths
/np {newpath} bind def
/cp {closepath} bind def
/sp {stroke} bind def
/fp {fill} bind def
/l {lineto} bind def
/rl {rlineto} bind def
/cs {gsave sstroke grestore clip} bind def
/ic {initclip} bind def
% line properties
/lw {setlinewidth} bind def
/sd {setdash} bind def
/nd {[] 0 setdash} bind def
% set color
/c {setrgbcolor} bind def
/bk {0 0 0} def
% text
/s {show} bind def
/sf {/Helvetica findfont e scalefont setfont} bind def 
% justification
/sw {stringwidth} def
/jcx {dup sw pop 2 div neg} bind def
/jc {jcx 0 rm} bind def
/jrx {dup sw pop neg} bind def
/jr {jrx 0 rm} bind def
% shortening
/nrs {neg rm show} bind def
/cbb {m lineto lineto lineto closepath clip} bind def
/bb {setlinewidth bk setrgbcolor nd newpath m lineto lineto lineto
     closepath stroke} bind def
/bbf {setlinewidth bk setrgbcolor nd newpath m lineto lineto lineto
      closepath gssave setrgbcolor fill grestore stroke} bind def
%%EndProlog
)foo";

    const Bounds bounds(width(), height(), padding());
    for (unsigned int p{0}; p != children().size(); ++p) {
      PSPage * page{children()[p]};
      page->finalize(*this, p + 1, bounds);
    }
    psout << "%%Trailer\n";
    psout << "%%EOF\n";
    psout.close();

    std::ostringstream mv;
    mv << "mv " << file_name << ext() << ".tmp " << file_name << ext();
    if (system(mv.str().c_str()) == -1) {
      std::cerr << "Problem moving ps file" << std::endl;
    }
    if (pdf() || png()) {
      std::ostringstream ps2pdf;
      ps2pdf << "ps2pdf -dAutoRotatePages=/None -dDEVICEWIDTHPOINTS=" << width()
             << " -dDEVICEHEIGHTPOINTS=" << height()
             << " " << file_name << ext() << " " << file_name << ".pdf";
      if (system(ps2pdf.str().c_str()) == -1) {
        std::cerr << "Problem creating pdf file" << std::endl;
      }
      if (png()) {
        std::ostringstream pdf2png;
        pdf2png << "convert " << file_name << ".pdf " << file_name << ".png";
        if (system(pdf2png.str().c_str()) == -1) {
          std::cerr << "Problem creating png file" << std::endl;
        }
      }
    }
    if (!ps()) unlink(file_name + ext());
    if (!pdf() && png()) unlink(file_name + ".pdf");
  } catch (Error & e) {
    std::cerr << "paa::Error:\n";
    std::cerr << e.what() << std::endl;
  } catch(std::exception & e) {
    std::cerr << e.what() << std::endl;
  } catch(...) {
    std::cerr << "Some exception was caught in finalize_doc."<< std::endl;
  }

  // Doc info
  unsigned int n_pages() const { return children().size(); }

  // Write to postsctipt file
  template <class Type>
  PSDocT & operator<<(const Type & val) {
    psout << val;
    return *this;
  }
  PSDocT & operator<<(std::ostream & (*)(std::ostream &)) {  // for endl
    psout << std::endl;
    return *this;
  }
  void annotate(const std::string & text) {
    for (unsigned int p{0}; p != children().size(); ++p)
      children()[p]->annotate(text, width() - 12, height() - 18.5);
  }

 private:
  std::string file_name{};
  std::string title{};
  std::ofstream psout{};
  bool finalized{false};
};

// Text justification
enum class Just {
  left, center, right
};

class Text {
 public:
  // Constructors
  static constexpr double default_size{12};
  static constexpr Just default_just{Just::left};
  // All arguments any order
  Text(const std::string & text__,
       const double size__,
       const Just just__,
       const std::string & color__ = "bk") :
      text_{text__}, size_{size__}, just_{just__}, color_{color__} { }
  Text(const std::string & text__,
       const double size__,
       const std::string & color__,
       const Just just__ = default_just) :
      Text{text__, size__, just__, color__} {}
  Text(const std::string & text__,
       const std::string & color__,
       const double size__,
       const Just just__ = default_just) :
      Text{text__, size__, just__, color__} {}
  Text(const std::string & text__,
       const std::string & color__,
       const Just just__,
       const double size__ = default_size) :
      Text{text__, size__, just__, color__} {}
  Text(const std::string & text__,
       const Just just__,
       const std::string & color__,
       const double size__ = default_size) :
      Text{text__, size__, just__, color__} {}
  Text(const std::string & text__,
       const Just just__,
       const double size__,
       const std::string & color__ = "bk") :
      Text{text__, size__, just__, color__} {}
  // Two arguments any order
  Text(const std::string & text__,
       const double size__) :
      Text{text__, size__, default_just, "bk"} {}
  Text(const std::string & text__,
       const Just just__) :
      Text{text__, default_size, just__, "bk"} {}
  Text(const std::string & text__,
       const std::string & color__) :
      Text{text__, default_size, default_just, color__} {}
  // One argument
  Text(const std::string & text__) :  // NOLINT
      Text{text__, default_size, default_just, "bk"} {}
  Text(const char * text__) :  // NOLINT
      Text{text__, default_size, default_just, "bk"} {}
  // Default constructor
  Text() : Text{"", default_size, default_just, "bk"} {}

  // Access
  operator bool() const { return text_.size(); }
  std::string text() const { return text_; }
  double size() const { return size_; }
  std::string color() const { return color_; }
  Just just() const { return just_; }

  // Modification
  std::string & text() { return text_; }
  double & size() { return size_; }
  std::string & color() { return color_; }
  Just & just() { return just_; }

  // Chained modification
  Text & text(const std::string & text__) {
    text_ = text__;
    return *this;
  }
  Text & size(const double & size__) {
    size_ = size__;
    return *this;
  }
  Text & color(const std::string & color__) {
    color_ = color__;
    return *this;
  }
  Text & just(const Just & just__) {
    just_ = just__;
    return *this;
  }

  // Render, passing in optional pre-existing size and color
  std::string ps(const double scale = 1.0,
                 const double doc_size = 0.0,
                 const std::string & doc_color = "") const {
    const double size__{size_ * scale};
    std::ostringstream result;
    result << std::setprecision(default_precision);
    if (size_ < doc_size || size_ > doc_size) {
      result << size__ << " sf ";
    }
    if (color_ != doc_color) {
      result << color_ << " c ";
    }
    result << "(" << text_ << ") ";
    if (just_ == Just::center) {
      result << "jc ";
    } else if (just_ == Just::right) {
      result << "jr ";
    }
    result << "s ";
    return result.str();
  }

 private:
  std::string text_;
  double size_;
  Just just_;
  std::string color_;
};

// Postscript page
template<class PSGraph>
class PSPageT :
      public Multiplexer<PSDocT<PSPageT<PSGraph> >, PSPageT<PSGraph>, PSGraph> {
 public:
  using PSPage = PSPageT<PSGraph>;
  using PSDoc = PSDocT<PSPage>;
  using Multi = Multiplexer<PSDoc, PSPage, PSGraph>;
  using Multi::add;
  using Multi::children;
  using Multi::finalize_doc;
  using Multi::manage;

  // factory
  template <class ... Args>
  static PSPageT * create(Args && ... args) {
    PSPageT * page{new PSPageT{std::forward<Args>(args)...}};
    ownp(page);
    return page;
  }

  // Construct
  explicit PSPageT(const std::string & title__ = "",
                   const std::string & layout__ = "") :
      Multi{this}, title_{title__, Just::center, 22}, layout_{layout__} { }
  explicit PSPageT(PSDoc & doc,
                   const std::string & title__ = "",
                   const std::string & layout__ = "") :
      PSPageT{title__, layout__} {
    add(&doc);
  }
  explicit PSPageT(const std::string & file_name,
                   const std::string & title__,
                   const std::string & layout__) :
      PSPageT{title__, layout__} {
    std::unique_ptr<PSDoc> doc_{std::make_unique<PSDoc>(file_name, title__)};
    manage(std::move(doc_));
  }
  virtual ~PSPageT() {
    if (VERBOSE) std::cerr << "Destroy page " << title_.text() << std::endl;
    finalize_doc();
  }
  PSPageT(PSPageT &&) = default;
  const Text & title() const { return title_; }
  Text & title() { return title_; }
  PSPageT & title(const Text & title__) {
    title_ = title__;
    return *this;
  }

  // Finalize
  virtual void finalize(PSDoc & doc, const unsigned int page, Bounds bounds) {
    // Start page
    doc << "%%Page: " << page << " " << page << "\n";
    doc << "save\n";

    // Title
    if (children().size() && title_) {
      doc << bounds.xc() << " " << bounds.yh() - 0.7 * title_.size()
          << " m " << title_.ps() << "\n";
      bounds.yh(bounds.yh() - 1.2 * title_.size());
    }

    // Adjust bounds and set layout
    const Layout layout{bounds, layout_};

    // Finalize graphs
    if (children().empty()) {
      doc << title_.text() << "\n";
    } else {
      for (unsigned int c{0}; c != children().size(); ++c) {
        children()[c]->finalize(doc, layout[c].bounds());
      }
    }

    // Add ps annotation
    if (ps_.size()) doc << ps_ << '\n';

    // Finish page
    doc << "restore\n";
    doc << "ps\n" << "%%EndPage: " << page << std::endl;
  }
  void annotate(const std::string & text, const double x, const double y) {
    std::ostringstream ps;
    ps << "10 sf 0 0 0 c " << x << " " << y << " m "
       << "(" << text << ") jr s";
    ps_ = ps.str();
  }

 private:
  Text title_;
  std::string layout_;
  std::string ps_{};
};

// Get components of labels from string like "Title;X Label;Y Label"
inline void parse_title(const std::string & input, std::string & title,
                        std::string & x_label, std::string & y_label) {
  std::istringstream in{input.c_str()};
  getline(in, title, ';');
  getline(in, x_label, ';');
  getline(in, y_label, ';');
}


template<class PSSeries>
class PSPartT : public Multiplexer<PSPageT<PSPartT<PSSeries> >,
                                   PSPartT<PSSeries>, PSSeries> {
 public:
  using PSPart = PSPartT<PSSeries>;
  using PSPage = PSPageT<PSPart>;
  using PSDoc = PSDocT<PSPage>;
  using Multi = Multiplexer<PSPage, PSPart, PSSeries>;

  PSPartT() : Multi{this} { }
  PSPartT(PSPartT &&) = default;
  virtual ~PSPartT() { }
  virtual void finalize(PSDoc & doc, Bounds bounds) = 0;
  PSPart & add_text(const double xf__, const double yf__, const Text & text__) {
    text_.emplace_back(std::make_pair(xf__, yf__), text__);
    return *this;
  }
  void write_text(PSDoc & doc, const Bounds & bounds, const double scale) {
    for (const std::pair<std::pair<double, double>, Text> & text_pos : text_) {
      const double x{text_pos.first.first};
      const double y{text_pos.first.second};
      Text text{text_pos.second};
      text.size() *= scale;
      doc << x * bounds.xw() + bounds.xl() << " "
          << y * bounds.yw() + bounds.yl() << " m "
          << text.ps() << "\n";
    }
  }

  bool log_x() const { return log_x_; }
  PSPartT & log_x(const bool log_x__) {
    log_x_ = log_x__;
    return *this;
  }
  bool log_y() const { return log_y_; }
  PSPartT & log_y(const bool log_y__) {
    log_y_ = log_y__;
    return *this;
  }
  PSPartT & log_log(const bool log__) {
    log_x_ = log__;
    log_y_ = log__;
    return *this;
  }

 protected:
  bool log_x_{false};
  bool log_y_{false};
  std::vector<std::pair<std::pair<double, double>, Text>> text_{};
};

template<class PSSeries>
class PSGraphT : public GraphSettings, public PSPartT<PSSeries> {
 public:
  using PSPart = PSPartT<PSSeries>;
  using PSPage = PSPageT<PSPart>;
  using PSDoc = PSDocT<PSPage>;
  using PSPart::Multi::add;
  using PSPart::Multi::children;
  using PSPart::Multi::finalize_doc;
  using PSPart::Multi::manage;
  using PSPart::write_text;

  // factory
  template <class ... Args>
  static PSGraphT * create(Args && ... args) {
    PSGraphT * graph__{new PSGraphT{std::forward<Args>(args)...}};
    ownp(graph__);
    return graph__;
  }

  // Construct
  explicit PSGraphT(const std::string & title__ = "",
                    const Bounds & range__ = Bounds{0.0}) :
      GraphSettings{graph_defaults},
      PSPart{},
    range_{range__} {
      std::string title___;
      std::string labels__[2];
      parse_title(title__, title___, labels__[0], labels__[1]);
      title_ = Text{title___, Just::center, graph_defaults.title_size()};
      labels_[0] = Text{labels__[0], Just::right, graph_defaults.label_size()};
      labels_[1] = Text{labels__[1], Just::right, graph_defaults.label_size()};
    }
  explicit PSGraphT(PSPage & page__) :
      PSGraphT{} {
    add(&page__);
  }
  explicit PSGraphT(PSPage & page__,
                    const std::string & title__,
                    const Bounds & range__ = Bounds{0.0}) :
      PSGraphT{title__, range__} {
    add(&page__);
  }
  explicit PSGraphT(PSPage & page__,
                    const Bounds & range__,
                    const std::string & title__ = "") :
      PSGraphT{title__, range__} {
    add(&page__);
  }
  explicit PSGraphT(PSDoc & doc) :
      PSGraphT{} {
    std::unique_ptr<PSPage> page_{std::make_unique<PSPage>(doc, "")};
    manage(std::move(page_));
  }
  PSGraphT(PSDoc & doc,
           const std::string & title__,
           const Bounds & range__ = Bounds{0.0}) :
      PSGraphT{title__, range__} {
    std::unique_ptr<PSPage> page_{std::make_unique<PSPage>(doc, "")};
    manage(std::move(page_));
  }
  PSGraphT(PSDoc & doc,
           const Bounds & range__,
           const std::string & title__ = "") :
      PSGraphT{title__, range__} {
    std::unique_ptr<PSPage> page_{std::make_unique<PSPage>(doc, "")};
    manage(std::move(page_));
  }
  PSGraphT(const std::string & file_name,
           const std::string & title__,
           const Bounds & range__ = Bounds{0.0}) :
      PSGraphT{title__, range__} {
    std::unique_ptr<PSPage> page_{std::make_unique<PSPage>(file_name, "", "")};
    manage(std::move(page_));
  }
  virtual ~PSGraphT() {
    if (VERBOSE) std::cerr << "Destroy graph " << title_.text() << std::endl;
    finalize_doc();
  }
  PSGraphT(PSGraphT &&) = default;

  // Graph info
  const Bounds & range() const { return range_; }
  const Text & title() const { return title_; }
  const Text & x_label() const { return labels_[0]; }
  const Text & y_label() const { return labels_[1]; }
  Bounds & range() { return range_; }
  Text & title() { return title_; }
  Text & x_label() { return labels_[0]; }
  Text & y_label() { return labels_[1]; }
  PSGraphT & range(const Bounds & range__) {
    range_ = range__;
    return *this;
  }
  PSGraphT & title(const Text & title__) {
    title_ = title__;
    return *this;
  }
  PSGraphT & x_label(const Text & x_label__) {
    labels_[0] = x_label__;
    return *this;
  }
  PSGraphT & y_label(const Text & y_label__) {
    labels_[1] = y_label__;
    return *this;
  }
  PSGraphT & y_label(const std::string & y_label__) {
    labels_[1].text() = y_label__;
    return *this;
  }
  PSGraphT & y_label(const char * y_label__) {
    labels_[1].text() = y_label__;
    return *this;
  }

  virtual void finalize(PSDoc & doc, Bounds bounds) {
    if (false && children().empty()) {
      return;
      std::cerr << "Empty graph" << std::endl;
      exit(1);
    }

    const Bounds saved_range{range_};

    // Set a scale for details of plot like marker size
    const double scale{pow(bounds.xw() * bounds.yw() /
                           doc.width() / doc.height(), 0.2)};

    // Graph labels
    do_labels(doc, bounds, scale);

    // Get graph range
    do_range(hist_ ? 0.0 : 0.01);

    // Axis label bounds adjustment
    const double tick_size__{tick_size() * scale};
    const double page_size{sqrt(doc.width() * doc.height())};
    const double x_page_scale{pow(bounds.xw() / page_size, 0.5)};
    const double y_page_scale{pow(bounds.yw() / page_size, 0.5)};
    const Axis x_axis{range().xl(), range().xh(),
          (this->log_x_ ? 3 : 7) * x_page_scale, this->log_x_};
    const Axis y_axis{range().yl(), range().yh(),
          (this->log_y_ ? 3 : 7) * y_page_scale, this->log_y_};
    // using Ticks = std::vector<std::pair<double, bool>>;
    const Ticks x_ticks{do_x_ticks_ ? x_axis.ticks() : Ticks()};
    const Ticks y_ticks{do_y_ticks_ ? y_axis.ticks() : Ticks()};
    if (y_ticks.size() && do_y_tick_labels_) {
      unsigned int max_size{min_tick_size};
      for (const std::pair<double, bool> tick : y_ticks) {
        if (!tick.second) continue;
        std::ostringstream sample;
        const double high_val{this->log_y_ ?
              pow(10, tick.first) : tick.first};
        sample << std::setprecision(10) << high_val;
        if (sample.str().size() > max_size)
          max_size = static_cast<unsigned int>(sample.str().size());
      }
      bounds.xl() += 0.65 * max_size * tick_size__;
    } else {
      if (do_y_tick_labels_) bounds.xl() += 0.7 * tick_size__;
    }
    if (do_x_tick_labels_) bounds.yl() += 1.3 * tick_size__;

    // Set x y scales
    double scales[2]{bounds.xw() / range().xw(), bounds.yw() / range().yw()};

    // Bounding box
    const double eb{border_width() * scale / 2};
    if (do_border_) do_bbox(doc, scale, bounds, eb);

    // Arbitrary ps
    do_ps(doc, scales, bounds, pre_ps_);

    // Axes and ticks and labels
    const double width1{0.5 * grid_width() * scale};
    const double width2{grid_width() * scale};
    const std::string dash1{minor_dash()};
    const std::string dash2{major_dash()};
    std::ostringstream grid2_out;
    grid2_out << std::setprecision(default_precision)
              << "np " << width2 << " lw " << dash2 << " sd\n";

    std::ostringstream ticks_out;
    ticks_out  << std::setprecision(default_precision)
               << "np bk c " << tick_size__ << " sf "
               << tick_width() << " lw [] 0 sd\n";

    doc << "/glx {" << "0 " << bounds.yw() + 2 * eb << " rl} def "
        << "/gly {" << bounds.xw() + 2 * eb << " 0" << " rl} def "
        << "/tmx {" << "0 " << tick_length() * scale << " rl} def "
        << "/tmy {" << tick_length() * scale << " 0" << " rl} def\n"
        << "np " << grid_color << " c " << width1 << " lw " << dash1 << " sd\n";

    unsigned int n_pname{0};
    for (const bool y : {false, true}) {
      for (const std::pair<double, bool> tick : (y ? y_ticks : x_ticks)) {
        std::ostringstream pname_out;
        pname_out << "p" << ++n_pname;
        const std::string pname{pname_out.str()};
        // Saved position command
        doc << "/" << pname << " {"
            << (y ? bounds.xl() - eb :
                (tick.first - range().xl()) * scales[0] + bounds.xl()) << " "
            << (y ? (tick.first - range().yl()) * scales[1] + bounds.yl() :
                bounds.yl() - eb) << " m} def\n";

        // Grid
        if (tick.second) {
          grid2_out << pname << " gl" << (y ? "y" : "x") << "\n";
        } else {
          doc << pname << " gl" << (y ? "y" : "x") << "\n";
        }

        // Ticks and labels
        if (tick.second) {
          ticks_out << pname << " tm" << (y ? "y" : "x") << " ";
          if ((!y && !do_x_tick_labels_) || (y && !do_y_tick_labels_)) {
            continue;
          }
          ticks_out << pname << " (" << std::setprecision(10)
                    << ((y && this->log_y_) || (!y && this->log_x_) ?
                        pow(10, tick.first) : tick.first)
                    << ") " << std::setprecision(default_precision);
          if (y) {
            ticks_out << "jrx " << 0.4 * tick_size__
                      << " sub " << 0.35 * tick_size__;
          } else {
            ticks_out << "jcx " << 1.1 * tick_size__;
          }
          ticks_out << " nrs\n";
        }
      }
    }
    doc << "sp\n"
        << grid2_out.str() << "sp [] 0 sd 0 0 0 c\n";

    // Last minute hist adjustments
    if (hist_ && children().size() == 2 &&
        typeid(children()[0]) == typeid(children()[1])) {
      prepare_hist_overlap();
    }

    // Arbitrary ps
    do_ps(doc, scales, bounds, mid_ps_);

    // Finalize series
    for (PSSeries * series : children()) {
      series->finalize(doc, bounds, range_, *this);
    }
    doc << "ic\n";
    doc << ticks_out.str() << "sp\n";

    // Auto Legend
    // const double lx{bounds.xl() + 7 * bounds.xw() / 10};
    // const double ly{bounds.yh() - bounds.yw() / 16};
    const double lx{bounds.xl() + 7 * bounds.xw() / 10};
    const double ly{bounds.yl() + bounds.yw() / 2};

    unsigned int c{0};
    std::ostringstream legend;
    const double legend_size__{legend_size() * scale};
    legend << "0 0 0 c " << legend_size__ << " sf "
           << border_width() * scale << " lw\n";
    for (PSSeries * series : children()) {
      const std::string text{series->title()};
      if (text.size()) {
        legend << "gs " << lx << " " << ly - 1.5 * c * legend_size__
               << " tr gs " << series->marker(scale).commands() << " sp gr "
               << 0.8 * legend_size__ << " " << -0.3 * legend_size__ << " m "
               << "(" << text << ") dup sw pop "
               << "/cx e def cx mx gt {/mx cx def} if "
               << "show gr\n";
        ++c;
      }
    }
    if (c) {
      const double lh{1.5 * c * legend_size__};
      doc << "/mx 0 def " << legend.str() << "np " << lx - 0.7 * legend_size__
          << " " << ly + 0.7 * legend_size__ << " m 0 " << -lh << " rl "
          << "mx " << legend_size__ * 1.8 << " add 0 rl 0 " << lh
          << " rl cp gs 0.9 0.9 0.9 c fp gr sp "
          << legend.str() << "\n";
    }

    // Arbitrary text
    write_text(doc, bounds, scale);

    // Arbitrary ps
    do_ps(doc, scales, bounds, ps_);

    // Restore range
    range_ = saved_range;
  }

  void prepare_hist_overlap();
  PSGraphT & hist(const bool hist__) {
    hist_ = hist__;
    return *this;
  }

  // Add raw postscript to graph
  PSGraphT & ps(const std::string & ps__) {
    ps_ += ps__ + "\n";
    return *this;
  }
  PSGraphT & mid_ps(const std::string & ps__) {
    mid_ps_ += ps__ + "\n";
    return *this;
  }
  PSGraphT & pre_ps(const std::string & ps__) {
    pre_ps_ += ps__ + "\n";
    return *this;
  }

  // Also make an eps version of this graph
  PSGraphT & eps(const std::string & file_name,
                 const std::string & title__ = "",
                 const unsigned int width__ = doc_defaults.width(),
                 const unsigned int height__ = doc_defaults.height(),
                 const double padding__ = doc_defaults.padding()) {
    std::unique_ptr<PSDoc> doc_{std::make_unique<PSDoc>(
        file_name, title__.size() ? title__ : file_name,
        width__, height__, padding__)};
    doc_->eps(true);
    std::unique_ptr<PSPage> page_{std::make_unique<PSPage>()};
    page_->manage(std::move(doc_));
    manage(std::move(page_));
    return *this;
  }

  PSGraphT & do_ticks(const bool do_ticks__) {
    do_x_ticks_ = do_ticks__;
    do_y_ticks_ = do_ticks__;
    return *this;
  }
  PSGraphT & do_x_ticks(const bool do_ticks__) {
    do_x_ticks_ = do_ticks__;
    return *this;
  }
  PSGraphT & do_y_ticks(const bool do_ticks__) {
    do_y_ticks_ = do_ticks__;
    return *this;
  }
  PSGraphT & do_x_tick_labels(const bool do_tick_labels__) {
    do_x_tick_labels_ = do_tick_labels__;
    return *this;
  }
  PSGraphT & do_y_tick_labels(const bool do_tick_labels__) {
    do_y_tick_labels_ = do_tick_labels__;
    return *this;
  }

  PSGraphT & do_border(const bool do_border__) {
    do_border_ = do_border__;
    return *this;
  }
  unsigned int min_tick_size{0};
  std::string grid_color{"0.6 0.6 0.6"};

 protected:
  void do_labels(PSDoc & doc, Bounds & bounds, const double scale,
                 const double extra_under_title = 0.0) {
    // Title
    if (title_) {
      Text title__{title_};
      doc << bounds.xc() << " " << bounds.yh() - 0.7 * title_.size() * scale
          << " m " << title_.ps(scale) << "\n";
      bounds.yh() -= 1.1 * title_.size() * scale;
    }

    // X Label
    if (labels_[0]) {
      doc << bounds.xh() << " " << bounds.yl() + 0.2 * labels_[0].size() * scale
          << " m " << labels_[0].ps(scale) << "\n";
      bounds.yl() += 1.1 * labels_[0].size() * scale;
    }

    // Maybe adjust top of graph
    bounds.yh() -= extra_under_title;

    // Y Label
    if (labels_[1]) {
      if (0)
      doc << bounds.xl() + 0.7 * labels_[1].size() * scale << " " << bounds.yh()
          << " m gs 90 ra " << labels_[1].ps(scale)
          << " gr\n";
      doc << bounds.xl() + 0.7 * labels_[1].size() * scale << " " << bounds.yh()
          << " m cxy 90 ra " << labels_[1].ps(scale)
          << "m -90 ra\n";
      bounds.xl() += 1.2 * labels_[1].size() * scale;
    }
  }

  void do_range(const double border = 0.01,
                const Bounds min_range =
                {unset(), nunset(), unset(), nunset()}) {
    Bounds new_range{range_};
    Bounds max_range{range_};
    for (PSSeries * series : children()) {
      series->get_range(range_, new_range, this->log_x_, this->log_y_);
      max_range.xl(min(max_range.xl(), new_range.xl()));
      max_range.xh(max(max_range.xh(), new_range.xh()));
      max_range.yl(min(max_range.yl(), new_range.yl()));
      max_range.yh(max(max_range.yh(), new_range.yh()));
    }
#if 0
    if (is_unset(range_.xl())) range_.xl(max_range.xl());
    if (is_unset(range_.xh())) range_.xh(max_range.xh());
    if (is_unset(range_.yl())) range_.yl(max_range.yl());
    if (is_unset(range_.yh())) range_.yh(max_range.yh());
#endif
    range_ = max_range;
    // Check for bad range problem and fake a good range if necessary
    if (is_unset(range_.xl()) || is_unset(range_.xh()) ||
        is_unset(range_.yl()) || is_unset(range_.yh())) {
      range_ = Bounds{0.001, 1, 0.001, 1};
    }
    if (de(range_.xl(), range_.xh())) {
      if (de(range_.xl(), 0)) {
        if (this->log_x_) {
          range_.xl() = 0.05;
          range_.xh() = 2;
        } else {
          range_.xl() -= 1;
          range_.xh() += 1;
        }
      } else {
        range_.xl() *= 0.9;
        range_.xh() *= 1.1;
      }
    }
    if (de(range_.yl(), range_.yh())) {
      if (de(range_.yl(), 0)) {
        if (this->log_y_) {
          range_.yl() = 0.05;
          range_.yh() = 2;
        } else {
          range_.yl() -= 1;
          range_.yh() += 1;
        }
      } else {
        range_.yl() *= 0.9;
        range_.yh() *= 1.1;
      }
    }

    if (range_.xl() > min_range.xl()) range_.xl() = min_range.xl();
    if (range_.xh() < min_range.xh()) range_.xh() = min_range.xh();
    if (range_.yl() > min_range.yl()) range_.yl() = min_range.yl();
    if (range_.yh() < min_range.yh()) range_.yh() = min_range.yh();
    if (this->log_x_) {
      range_.xl(log10(range_.xl()));
      range_.xh(log10(range_.xh()));
    }
    if (this->log_y_) {
      if (hist_) range_.yh() *= 1.5;
      range_.yl(log10(range_.yl()));
      range_.yh(log10(range_.yh()));
    }
    if (border > 0.0) {
      range_.xl() -= range_.xw() * border;
      range_.xh() += range_.xw() * border;
      range_.yl() -= range_.yw() * border;
      range_.yh() += range_.yw() * border;
    }
  }

  void do_bbox(PSDoc & doc, const double scale, const Bounds & bounds,
               const double eb) {
    // Draw bounding box
    const bool do_fill{fill() != "1 1 1"};
    if (do_fill) doc << fill() << " ";
    doc << bounds.xl() << " " << bounds.yh() << " "
        << bounds.xh() << " " << bounds.yh() << " "
        << bounds.xh() << " " << bounds.yl() << " "
        << bounds.xl() << " " << bounds.yl() << " "
        << border_width() * scale << (do_fill ? " bbf " : " bb ")
        << hist_width() * scale << " lw\n";

    doc << "/c2g {"
        << bounds.xl() + eb << " " << bounds.yh() - eb << " "
        << bounds.xh() - eb << " " << bounds.yh() - eb << " "
        << bounds.xh() - eb << " " << bounds.yl() + eb << " "
        << bounds.xl() + eb << " " << bounds.yl() + eb << " "
        << "cbb } def c2g\n";
  }

  void do_ps(PSDoc & doc, double * scales, const Bounds & bounds,
             const std::string & ps__) {
    if (&ps__ == &pre_ps_ && (pre_ps_.size() || mid_ps_.size() || ps_.size())) {
      doc << "/pfc { /y e def "
          << doc.width() << " mul " << " y " << doc.height() << " mul "
          << "} def\n";
      doc << "/xc { " << (this->log_x() ? "log " : "")
          << range().xl() << " sub "
          << scales[0] << " mul " << bounds.xl() << " add "
          << "} def\n";
      doc << "/yc { " << (this->log_y() ? "log " : "")
          << range().yl() << " sub "
          << scales[1] << " mul " << bounds.yl() << " add "
          << "} def\n";
      doc << "/gc { /y e " << (this->log_y() ? "log " : "")
          << "def " << (this->log_x() ? "log " : "")
          << range().xl() << " sub "
          << scales[0] << " mul " << bounds.xl() << " add "
          << " y " << range().yl() << " sub "
          << scales[1] << " mul " << bounds.yl() << " add "
          << "} def\n";
      doc << "/gfc { /y e def "
          << bounds.xw() << " mul " << bounds.xl() << " add "
          << " y " << bounds.yw() << " mul " << bounds.yl() << " add "
          << "} def\n";
      doc << "/xfc { "
          << bounds.xw() << " mul " << bounds.xl() << " add "
          << "} def\n";
      doc << "/yfc { "
          << bounds.yw() << " mul " << bounds.yl() << " add "
          << "} def\n";
    }
    if (ps__.size()) doc << ps__ << "\n";
  }

  Bounds range_;
  Text title_{};
  Text labels_[2];
  bool hist_{false};
  bool do_x_ticks_{doc_defaults.ticks()};
  bool do_y_ticks_{doc_defaults.ticks()};
  bool do_x_tick_labels_{true};
  bool do_y_tick_labels_{true};
  bool do_border_{true};
  std::string ps_{};
  std::string mid_ps_{};
  std::string pre_ps_{};
};

// Histogram class
template<class PSSeries, class Int>
class PSHistT : public PSGraphT<PSSeries> {
};

std::string dot() {
  return "-0.25 0 m 0.25 0 l";
}
double default_marker_size() {
  return 5.0;
}
std::string circle() {
  std::ostringstream result;
  result << std::setprecision(default_precision)
         <<  "0 0 " << default_marker_size() << " 0 360 arc ";
  return result.str();
}
// constexpr double twopi() { return 2 * acos(-1); }
constexpr double twopi() { return 2 * 3.14159265358979323846; }
std::string regular_polygon(const unsigned int sides, const double angle) {
  std::ostringstream result;
  result << std::setprecision(default_precision);
  for (unsigned int i{0}; i != sides; ++i) {
    const double radians{twopi() * (1.0 * i / sides + angle / 360)};
    result << default_marker_size() * cos(radians) << " "
           << default_marker_size() * sin(radians) << " ";
    result << (i ? "l" : "m") << " ";
  }
  result << "cp ";
  return result.str();
}
std::string triangle() { return regular_polygon(3, 90); }
std::string itriangle() { return regular_polygon(3, -90); }
std::string square() { return regular_polygon(4, 45); }
std::string diamond() { return regular_polygon(4, 0); }
std::string pentagon() { return regular_polygon(5, 90); }
std::string hexagon() { return regular_polygon(6, 0); }
std::string octagon() { return regular_polygon(8, 360.0 / 16); }
std::string plus() {
  std::ostringstream result;
  const double s{default_marker_size()};
  const double t{0.3 * s};
  result << std::setprecision(default_precision)
         << -s << " " << -t << " m "
         << -t << " " << -t << " l "
         << -t << " " << -s << " l "
         << t << " " << -s << " l "
         << t << " " << -t << " l "
         << s << " " << -t << " l "
         << s << " " << t << " l "
         << t << " " << t << " l "
         << t << " " << s << " l "
         << -t << " " << s << " l "
         << -t << " " << t << " l "
         << -s << " " << t << " l "
         << "cp ";
  return result.str();
}
std::string xshape() {
  return std::string("45 ra ") + plus();
}
std::string star() {
  std::ostringstream result;
  result << std::setprecision(default_precision);
  for (unsigned int i{0}; i != 10; ++i) {
    const double dist{default_marker_size() * ((i % 2) ? 0.382 : 1.0)};
    result << dist * cos(twopi() / 10 * (i + 0.5)) << " "
           << dist * sin(twopi() / 10 * (i + 0.5)) << " ";
    result << (i ? "l" : "m") << " ";
  }
  result << "cp ";
  return result.str();
}

class Marker {
 public:
  // Default values
  static constexpr double default_scale() { return 0.25; }
  static constexpr double default_line_width() { return 0.5;}
  static std::string default_color() { return "bk"; }
  static std::string default_fill_color() { return "bk"; }

  // Construct a Marker
  Marker(const std::string & path__ = circle(),
         const double scale__ = default_scale(),
         const std::string & color__ = default_color(),
         const double line_width__ = default_line_width(),
         const bool fill__ = true,
         const std::string & fill_color__ = "",  // default_fill_color(),
         const double rotation__ = 0.0) :
      path_{path__}, scale_{scale__}, color_{color__},
    line_width_{line_width__}, fill_{fill__},
      fill_color_{fill_color__.size() ? fill_color__ : color__},
    rotation_{rotation__} { }

  // Ordering
  bool operator<(const Marker & rhs) const {
    if (scale_ >= rhs.scale_ && scale_ <= rhs.scale_) {
      if (line_width_ >= rhs.line_width_ && line_width_ <= rhs.line_width_) {
        if (fill_ == rhs.fill_) {
          if (rotation_ >= rhs.rotation_ && rotation_ <= rhs.rotation_) {
            if (color_ == rhs.color_) {
              if (fill_color_ == rhs.fill_color_) {
                return path_ < rhs.path_;
              } else {
                return fill_color_ < rhs.fill_color_;
              }
            } else {
              return color_ < rhs.color_;
            }
          } else {
            return rotation_ < rhs.rotation_;
          }
        } else {
          return fill_ < rhs.fill_;
        }
      } else {
        return line_width_ < rhs.line_width_;
      }
    } else {
      return scale_ < rhs.scale_;
    }
  }

  // Marker postscript commands
  std::string setup_commands() const {
    std::ostringstream result;
    result << std::setprecision(default_precision);
    if (color_ != default_color()) {
      result << color_ << " c ";
    }
    if (line_width_ < default_line_width() ||
        line_width_ > default_line_width()) {
      result << line_width_ << " lw ";
    }
    return result.str();
  }
  std::string draw_commands() const {
    std::ostringstream result;
    result << std::setprecision(default_precision);
    if (rotation_ < 0 || rotation_ > 0) {
      result << rotation_ << " ra ";
    }
    if (scale_ < 1 || scale_ > 1) {
      result << scale_ << " sc ";
    }
    result << path_;
    if (fill_) {
      result << "gs ";
      if (fill_color_ != color_) {
        result << fill_color_ << " c ";
      }
      result << "fp gr ";
    }
    return result.str();
  }
  std::string commands() const {
    return setup_commands() + draw_commands();
  }
  operator std::string() const {
    return commands();
  }

  // Retrieve values
  std::string path() const { return path_; }
  double scale() const { return scale_; }
  std::string color() const { return color_; }
  double line_width() const { return line_width_; }
  std::string fill_color() const { return fill_color_; }
  bool fill() const { return fill_; }
  double rotation() const { return rotation_; }

  // Modify values
  Marker & path(const std::string & path__) {
    path_ = path__;
    return *this;
  }
  Marker & scale(const double scale__) {
    scale_ = scale__;
    return *this;
  }
  Marker & color(const std::string & color__) {
    color_ = color__;
    return *this;
  }
  Marker & line_width(const double line_width__) {
    line_width_ = line_width__;
    return *this;
  }
  Marker & fill_color(const std::string & fill_color__) {
    fill_color_ = fill_color__;
    return *this;
  }
  Marker & fill(const bool fill__) {
    fill_ = fill__;
    return *this;
  }
  Marker & rotation(const double rotation__) {
    rotation_ = rotation__;
    return *this;
  }
  Marker & colors(const std::string & color__) {
    fill_ = true;
    color_ = color__;
    fill_color_ = color__;
    return *this;
  }

 private:
  std::string path_;
  double scale_;
  std::string color_;
  double line_width_;
  bool fill_;
  std::string fill_color_;
  double rotation_;
};

class PSSeries : public Multiplexer<PSPartT<PSSeries>, PSSeries, PSSeries> {
 public:
  using PSPart = PSPartT<PSSeries>;
  using PSPage = PSPageT<PSPart>;
  using PSDoc = PSDocT<PSPage>;
  using Multi = Multiplexer<PSPart, PSSeries, PSSeries>;
  using Multi::add;
  using Multi::manage;

  PSSeries() : Multi{this} { }
  explicit PSSeries(PSPart & graph__) : Multi{this, &graph__} { }
  PSSeries(PSSeries &&) = default;
  virtual ~PSSeries() { }
  virtual void get_range(const Bounds & range, Bounds & newrange,
                         const bool log_x, const bool log_y) = 0;
  virtual void finalize(PSDoc & doc,
                        const Bounds & bounds,
                        const Bounds & range,
                        const PSPart & graph) = 0;
  virtual double y_val(const uint64_t i) const = 0;
  virtual uint64_t size() const = 0;
  virtual std::string title() const {
    return "";
  }
  virtual Marker marker(const double scale__ = 1.0) const {
    return Marker(circle(), scale__);
  }
  PSPart & graph() {
    return *parents().front();
  }
};

std::string mix_colors(const std::string & c1, const std::string & c2) {
  double rgbs[2][3];
  // Parse colors
  for (const unsigned int ci : {0, 1}) {
    const std::string & color{ci ? c2 : c1};
    std::istringstream cs{color.c_str()};
    for (const unsigned int cj : {0, 1, 2}) {
      cs >> rgbs[ci][cj];
    }
    if (!cs) throw Error("Problem parsing color") << color;
  }

  double rgb[3];
  std::ostringstream c3;
  c3 << std::setprecision(3);
  for (const unsigned int cj : {0, 1, 2}) {
    rgb[cj] = 0;
    for (const unsigned int ci : {0, 1}) {
      rgb[cj] += 0.5 * rgbs[ci][cj] * rgbs[ci][cj];
    }
    rgb[cj] = std::min(255.0, sqrt(rgb[cj]));
    if (cj) c3 << " ";
    c3 << rgb[cj];
  }

  if (0) std::cout << "combining colors " << c1 << " and " << c2
                   << " gives " << c3.str() << std::endl;
  return c3.str();
}

template <class ValType, class CountType>
class PSHSeries : public PSSeries {
 public:
  using PSGraph = PSGraphT<PSSeries>;

  // factory
  template <class ... Args>
  static PSHSeries * create(Args && ... args) {
    PSHSeries * xyseries__{new PSHSeries{std::forward<Args>(args)...}};
    paa::ownp(xyseries__);
    return xyseries__;
  }

  PSHSeries(PSPart & hist__,
            const uint64_t n_bins__ = 100,
            const std::string & color__ = "0 0 0",
            const bool normalize__ = true,
            const std::string & title__ = "") :
      PSSeries{hist__}, n_bins_{n_bins__},
    color_{color__}, normalize_{normalize__}, title_{title__} {
      static_cast<PSGraph &>(hist__).hist(true);
    }
  PSHSeries(PSDoc & doc__,
            const std::string & title__,
            const Bounds & range__ = Bounds{},
            const uint64_t n_bins__ = 100,
            const std::string & color__ = "1 0 0") :
      n_bins_{n_bins__},
    color_{color__}, normalize_{false}, title_{""} {
      std::unique_ptr<PSGraph> graph_{std::make_unique<PSGraph>(
          doc__, title__, range__)};
      graph_->hist(true);
      manage(std::move(graph_));
    }
  PSHSeries(PSDoc & doc__,
            const std::string & title__,
            const uint64_t n_bins__ = 100,
            const Bounds & range__ = Bounds{},
            const std::string & color__ = "1 0 0") :
      n_bins_{n_bins__},
    color_{color__}, normalize_{false}, title_{""} {
      std::unique_ptr<PSGraph> graph_{std::make_unique<PSGraph>(
          doc__, title__, range__)};
      graph_->hist(true);
      manage(std::move(graph_));
    }
  PSHSeries(PSPage & page__,
            const std::string & title__,
            const Bounds & range__,
            const uint64_t n_bins__ = 100,
            const std::string & color__ = "1 0 0",
            const bool normalize__ = false) :
      n_bins_{n_bins__},
    color_{color__}, normalize_{normalize__}, title_{""} {
      std::unique_ptr<PSGraph> graph_{std::make_unique<PSGraph>(
          page__, title__, range__)};
      graph_->hist(true);
      manage(std::move(graph_));
    }
  PSHSeries(PSPage & page__,
            const std::string & title__,
            const uint64_t n_bins__ = 100,
            const Bounds & range__ = Bounds{0.0},
            const std::string & color__ = "1 0 0",
            const bool normalize__ = false) :
      n_bins_{n_bins__},
    color_{color__}, normalize_{normalize__}, title_{""} {
      std::unique_ptr<PSGraph> graph_{std::make_unique<PSGraph>(
          page__, title__, range__)};
      graph_->hist(true);
      manage(std::move(graph_));
    }
  template<class Hist>
  PSHSeries(PSPage & page__,
            const Hist & histogram,
            const std::string & title__,
            const std::string & color__ = "1 0 0",
            const bool normalize__ = false) :
      h_(histogram.n_bins()), n_bins_{histogram.n_bins()},
      color_{color__}, normalize_{normalize__}, title_{""}, filled{true} {
      std::unique_ptr<PSGraph> graph_{std::make_unique<PSGraph>(
          page__, title__,
          Bounds{1.0 * histogram.min(), 1.0 * histogram.max()})};
      graph_->hist(true);
      manage(std::move(graph_));
      for (uint64_t b{0}; b != n_bins_; ++b) h_[b] = histogram[b];
    }

  explicit PSHSeries(std::vector<PSSeries *> others) {
    const PSSeries * first{others.front()};
    const PSSeries * second{others.back()};
    n_bins_ = first->size();
    // color_ = "0.8 0 0.8";
    color_ = mix_colors(first->marker().fill_color(),
                        second->marker().fill_color());
    // color_ = mix_colors("1 0 0", "0 0 1");
    normalize_ = false;
    h_.resize(n_bins_);
#if 0
    if (first->title().size() || second->title().size()) {
      title_ = "Area in Common";
    }
#endif
    for (uint64_t i{0}; i != h_.size(); ++i) {
      h_[i] = std::min(first->y_val(i), second->y_val(i));
    }
  }
  virtual ~PSHSeries() {
    if (VERBOSE) std::cerr << "Destroy PSHSeries " << std::endl;
    this->finalize_doc();
  }
  PSHSeries(PSHSeries &&) = default;
  void add_point(const ValType x__) {
    x_.push_back(x__);
  }

  void add_value(const ValType x__, const uint64_t count) {
    for (uint64_t i{0}; i != count; ++i) {
      x_.push_back(x__);
    }
  }

  void get_range(const Bounds & range, Bounds & new_range,
                 const bool, const bool log_y) {
    if (x_.empty() && !filled) return;

    const bool get_range_low{is_unset(range.xl())};
    const bool get_range_high{is_unset(range.xh())};
    const bool do_get_range{get_range_low || get_range_high};

    if (do_get_range) {
      const auto minmax = std::minmax_element(x_.begin(), x_.end());
      if (get_range_low) new_range.xl(*minmax.first);
      if (get_range_high) {
        new_range.xh(*minmax.second);
        if (n_bins_)
          new_range.xh(*minmax.second + new_range.xw() / (100 * n_bins_));
      }
    }

    if (!filled) {
      h_.clear();
      if (n_bins_ == 0) {
        new_range.xl() -= 1;
        new_range.xh() += 2;
        n_bins_ = new_range.xh() - new_range.xl();
      }
      h_.resize(n_bins_);

      for (const ValType val : x_) {
        const uint64_t bin{
          static_cast<uint64_t>((
              val - new_range.xl()) * n_bins_ / new_range.xw())};
        if (bin >= 0 && bin < n_bins_) {
          ++h_[bin];
        }
      }
    }
    for (const CountType val : h_) {
      const double nval{1.05 *
            (normalize_ ? 1.0 * val / x_.size() : 1.0 * val)};
      if (new_range.yh() < nval) {
        new_range.yh() = nval;
      }
      if (new_range.yl() > nval) {
        new_range.yl() = nval;  // useless see last line
      }
    }
    if (normalize_ && new_range.yh() > 1) {
      new_range.yh(1);
    }
    new_range.yl(log_y ? 0.5 : 0);
    if (log_y && new_range.yh() <= 1) new_range.yh(10);
  }

  double y_val(const uint64_t i) const {
    return normalize_ ? 1.0 * h_[i] / x_.size() : 1.0 * h_[i];
  }
  uint64_t size() const { return h_.size(); }
  virtual void finalize(PSDoc & doc,
                const Bounds & bounds,
                const Bounds & range,
                const PSPart & part) {
    if (x_.empty() && !filled) return;

    const double scales[2]{bounds.xw() / range.xw(), bounds.yw() / range.yw()};
    const double binw{bounds.xw() / n_bins_};
    if (h_.size()) {
      doc << "bk c np /h {" << bounds.yl() << " m "
          << binw << " 0 rl 0 e rl "
          << -binw << " 0 rl cp ";
      if (color_.size()) doc << "gs " << color_ << " c fp gr ";
      doc << "sp} def\n";
      for (uint64_t i{0}; i != h_.size(); ++i)
        if (h_[i] > 0)
          doc << ((part.log_y() ? log10(y_val(i)) : y_val(i)) -
                  range.yl()) * scales[1]
              << " " << bounds.xl() + i * binw << " h"
              << "\n";
      doc << "sp\n";
    }
  }

  std::string color() const { return color_; }
  std::string title() const { return title_; }
  Marker marker(const double scale__ = 1.0) const {
    return Marker{square(), 1.5 * scale__, "bk", 1.0, true, color_};
  }

  PSGraph & graph() {
    return dynamic_cast<PSGraph&>(*parents().front());
  }

 private:
  std::vector<ValType> x_{};
  std::vector<CountType> h_{};
  uint64_t n_bins_{};
  std::string color_{};
  bool normalize_{};
  std::string title_{};
  bool filled{false};
};

class PSXYSeries : public PSSeries {
 public:
  using PSGraph = PSGraphT<PSSeries>;
  // using PSPart = PSPartT<PSSeries>;
  using PSPage = PSPageT<PSPart>;
  using PSDoc = PSDocT<PSPage>;
  using PSSeries::add;
  using PSSeries::manage;

  // factory
  template <class ... Args>
  static PSXYSeries * create(Args && ... args) {
    PSXYSeries * xyseries__{new PSXYSeries{std::forward<Args>(args)...}};
    paa::ownp(xyseries__);
    return xyseries__;
  }

  // Construct
  explicit PSXYSeries(const Marker & marker__ = Marker()) :
      marker_{marker__},
      draw_commands_{marker_.draw_commands()},
      setup_commands_{marker_.setup_commands()} { }
  PSXYSeries(PSPart & graph__,
             const std::string & draw_commands__,
             const std::string & setup_commands__ = "") :
      PSSeries{graph__},
    draw_commands_{draw_commands__},
    setup_commands_{setup_commands__} { }
  explicit PSXYSeries(PSPart & graph__) :
      PSXYSeries{graph__, Marker().draw_commands(),
        Marker().setup_commands()} { }
  PSXYSeries(PSPart & graph__, const Marker & marker__,
             const std::string & title__ = "") :
      PSSeries{graph__}, marker_{marker__},
      draw_commands_{marker__.draw_commands()},
    setup_commands_{marker__.setup_commands()},
    title_{title__} { }
  explicit PSXYSeries(PSDoc & doc) :
      marker_{Marker()},
      draw_commands_{marker_.draw_commands()},
      setup_commands_{marker_.setup_commands()} {
      std::unique_ptr<PSGraph> graph_{std::make_unique<PSGraph>(doc)};
      manage(std::move(graph_));
    }
  template <class Hist>
  explicit PSXYSeries(const Hist & hist, PSDoc & doc,
                      const std::string & title__ = "") :
      marker_{Marker()},
      draw_commands_{marker_.draw_commands()},
      setup_commands_{marker_.setup_commands()} {
      std::unique_ptr<PSGraph> graph_{std::make_unique<PSGraph>(doc, title__)};
      manage(std::move(graph_));
      for (uint64_t bin{0}; bin != hist.n_bins(); ++bin)
        add_point(hist.mid(bin), hist[bin]);
    }
  PSXYSeries(PSDoc & doc,
             const std::string & title__) :
      PSXYSeries(doc, title__, Bounds{0.0}, Marker{}) {}
  PSXYSeries(PSDoc & doc,
             const std::string & title__,
             const Bounds & range__,
             const Marker & marker__ = Marker()) :
      marker_{marker__},
      draw_commands_{marker__.draw_commands()},
      setup_commands_{marker__.setup_commands()} {
      std::unique_ptr<PSGraph> graph_{
        std::make_unique<PSGraph>(doc, title__, range__)};
      manage(std::move(graph_));
    }
  PSXYSeries(PSDoc & doc,
           const std::string & title__,
           const Marker & marker__,
           const Bounds & range__ = Bounds{0.0}) :
      PSXYSeries(doc, title__, range__, marker__) {}
  explicit PSXYSeries(const std::string & file_name) :
    draw_commands_{Marker().draw_commands()},
    setup_commands_{Marker().setup_commands()} {
      std::unique_ptr<PSGraph> graph_{std::make_unique<PSGraph>(file_name)};
      manage(std::move(graph_));
    }
  PSXYSeries(PSPage & page,
             const std::string & title__) :
      marker_{Marker()},
      draw_commands_{Marker().draw_commands()},
      setup_commands_{Marker().setup_commands()} {
      std::unique_ptr<PSGraph> graph_{
        std::make_unique<PSGraph>(page, title__, Bounds{0.0})};
      manage(std::move(graph_));
    }
  PSXYSeries(PSPage & page,
             const std::string & title__,
             const Bounds & range__,
             const Marker & marker__ = Marker()) :
      marker_{marker__},
      draw_commands_{marker__.draw_commands()},
      setup_commands_{marker__.setup_commands()} {
      std::unique_ptr<PSGraph> graph_{
        std::make_unique<PSGraph>(page, title__, range__)};
      manage(std::move(graph_));
    }
  PSXYSeries(PSPage & page,
             const std::string & title__,
             const Marker & marker__,
             const Bounds & range__ = Bounds{0.0}) :
      marker_{marker__},
      draw_commands_{marker__.draw_commands()},
      setup_commands_{marker__.setup_commands()} {
      std::unique_ptr<PSGraph> graph_{
        std::make_unique<PSGraph>(page, title__, range__)};
      manage(std::move(graph_));
    }
  PSXYSeries(const std::string & file_name,
             const std::string & title__,
             const Bounds & range__ = Bounds{0.0},
             const Marker & marker__ = Marker()) :
      marker_{marker__},
      draw_commands_{marker_.draw_commands()},
      setup_commands_{marker_.setup_commands()} {
        std::unique_ptr<PSGraph> graph_{
          std::make_unique<PSGraph>(file_name, title__, range__)};
      manage(std::move(graph_));
    }
  PSXYSeries(const std::string & file_name,
             const std::string & title__,
             const Marker & marker__,
             const Bounds & range__ = Bounds{0.0}) :
      PSXYSeries{file_name, title__, range__, marker__} { }
  PSXYSeries(PSXYSeries && other) = default;
  virtual ~PSXYSeries() {
    if (VERBOSE) std::cerr << "Destroy PSXYSeries " << std::endl;
    this->finalize_doc();
  }

  PSGraph & graph() {
    return dynamic_cast<PSGraph&>(*parents().front());
  }

  // Add data point
  void add_point(const double x_, const double y_) {
    x.push_back(x_);
    y.push_back(y_);
  }
  template <class Position>
  void add_point(const Position pos) {
    x.push_back(pos.x());
    y.push_back(pos.y());
  }

  // Set X Y limits
  virtual void get_range(const Bounds & range, Bounds & new_range,
                         const bool log_x, const bool log_y) {
    if (x.empty()) return;

    // Get x y ranges
    if (is_unset(range.xl()) && is_unset(range.xh())) {
      set_x_limits(new_range.xl(), new_range.xh(), log_x);
    } else if (is_unset(range.xl())) {
      set_xl(new_range.xl(), log_x);
    } else if (is_unset(range.xh())) {
      set_xh(new_range.xh());
    }
    if (is_unset(range.yl()) && is_unset(range.yh())) {
      set_y_limits(new_range.yl(), new_range.yh(), log_y);
    } else if (is_unset(range.yl())) {
      set_yl(new_range.yl(), log_y);
    } else if (is_unset(range.yh())) {
      set_yh(new_range.yh());
    }
  }
  void set_x_limits(double & xl, double & xh, const bool log_x) const {
    if (x.empty()) return;
    if (log_x) {
      for (const double val : x) {
        if (xl > val && val > 0) {
          xl = val;
        }
        if (xh < val) {
          xh = val;
        }
      }
    } else {
      const auto minmax = minmax_element(x.begin(), x.end());
      xl = std::min(*minmax.first, xl);
      xh = std::max(*minmax.second, xh);
    }
  }
  void set_y_limits(double & yl, double & yh, const bool log_y) const {
    if (y.empty()) return;
    if (log_y) {
      for (const double val : y) {
        if (yl > val && val > 0) {
          yl = val;
        }
        if (yh < val) {
          yh = val;
        }
      }
    } else {
      const auto minmax = minmax_element(y.begin(), y.end());
      yl = std::min(*minmax.first, yl);
      yh = std::max(*minmax.second, yh);
    }
  }
  void set_xl(double & xl, const bool log_x) const {
    if (x.empty()) return;
    if (log_x) {
      for (const double val : x) {
        if (xl > val && val > 0) {
          xl = val;
        }
      }
    } else {
      xl = *min_element(x.begin(), x.end());
    }
  }
  void set_xh(double & xh) const {
    if (x.empty()) return;
    xh = *max_element(x.begin(), x.end());
  }
  void set_yl(double & yl, const bool log_y) const {
    if (y.empty()) return;
    if (log_y) {
      for (const double val : y) {
        if (yl > val && val > 0) {
          yl = val;
        }
      }
    } else {
      yl = *min_element(y.begin(), y.end());
    }
  }
  void set_yh(double & yh) const {
    if (y.empty()) return;
    yh = *max_element(y.begin(), y.end());
  }

  double y_val(const uint64_t i) const { return y[i]; }
  uint64_t size() const { return x.size(); }

  // Finalize
  virtual void finalize(PSDoc & doc, const Bounds & bounds,
                const Bounds & range, const PSPart & part) {
    if (x.empty()) return;

    const double scales[2]{bounds.xw() / range.xw(),
          bounds.yw() / range.yw()};
    const double scale{pow(bounds.xw() * bounds.yw() /
                           doc.width() / doc.height(), 0.2)};
    if (setup_commands_.size()) doc << setup_commands_  << "\n";
    doc << "/p {gs tr np " << scale << " sc "
        << draw_commands_ << "sp gr} def\n";
    std::ostringstream lines_out;
    lines_out << std::setprecision(default_precision);
    double x_start{0.0};
    double last_x{0.0};
    double last_y{0.0};
    if (do_lines()) lines_out << "0.5 lw np\n";
    if (e.size() && error_bars) {
      if (part.log_y()) throw Error("No error bars allowed on log graphs");
      const double bar{16 * scale};
      doc << "/ebar {/EE e def /EY e def /EX e def gs 2 " << scale
          << " mul lw "
          << "gs np EX " << bar / 2 << " sub EY EE sub m " << bar
          << " 0 rl EX " << bar / 2 << " sub EY EE add m " << bar
          << " 0 rl EX EY EE sub m 0 EE 2 mul rl sp gr} def\n";
    }
    for (unsigned int i{0}; i != x.size(); ++i) {
      if (part.log_x() && x[i] <= 0) continue;
      if (part.log_y() && y[i] <= 0) continue;
      const double xv{part.log_x() ? log10(x[i]) : x[i]};
      const double yv{part.log_y() ? log10(y[i]) : y[i]};
      if (range.includes(xv, yv) || e.size()) {
        if (e.size() && e[i] > 0 && error_bars)
          doc << bounds.xl() + (xv - range.xl()) * scales[0] << " "
              << bounds.yl() + (yv - range.yl()) * scales[1] << " "
              << e[i] * scales[1] << " ebar\n";
        doc << bounds.xl() + (xv - range.xl()) * scales[0] << " "
            << bounds.yl() + (yv - range.yl()) * scales[1] << " " << "p\n";
        if (do_lines()) {
          if (i) {
            if (yv < last_y || yv > last_y || i + 1 == x.size()) {
              lines_out << bounds.xl() + (last_x - range.xl()) * scales[0]
                        << " "
                        << bounds.yl() + (last_y - range.yl()) * scales[1]
                        << " m "
                        << bounds.xl() + (xv - range.xl()) * scales[0]
                        << " "
                        << bounds.yl() + (yv - range.yl()) * scales[1]
                        << " l\n";
              if (x_start < last_x || x_start > last_x || i + 1 == x.size()) {
                lines_out << bounds.xl() + (x_start - range.xl()) * scales[0]
                          << " "
                          << bounds.yl() + (last_y - range.yl()) * scales[1]
                          << " m "
                          << bounds.xl() + (last_x - range.xl()) * scales[0]
                          << " "
                          << bounds.yl() + (last_y - range.yl()) * scales[1]
                          << " l\n";
              }
              x_start = xv;
              last_y = yv;
            }
          } else {
            x_start = xv;
            last_y = yv;
          }
          last_x = xv;
        }
      }
    }
    if (do_lines()) {
      lines_out << "sp\n";
      doc << lines_out.str();
    }
  }
  std::string title() const { return title_; }
  Marker marker(const double scale__ = 1.0) const {
    if (0) std::cerr << scale__ << std::endl;
    return marker_;
  }
  bool do_lines() const { return do_lines_; }
  PSXYSeries & do_lines(const bool do_lines__) {
    do_lines_ = do_lines__;
    return *this;
  }

  bool error_bars{true};

 protected:
  std::vector<double> x{};
  std::vector<double> y{};
  std::vector<double> e{};  // only used by XYESeries class below
  Marker marker_{};
  std::string draw_commands_{};
  std::string setup_commands_{};
  std::string title_{};
  bool do_lines_{false};
};

class PSXYESeries : public PSXYSeries {
 public:
  using PSSeries = PSXYSeries::PSSeries;
  using PSGraph = PSGraphT<PSSeries>;
  using PSPage = PSPageT<PSPart>;
  using PSDoc = PSDocT<PSPage>;
  using PSSeries::add;
  using PSSeries::manage;

  // factory
  template <class ... Args>
  static PSXYESeries * create(Args && ... args) {
    PSXYESeries * xyseries__{new PSXYESeries{std::forward<Args>(args)...}};
    paa::ownp(xyseries__);
    return xyseries__;
  }
  using PSXYSeries::PSXYSeries;

  void add_point_error(const double x_, const double y_, const double e_) {
    x.push_back(x_);
    y.push_back(y_);
    e.push_back(e_);
  }

  // Set X Y limits
  virtual void get_range(const Bounds & range, Bounds & new_range,
                         const bool log_x, const bool log_y) {
    // Calculate averages, errors
    if (e.empty()) {
      std::unordered_map<double, RunningMean> means;
      for (uint64_t i{0}; i != size(); ++i) means[x[i]] += y[i];
      x.clear();
      y.clear();
      std::set<double> values;
      for (const auto & vals : means) values.insert(vals.first);
      for (const double value : values) {
        const RunningMean & mean{means.at(value)};
        x.push_back(value);
        y.push_back(mean.mean());
        e.push_back(seom ? mean.seom() : mean.stdev());
      }
    } else {
      errors_range = true;
    }
    PSXYSeries::get_range(range, new_range, log_x, log_y);
    if (errors_range) {
      std::pair<double, double> erange{unset(), nunset()};
      for (uint64_t i{0}; i != y.size(); ++i) {
        erange.first = min(erange.first, y[i] - e[i]);
        erange.second = max(erange.second, y[i] + e[i]);
      }
      new_range.yl(erange.first);
      new_range.yh(erange.second);
    }
  }
  bool seom{false};
  bool errors_range{false};
};


class PSXYDSeries : public PSSeries {
 public:
  using PSGraph = PSGraphT<PSSeries>;
  // using PSPart = PSPartT<PSSeries>;
  using PSPage = PSPageT<PSPart>;
  using PSDoc = PSDocT<PSPage>;
  using PSSeries::add;
  using PSSeries::manage;

  // Construct
  PSXYDSeries(PSPart & graph__,
              const unsigned int x_bins__ = 100,
              const unsigned int y_bins__ = 100) :
      PSSeries{graph__},
    x_bins_{x_bins__}, y_bins_{y_bins__} { }
  PSXYDSeries(PSXYDSeries && other) = default;
  virtual ~PSXYDSeries() {
    if (VERBOSE) std::cerr << "Destroy PSXYDSeries " << std::endl;
    this->finalize_doc();
  }

  // Add data point
  void add_point(const double x_, const double y_) {
    x.push_back(x_);
    y.push_back(y_);
  }

  // Set X Y limits
  void get_range(const Bounds & range, Bounds & new_range,
                 const bool, const bool) {
    // Get x y ranges
    if (is_unset(range.xl()) && is_unset(range.xh())) {
      set_x_limits(new_range.xl(), new_range.xh());
    } else if (is_unset(range.xl())) {
      set_xl(new_range.xl());
    } else if (is_unset(range.xh())) {
      set_xh(new_range.xh());
    }
    if (is_unset(range.yl()) && is_unset(range.yh())) {
      set_y_limits(new_range.yl(), new_range.yh());
    } else if (is_unset(range.yl())) {
      set_yl(new_range.yl());
    } else if (is_unset(range.yh())) {
      set_yh(new_range.yh());
    }
  }
  void set_x_limits(double & xl, double & xh) const {
    if (x.empty()) return;
    const auto minmax = minmax_element(x.begin(), x.end());
    xl = std::min(*minmax.first, xl);
    xh = std::max(*minmax.second, xh);
  }
  void set_y_limits(double & yl, double & yh) const {
    if (y.empty()) return;
    const auto minmax = minmax_element(y.begin(), y.end());
    yl = std::min(*minmax.first, yl);
    yh = std::max(*minmax.second, yh);
  }
  void set_xl(double & xl) const {
    if (x.empty()) return;
    xl = *min_element(x.begin(), x.end());
  }
  void set_xh(double & xh) const {
    if (x.empty()) return;
    xh = *max_element(x.begin(), x.end());
  }
  void set_yl(double & yl) const {
    if (y.empty()) return;
    yl = *min_element(y.begin(), y.end());
  }
  void set_yh(double & yh) const {
    if (y.empty()) return;
    yh = *max_element(y.begin(), y.end());
  }

  double y_val(const uint64_t i) const { return y[i]; }
  uint64_t size() const { return x.size(); }

  // Finalize
  virtual void finalize(PSDoc & doc, const Bounds & bounds,
                const Bounds & range, const PSPart &) {
    const double scales[2]{bounds.xw() / range.xw(),
          bounds.yw() / range.yw()};
    std::vector<std::vector<unsigned int>> hist(
        x_bins_, std::vector<unsigned int>(y_bins_));
    const double bw{range.xw() / x_bins_};
    const double bh{range.yw() / y_bins_};
    for (unsigned int i{0}; i != x.size(); ++i) {
      const double xv{x[i]};
      const double yv{y[i]};
      const unsigned int xb{static_cast<unsigned int>(
          (xv - range.xl()) / bw)};
      const unsigned int yb{static_cast<unsigned int>(
          (yv - range.yl()) / bh)};
      if (xb < x_bins_ && yb < y_bins_)
        ++hist[xb][yb];
    }
    unsigned int max_count{0};
    for (unsigned int xi{0}; xi != x_bins_; ++xi) {
      for (unsigned int yi{0}; yi != y_bins_; ++yi) {
        max_count = max(max_count, hist[xi][yi]);
      }
    }
    doc << "/hb {gs /y e def /x e def c np "
        << "/xh " << bw * scales[0] << " x add def "
        << "/yh " << bh * scales[1] << " y add def "
        << "x y m xh y l xh yh l x yh l cp fp "
        << "sp gr} def\n";

    // Histogram equalization of counts
    std::vector<unsigned int> counts;
    counts.reserve(x_bins_ * y_bins_);
    for (unsigned int xi{0}; xi != x_bins_; ++xi) {
      for (unsigned int yi{0}; yi != y_bins_; ++yi) {
        counts.push_back(hist[xi][yi]);
      }
    }
    sort(counts.begin(), counts.end());
    std::vector<unsigned int> levels;
    const unsigned int n_levels{100};
    for (unsigned int l{0}; l != n_levels; ++l) {
      levels.push_back(counts[counts.size() * 1.0 * l / n_levels]);
    }
    for (unsigned int xi{0}; xi != x_bins_; ++xi) {
      for (unsigned int yi{0}; yi != y_bins_; ++yi) {
        const unsigned int val{hist[xi][yi]};
        const unsigned int level{static_cast<unsigned int>(
            lower_bound(levels.begin(), levels.end(), val) - levels.begin())};
        const std::string cf{std::to_string(
            pow(1 - 1.0 * level / n_levels, 0.5))};
        const std::string color{"1 " + cf + " " + cf};
        doc << color << " "
              << xi * bounds.xw() / x_bins_ + bounds.xl() << " "
              << yi * bounds.yw() / y_bins_ + bounds.yl() << " hb\n";
      }
    }

    if (0) {
      for (unsigned int xi{0}; xi != x_bins_; ++xi) {
        for (unsigned int yi{0}; yi != y_bins_; ++yi) {
          const double frac{1.0 * hist[xi][yi] / max_count};
          const std::string cf{std::to_string(pow(1 - frac, 5))};
          const std::string color{"1 " + cf + " " + cf};
          doc << color << " "
              << xi * bounds.xw() / x_bins_ + bounds.xl() << " "
              << yi * bounds.yw() / y_bins_ + bounds.yl() << " hb\n";
        }
      }
    }
  }

 protected:
  std::vector<double> x{};
  std::vector<double> y{};
  unsigned int x_bins_{};
  unsigned int y_bins_{};
};


class PSXYMSeries : public PSXYSeries {
 public:
  using PSDoc = PSDocT<PSPage>;

  // factory
  template <class ... Args>
  static PSXYMSeries * create(Args && ... args) {
    PSXYMSeries * xyseries__{new PSXYMSeries{std::forward<Args>(args)...}};
    paa::ownp(xyseries__);
    return xyseries__;
  }

  // using PSXYSeries::PSXYSeries;
  explicit PSXYMSeries(PSPart & graph__) :
      PSXYSeries{graph__} { }
  PSXYMSeries(PSPart & graph__,
              const std::string & draw_commands__,
              const std::string & setup_commands__ = "") :
      PSXYSeries{graph__, draw_commands__, setup_commands__} { }
  PSXYMSeries(PSDoc & doc,
              const std::string & title__,
              const Bounds & range__ = Bounds{0.0},
              const Marker & marker__ = Marker()) :
      PSXYSeries{doc, title__, range__, marker__} {}
  PSXYMSeries(PSPage & page,
              const std::string & title__,
              const Marker & marker__ = Marker(),
              const Bounds & range__ = Bounds{0.0}) :
      PSXYSeries{page, title__, marker__, range__} {}
  PSXYMSeries(PSXYMSeries && other) = default;
  virtual ~PSXYMSeries() {
    if (VERBOSE) std::cerr << "Destroy PSXYMSeries " << std::endl;
    this->finalize_doc();
  }

  // Add data point
  void add_point(const double x_, const double y_) {
    x.push_back(x_);
    y.push_back(y_);
    m.push_back(nullptr);
  }
  void add_point(const double x_, const double y_, const Marker & m_) {
    x.push_back(x_);
    y.push_back(y_);
    m.push_back(&*markers.insert(m_).first);
  }
  template <class Position>
  void add_point(const Position pos, const Marker & m_) {
    x.push_back(pos.x());
    y.push_back(pos.y());
    m.push_back(&*markers.insert(m_).first);
  }

  // Finalize
  virtual void finalize(PSDoc & doc, const Bounds & bounds,
                const Bounds & range, const PSPart & part) {
    const double scales[2]{bounds.xw() / range.xw(),
          bounds.yw() / range.yw()};
    const double scale{pow(bounds.xw() * bounds.yw() /
                           doc.width() / doc.height(), 0.2)};
    if (setup_commands_.size()) doc << setup_commands_  << "\n";
    doc << "/p {gs tr np " << scale << " sc "
        << draw_commands_ << "sp gr} def\n";
    for (unsigned int i{0}; i != x.size(); ++i) {
      if (part.log_x() && x[i] <= 0) continue;
      if (part.log_y() && y[i] <= 0) continue;
      const double xv{part.log_x() ? log10(x[i]) : x[i]};
      const double yv{part.log_y() ? log10(y[i]) : y[i]};
      if (range.includes(xv, yv)) {
        doc << bounds.xl() + (xv - range.xl()) * scales[0] << " "
            << bounds.yl() + (yv - range.yl()) * scales[1] << " ";
        if (m[i]) {
          doc << "gs " << m[i]->setup_commands() << "tr np "
              << scale << " sc " << m[i]->draw_commands()
              << "sp gr";
        } else {
          doc << "p";
        }
        doc << "\n";
      }
    }
  }

 private:
  std::vector<const Marker *> m{};
  std::set<Marker> markers{};
};

class HistEq {
 public:
  HistEq(std::vector<double> values, const uint64_t n_choices) {
    sort(values.begin(), values.end());
    for (uint64_t bin{0}; bin != n_choices; ++bin) {
      const double value{values[static_cast<uint64_t>(
          1.0 * bin * values.size() / n_choices)]};
      stops.push_back(value);
    }
  }
  uint64_t operator()(const double value) const {
    return lower_bound(stops.begin(), stops.end(), value) - stops.begin();
  }
 private:
  std::vector<double> stops{};
};

class PSXYSSeries : public PSXYSeries {
 public:
  using PSDoc = PSDocT<PSPage>;

  // factory
  template <class ... Args>
  static PSXYSSeries * create(Args && ... args) {
    PSXYSSeries * xyseries__{new PSXYSSeries{std::forward<Args>(args)...}};
    paa::ownp(xyseries__);
    return xyseries__;
  }

  PSXYSSeries() : PSXYSeries{} { }
  explicit PSXYSSeries(PSPart & graph__, const double side__ = 1.0) :
      PSXYSeries{graph__},
      side_{side__} { }
  PSXYSSeries(PSXYSSeries && other) = default;

  virtual ~PSXYSSeries() {
    if (VERBOSE) std::cerr << "Destroy PSXYSSeries " << std::endl;
    this->finalize_doc();
  }

  // Add data point
  void set_value(const double x_, const double y_, const double z) {
    add_point(x_, y_);
    values.push_back(z);
  }
  void set_value_symmetric(const double x_, const double y_, const double z) {
    set_value(x_, y_, z);
    if (x_ < y_ || x_ > y_) set_value(y_, x_, z);
  }

  virtual void get_range(const Bounds & range, Bounds & new_range,
                         const bool log_x, const bool log_y) {
    PSXYSeries::get_range(range, new_range, log_x, log_y);
    new_range.xh() += side_;
    new_range.yh() += side_;
  }

  // Finalize
  virtual void finalize(PSDoc & doc, const Bounds & bounds,
                        const Bounds & range, const PSPart &) {
    const double scales[2]{bounds.xw() / range.xw(),
          bounds.yw() / range.yw()};

    if (histeq_) hist_eq = std::make_unique<HistEq>(values, n_values);

    std::function<double(double)> identity_color{
      [] (const double value) {
        return value;
      }};
    const auto minmax = minmax_element(values.begin(), values.end());
    std::function<double(double)> expanded_color{
      [minmax](const double value) {
        return (value - *minmax.first) / (*minmax.second - *minmax.first);
      }};
    std::function<double(double)> histeq_color{
        [this](const double value) {
          return 1.0 * ((*hist_eq)(value)) / n_values;
        }};
    color = [this, identity_color, expanded_color, histeq_color]
        (const double value) {
      return bound((histeq_ ? histeq_color :
                    (expanded_ ? expanded_color : identity_color))(value),
                   0, 1);
    };
    doc << "gs "
        << bounds.xl() << " " << bounds.yl() << " tr "
        << scales[0] << " " << scales[1] << " scale "
        << -range.xl() << " " << -range.yl() << " tr\n";
    doc << "/sz " << side_ << " def "
        << "/p { np 0 0 c m "
        << "sz 0 rl 0 sz rl sz neg 0 rl 0 sz neg rl cp fp } def\n";
    for (unsigned int i{0}; i != x.size(); ++i) {
      if (false && values.size() < 200) {
        std::cerr << i << " " << x[i] << " " << y[i]
                  << " " << values[i] << " " << color(values[i]) << std::endl;
      }
      doc << x[i] << " " << y[i] << " " << color(values[i]) << " p\n";
    }
    doc << "gr\n";
  }
  bool histeq() const { return histeq_; }
  PSXYSSeries & histeq(const bool histeq__) {
    histeq_ = histeq__;
    return *this;
  }
  bool expanded() const { return expanded_; }
  PSXYSSeries & expanded(const bool expanded__) {
    expanded_ = expanded__;
    return *this;
  }
  PSXYSSeries & side(const double side__) {
    side_ = side__;
    return *this;
  }

 protected:
  double side_{1.0};
  std::function<double(double)> color{nullptr};
  const unsigned int n_values{100};
  std::vector<double> values{};
  std::unique_ptr<HistEq> hist_eq{nullptr};
  bool histeq_{false};
  bool expanded_{false};
};

class PSHeat : public PSXYSSeries {
 public:
  using PSDoc = PSDocT<PSPage>;
  using PSGraph = PSGraphT<PSSeries>;
  using PSPage = PSPageT<PSPart>;

  // factory
  template <class ... Args>
  static PSHeat * create(Args && ... args) {
    PSHeat * xyseries__{new PSHeat{std::forward<Args>(args)...}};
    paa::ownp(xyseries__);
    return xyseries__;
  }

  explicit PSHeat(PSDoc & doc__, const std::string & title__) :
      PSXYSSeries{},
      page{doc__, "", "2 (0.75) 1 ^1 0 1 4 (0.45 0.05 0.45)^"},
      graph{page, title__},
      hist{page, ";Value;N", {0, 1}, n_values},
      legend{page, ";Value;Color"},
      legend_series{legend, 1.0 / n_values},
      color_hist{page, ";Color;N", {0, 1}, n_values},
      color_legend{page, ";Color;Color"},
      color_legend_series{color_legend, 1.0 / n_values} {
        page.title(graph.title());
        graph.title("");
        add(&graph);
        // legend.do_y_ticks(false);  // broken...
        // color_legend.do_y_ticks(false);
      }
  PSHeat(PSHeat && other) = default;

  virtual ~PSHeat() {
    if (VERBOSE) std::cerr << "Destroy PSHeat " << std::endl;
    this->finalize_doc();
  }

  virtual void finalize(PSDoc & doc, const Bounds & bounds,  // NOLINT
                        const Bounds & range, const PSPart & part) {  // NOLINT
    PSXYSSeries::finalize(doc, bounds, range, part);
    if (false) std::cerr << page.title().text() << std::endl;
    if (0) std::cerr << "Display States: "
                     << " " << expanded()
                     << " " << histeq()
                     << " " << legend_series.expanded()
                     << " " << legend_series.histeq()
                     << " " << color_legend_series.expanded()
                     << " " << color_legend_series.histeq() << std::endl;
    for (uint64_t vb{0}; vb <= n_values; ++vb) {
      const double value{1.0 * vb / n_values};
      const double color_{color(value)};
      legend_series.set_value(value, 0, color_);
      // std::cerr << value << " " << color_ << std::endl;
    }
    for (uint64_t cb{0}; cb <= n_values; ++cb) {
      const double color_{1.0 * cb / n_values};
      color_legend_series.set_value(color_, 0, color_);
    }
    for (const double value : values) {
      hist.add_point(value);
      color_hist.add_point(color(value));
    }
  }

  // Add raw postscript to graph
  PSHeat & ps(const std::string & ps__) {
    graph.ps(ps__);
    return *this;
  }

 private:
  PSPage page;
  PSGraph graph;
  PSHSeries<double, uint64_t> hist;
  PSGraph legend;
  PSXYSSeries legend_series;
  PSHSeries<double, uint64_t> color_hist;
  PSGraph color_legend;
  PSXYSSeries color_legend_series;
};

template <class PSSeries>
void PSGraphT<PSSeries>::prepare_hist_overlap() {
  std::unique_ptr<PSHSeries<int, double>> overlap_{
    std::make_unique<PSHSeries<int, double>>(children())};
  manage(std::move(overlap_));
}

class PSShade : public PSXYSeries {
 public:
  using PSDoc = PSDocT<PSPage>;

  // factory
  template <class ... Args>
  static PSShade * create(Args && ... args) {
    PSShade * xyseries__{new PSShade{std::forward<Args>(args)...}};
    paa::ownp(xyseries__);
    return xyseries__;
  }

  PSShade() : PSXYSeries{} { }
  explicit PSShade(PSGraph & graph__,
                   const std::string & target_color__ = "1 0 0") :
      PSXYSeries{graph__}, target_color{target_color__} {}
  PSShade(PSDoc & doc, const std::string & title__,
          const std::string & target_color__ = "1 0 0") :
      PSXYSeries{doc, title__, Bounds{0.0}, Marker{}},
      target_color{target_color__} {}
  PSShade(PSPage & page, const std::string & title__,
          const std::string & target_color__ = "1 0 0") :
      PSXYSeries{page, title__, Bounds{0.0}, Marker{}},
      target_color{target_color__} {}
  PSShade(PSShade && other) = default;

  virtual ~PSShade() {
    if (VERBOSE) std::cerr << "Destroy PSShade " << std::endl;
    this->finalize_doc();
  }

  // Finalize
  virtual void finalize(PSDoc & doc, const Bounds & bounds,
                        const Bounds & range, const PSPart &) {
    const double scales[2]{bounds.xw() / range.xw(),
          bounds.yw() / range.yw()};

    // Make histogram of points / values
    const uint64_t n_points{x.size()};
    const uint64_t n_bins{static_cast<uint64_t>(sqrt(n_points))};
    const uint64_t n_bins_side{static_cast<uint64_t>(
        std::min(2000.0, std::max(500.0, pow(n_bins, 0.75))))};
    // std::cerr << "Side Bins " << n_bins_side << std::endl;
    std::vector<std::vector<uint64_t>> hist(n_bins_side,
                                            std::vector<uint64_t>(n_bins_side));
    for (uint64_t n{0}; n != n_points; ++n) {
      const uint64_t x_bin(n_bins_side * (x[n] - range.xl()) / range.xw());
      const uint64_t y_bin(n_bins_side * (y[n] - range.yl()) / range.yw());
      ++hist[x_bin][y_bin];
    }

    uint64_t minb{100000000000000000};
    uint64_t maxb{0};
    for (uint64_t xi{0}; xi != n_bins_side; ++xi) {
      for (uint64_t yi{0}; yi != n_bins_side; ++yi) {
        if (minb > hist[xi][yi]) minb = hist[xi][yi];
        if (maxb < hist[xi][yi]) maxb = hist[xi][yi];
      }
    }
    std::vector<double> xb;
    std::vector<double> yb;
    std::vector<double> values;
    for (uint64_t xi{0}; xi != n_bins_side; ++xi) {
      const double xv{range.xl() + xi * range.xw() / n_bins_side};
      for (uint64_t yi{0}; yi != n_bins_side; ++yi) {
        const double yv{range.yl() + yi * range.yw() / n_bins_side};
        if (hist[xi][yi]) {
          const double zv{1.0 * hist[xi][yi] / maxb};
          xb.push_back(xv);
          yb.push_back(yv);
          values.push_back(zv);
        }
      }
    }
    if (values.empty()) return;
    if (histeq_) hist_eq = std::make_unique<HistEq>(values, n_colors);

    std::function<double(double)> identity_color{
      [] (const double value) {
        return value;
      }};
    const auto minmax = minmax_element(values.begin(), values.end());
    std::function<double(double)> expanded_color{
      [minmax](const double value) {
        return (value - *minmax.first) / (*minmax.second - *minmax.first);
      }};
    std::function<double(double)> histeq_color{
        [this](const double value) {
          return 1.0 * ((*hist_eq)(value)) / n_colors;
        }};
    color = [this, identity_color, expanded_color, histeq_color]
        (const double value) {
      return bound((histeq_ ? histeq_color :
                    (expanded_ ? expanded_color : identity_color))(value),
                   0, 1);
    };
    const double x_width{range.xw() / n_bins_side};
    const double y_width{range.yw() / n_bins_side};
    std::istringstream color_stream{target_color.c_str()};
    double r;
    double g;
    double b;
    color_stream >> r >> g >> b;
    const double total_color{r + g + b};
    const double max_color{std::max(r, std::max(g, b))};
    const bool single_color{total_color >= 1 && total_color <= 1 &&
          max_color >= 1 && max_color <= 1};
    const unsigned int max_index{r > g ? (r > b ? 0U : 2U) : (g > b ? 1U : 2U)};

    doc << "gs "
        << bounds.xl() << " " << bounds.yl() << " tr "
        << scales[0] << " " << scales[1] << " scale "
        << range.xl() << " neg  " << range.yl() << " neg tr\n";
    doc << "/xsz " << x_width << " def "
        << "/ysz " << y_width << " def "
        << "/p { np ";
    if (single_color) {
      if (max_index == 0) {
        doc << "0 0";
      } else if (max_index == 1) {
        doc << "0 e 0";
      } else {
        doc << "0 e 0 e";
      }
    }
    doc << " c m "
        << "xsz 0 rl 0 ysz rl xsz neg 0 rl 0 ysz neg rl cp fp } def\n";
    for (unsigned int i{0}; i != xb.size(); ++i) {
      const double c{color(values[i])};
      doc << xb[i] << " " << yb[i] << " ";
      if (single_color) {
        doc << c << " p\n";
      } else {
        doc << r * c << " " << g * c << " " << b * c << " p\n";
      }
    }
    doc << "gr\n";
  }
  bool histeq() const { return histeq_; }
  PSShade & histeq(const bool histeq__) {
    histeq_ = histeq__;
    return *this;
  }
  bool expanded() const { return expanded_; }
  PSShade & expanded(const bool expanded__) {
    expanded_ = expanded__;
    return *this;
  }

 protected:
  std::function<double(double)> color{nullptr};
  const unsigned int n_colors{1000};
  std::unique_ptr<HistEq> hist_eq{nullptr};
  bool histeq_{true};
  bool expanded_{false};
  std::string target_color{"1 0 0"};
};

using PSGraph = PSXYSeries::PSGraph;
using PSPage = PSGraph::PSPage;
using PSDoc = PSPage::PSDoc;

#if 1

class MarkerGraph {
 public:
  explicit MarkerGraph(PSDoc & doc) :
      colors{"0.75 0.75 0.75", "0.5 0.5 0.5", "0.25 0.25 0.25", "0 0 0",
        "1 0 0", "0.75 0 0", "0.5 0 0", "0.25 0 0",
        "0 1 0", "0 0.75 0", "0 0.5 0", "0 0.25 0",
        "0 0 1", "0 0 0.75", "0 0 0.5", "0 0 0.25"},
    shapes{circle, triangle, itriangle, square, diamond, pentagon,
          hexagon, octagon, plus, xshape, star},
    graph{doc, "Markers;Shape;Color",
          Bounds{0.0, shapes.size() * 3 + 1.0, 0.0, colors.size() + 1.0}} {
      series.reserve(shapes.size() * 3);
      auto shape = shapes.begin();
      Marker marker{"", 1.5};
      for (unsigned int s{0}; s!= shapes.size(); ++s) {
        marker.path((*shape++)());
        for (const unsigned int fill : { 0, 1, 2 }) {
          series.emplace_back(graph);
          for (unsigned int c{0}; c != colors.size(); ++c) {
            marker.colors(colors[c]).fill(fill);
            if (fill == 2) marker.color("bk");
            series.back().add_point(3 * s + fill + 1, c + 1, marker);
          }
        }
      }
    }

 private:
  std::vector<std::string> colors{};
  std::vector<std::string (*)()> shapes{};
  PSGraph graph{};
  std::vector<PSXYMSeries> series{};
};

#endif
#if 0

// Test layout
class FormatTester {
 public:
  FormatTester() {
    std::vector<std::string> formats{"", "2 2", "2 2 =1 1 2 1=",
          "2 2 (0.3) +1 1 2 2+", "2 2 (0.3) +1 1 2 2+ +0 1 1 2 =0 1 2 1=+"};
    for (const std::string & format : formats) {
      std::cout << std::endl << "format: (" << format << ")" << std::endl;
      Layout layout{Bounds{0, 1, 0, 10}, format};
      layout.describe();
    }
    Layout layout{Bounds{0, 1, 0, 10}, formats.back()};
    Bounds bounds{layout.bounds(0, 1, 0, 1, 1, 0)};
    sout << bounds << " = " << layout[3].bounds() << std::endl;
  }
};
class LayoutTester {
 public:
  explicit LayoutTester(PSDoc & ps) :
      page{ps, "Layout Tester",
        "1 2 (0.3) =0 1 2 (0.6) 1 +1 0 2 2 -0 0 2 1- -1 1 2 2 |0 0 3 3|-+="},
    graph{page, "Layout;X;Y"},
    series{graph},
    shapes{circle, circle, circle, circle, triangle, itriangle, square,
          diamond, pentagon, hexagon, octagon, plus, xshape, star,
          circle, circle, circle} {
      for (unsigned int s{0}; s != 17; ++s) {
        series.add_point(s, s * s, Marker{shapes[s]()});
        page.add(&graph);
      }
    }

 private:
  PSPage page{};
  PSGraph graph{};
  PSSeries series{};
  std::vector<std::string (*)()> shapes{};
};

#endif



const std::string & get_color(const uint64_t index) {
  static const std::vector<std::string> distinct_colors{
    "0.898039 0 0", "0.145098 0 0.619608",
        "0 0.717647 0", "0.898039 0.745098 0",
        "0.0235294 0.337255 0.576471", "0.717647 0.866667 0",
        "0.898039 0.513725 0", "0.584314 0 0.584314",
        "0.972549 0.470588 0.972549", "0 0.0941176 0",
        "0 0.941176 0.533333", "0.564706 0.627451 0.533333",
        "0.972549 0.972549 0.627451", "0 0.658824 0.972549",
        "0.439216 0.313725 0.972549", "0.972549 0.0313725 0.972549",
        "0.470588 0.282353 0.188235", "0.972549 0.25098 0.470588",
        "0.470588 0.972549 0.376471", "0 0.156863 0.972549",
        "0.439216 0.596078 0", "0.12549 0.627451 0.376471",
        "0.972549 0.596078 0.470588", "0.627451 0.658824 0.972549",
        "0.282353 0.972549 0", "0.0627451 0.407843 0.0941176",
        "0.439216 0 0", "0.0313725 0.972549 0.909804",
        "0.376471 0 0.941176", "0.596078 0.313725 0.596078",
        "0.972549 0.784314 0.972549", "0.282353 0.784314 0.721569",
        "0.12549 0.156863 0.313725", "0.972549 0.972549 0.219608",
        "0.282353 0.501961 0.752941", "0.658824 0.878431 0.690196",
        "0.815686 0.25098 0.12549", "0.784314 0.25098 0.909804",
        "0 0.972549 0.219608", "0.878431 0 0.376471",
        "0.690196 0.470588 0.282353", "0.784314 0.784314 0.345098",
        "0 0.407843 0.909804", "0.345098 0.439216 0.439216",
        "0.282353 0.219608 0.721569", "0 0.627451 0.658824",
        "0.407843 0 0.313725", "0.376471 0.752941 0.25098",
        "0.784314 0.533333 0.721569", "0.784314 0.972549 0.972549",
        "0.690196 0 0.909804", "0.658824 0.156863 0.376471",
        "0.627451 0.407843 0", "0.313725 0.658824 0.972549",
        "0.847059 0.0941176 0.690196", "0.25098 0.25098 0",
        "0.878431 0.752941 0.658824", "0.376471 0.219608 0.439216",
        "0.972549 0.407843 0.25098", "0.407843 0.972549 0.658824",
        "0.658824 0.439216 0.941176", "0.690196 0.658824 0.12549",
        "0.658824 0 0.156863", "0.533333 0.156863 0.815686",
        "0.0627451 0.407843 0.345098", "0.25098 0.847059 0.470588",
        "0.188235 0.596078 0.12549", "0.188235 0 0.156863",
        "0.721569 0.972549 0.25098", "0 0.156863 0.721569",
        "0.972549 0.407843 0.658824", "0.470588 0.815686 0",
        "0 0.815686 0.752941", "0.596078 0.188235 0",
        "0.470588 0.564706 0.25098", "0.0941176 0.784314 0.219608",
        "0 0 0.407843", "0.313725 0.627451 0.564706",
        "0.533333 0.533333 0.752941", "0.0941176 0 0.847059",
        "0.25098 0.847059 0.972549", "0.0313725 0.909804 0",
        "0.784314 0.407843 0.501961", "0.25098 0.439216 0.972549",
        "0.909804 0.627451 0.219608", "0.470588 0.784314 0.878431",
        "0.313725 0.439216 0.0627451", "0.564706 0.815686 0.470588",
        "0.972549 0.12549 0.188235", "0.752941 0.972549 0.501961",
        "0.25098 0.941176 0.25098", "0.501961 0.972549 0.12549",
        "0.972549 0.627451 0.815686", "0.376471 0.0313725 0.690196",
        "0.25098 0.188235 0.941176", "0 0.752941 0.501961",
        "0.156863 0.972549 0.721569", "0.596078 0.815686 0.219608",
        "0.972549 0.847059 0.439216", "0.12549 0.501961 0.533333",
        "0.219608 0.313725 0.219608", "0.25098 0.752941 0.0313725",
        "0.156863 0.188235 0.533333", "0.972549 0.25098 0.784314",
        "0.972549 0.345098 0", "0.752941 0.784314 0.878431",
        "0.439216 0.345098 0.752941", "0.564706 0.345098 0.376471",
        "0.156863 0.658824 0.815686", "0.376471 0.12549 0.156863",
        "0.470588 0.533333 0.972549", "0.972549 0.972549 0.878431",
        "0 0.219608 0.156863", "0.752941 0.596078 0.439216",
        "0.878431 0.972549 0", "0.501961 0.156863 0.596078",
        "0.658824 0.690196 0.721569", "0 0.533333 0.219608",
        "0.784314 0.219608 0.533333", "0.784314 0.345098 0.721569",
        "0.12549 0.345098 0.752941", "0.815686 0.596078 0.972549",
        "0.972549 0.784314 0.219608", "0.282353 0.0627451 0.470588",
        "0.282353 0.345098 0.596078", "0.752941 0.313725 0.313725",
        "0.815686 0.972549 0.752941", "0.501961 0.470588 0.564706",
        "0.25098 0.470588 0.25098", "0.407843 0.658824 0.752941",
        "0 0.470588 0.690196", "0 0.815686 0.941176",
        "0.564706 0.690196 0.345098", "0.501961 0.784314 0.658824",
        "0.627451 0 0.376471", "0.219608 0.0941176 0",
        "0.972549 0 0.596078", "0.784314 0.12549 0",
        "0.847059 0.156863 0.345098", "0.313725 0.972549 0.847059",
        "0.627451 0.878431 0.972549", "0.815686 0.407843 0.12549",
        "0 0.533333 0", "0.0627451 0.878431 0.376471",
        "0.439216 0.313725 0", "0.219608 0.313725 0.407843",
        "0.345098 0.627451 0.376471", "0.752941 0.0627451 0.501961",
        "0 0.0627451 0.219608", "0 0.25098 0.407843",
        "0.0313725 0.282353 0", "0.941176 0.188235 0.972549",
        "0.627451 0.313725 0.815686", "0.721569 0.752941 0.533333",
        "0.847059 0.878431 0.156863", "0.690196 0.188235 0.721569",
        "0.596078 0.156863 0.188235", "0 0.0941176 0.564706",
        "0.501961 0.439216 0.156863", "0.972549 0.188235 0",
        "0.658824 0.156863 0.972549", "0.533333 0 0.784314",
        "0.721569 0.533333 0", "0.815686 0.0627451 0.188235",
        "0.501961 0.690196 0.12549", "0.0941176 0.282353 0.941176",
        "0.0941176 0.533333 0.941176", "0.439216 0.156863 0.972549",
        "0.564706 0.972549 0.564706", "0.564706 0.972549 0.815686",
        "0.658824 0.439216 0.658824", "0.12549 0.847059 0.627451",
        "0.721569 0 0.721569", "0.847059 0.376471 0.972549",
        "0.815686 0.0941176 0.941176", "0.972549 0.847059 0.784314",
        "0.25098 0.972549 0.564706", "0.219608 0.0941176 0.784314",
        "0.219608 0.752941 0.345098", "0.658824 0.313725 0.156863",
        "0.470588 0.0941176 0.439216", "0.188235 0.878431 0.0941176",
        "0.847059 0.501961 0.345098", "0.407843 0.752941 0.439216",
        "0.188235 0 0.345098", "0.282353 0.345098 0.847059",
        "0.188235 0.156863 0.156863", "0.376471 0.345098 0.313725",
        "0.878431 0.972549 0.376471", "0.596078 0.533333 0.407843",
        "0.376471 0.878431 0.156863", "0.972549 0.439216 0.439216",
        "0.909804 0.282353 0.282353", "0.941176 0 0.784314",
        "0.878431 0.721569 0.470588", "0 0.658824 0.156863",
        "0.313725 0.12549 0.313725", "0.188235 0.0627451 0.972549",
        "0.188235 0.721569 0.501961", "0.815686 0.658824 0.784314",
        "0.847059 0.878431 0.596078", "0.941176 0.12549 0.533333",
        "0.627451 0.752941 0", "0.847059 0.627451 0.596078",
        "0.407843 0.878431 0.533333", "0.188235 0.972549 0.407843",
        "0.878431 0.439216 0.815686", "0 0.533333 0.407843",
        "0 0.721569 0.345098", "0 0.345098 0.219608",
        "0.407843 0.156863 0", "0.972549 0.627451 0.0313725",
        "0 0.564706 0.815686", "0.439216 0.345098 0.533333",
        "0.878431 0.501961 0.564706", "0.407843 0.439216 0.878431",
        "0.156863 0.815686 0.847059", "0.501961 0.0313725 0.156863",
        "0.721569 0 0", "0.627451 0.909804 0.376471",
        "0.156863 0.690196 0.658824", "0.156863 0.721569 0.972549",
        "0.533333 0.0313725 0.972549", "0.752941 0.658824 0.282353",
        "0.407843 0.878431 0.784314", "0.690196 0.564706 0.878431",
        "0.627451 0.972549 0", "0.909804 0.282353 0.627451",
        "0.784314 0.533333 0.156863", "0.345098 0.501961 0.596078",
        "0.470588 0.188235 0.313725", "0.815686 0.752941 0.156863",
        "0.470588 0.847059 0.313725", "0.972549 0.847059 0.0627451",
        "0.470588 0.658824 0.909804", "0.470588 0.470588 0",
        "0.564706 0.564706 0.0941176", "0.596078 0.282353 0.972549",
        "0.972549 0 0.156863", "0.972549 0.12549 0.815686",
        "0 0 0.721569", "0 0 0.972549",
        "0.345098 0.627451 0.188235", "0 0.815686 0.0941176",
        "0.0313725 0.627451 0.501961", "0.690196 0.564706 0.596078",
        "0.721569 0.282353 0", "0.972549 0.313725 0.941176",
        "0.12549 0.219608 0.815686", "0.12549 0.470588 0.815686",
        "0.0313725 0 0.0941176", "0.188235 0.972549 0.972549",
        "0.376471 0.12549 0.815686", "0.407843 0 0.533333",
        "0.156863 0.564706 0.658824", "0.470588 0.972549 0.972549"
        };
  return distinct_colors[index % distinct_colors.size()];
}

const std::string & get_hex_color(const uint64_t index) {
  static const std::vector<std::string> distinct_colors{
    "#e50000", "#25009e", "#00b700", "#e5be00",
        "#065693", "#b7dd00", "#e58300", "#950095",
        "#f878f8", "#001800", "#00f088", "#90a088",
        "#f8f8a0", "#00a8f8", "#7050f8", "#f808f8",
        "#784830", "#f84078", "#78f860", "#0028f8",
        "#709800", "#20a060", "#f89878", "#a0a8f8",
        "#48f800", "#106818", "#700000", "#08f8e8",
        "#6000f0", "#985098", "#f8c8f8", "#48c8b8",
        "#202850", "#f8f838", "#4880c0", "#a8e0b0",
        "#d04020", "#c840e8", "#00f838", "#e00060",
        "#b07848", "#c8c858", "#0068e8", "#587070",
        "#4838b8", "#00a0a8", "#680050", "#60c040",
        "#c888b8", "#c8f8f8", "#b000e8", "#a82860",
        "#a06800", "#50a8f8", "#d818b0", "#404000",
        "#e0c0a8", "#603870", "#f86840", "#68f8a8",
        "#a870f0", "#b0a820", "#a80028", "#8828d0",
        "#106858", "#40d878", "#309820", "#300028",
        "#b8f840", "#0028b8", "#f868a8", "#78d000",
        "#00d0c0", "#983000", "#789040", "#18c838",
        "#000068", "#50a090", "#8888c0", "#1800d8",
        "#40d8f8", "#08e800", "#c86880", "#4070f8",
        "#e8a038", "#78c8e0", "#507010", "#90d078",
        "#f82030", "#c0f880", "#40f040", "#80f820",
        "#f8a0d0", "#6008b0", "#4030f0", "#00c080",
        "#28f8b8", "#98d038", "#f8d870", "#208088",
        "#385038", "#40c008", "#283088", "#f840c8",
        "#f85800", "#c0c8e0", "#7058c0", "#905860",
        "#28a8d0", "#602028", "#7888f8", "#f8f8e0",
        "#003828", "#c09870", "#e0f800", "#802898",
        "#a8b0b8", "#008838", "#c83888", "#c858b8",
        "#2058c0", "#d098f8", "#f8c838", "#481078",
        "#485898", "#c05050", "#d0f8c0", "#807890",
        "#407840", "#68a8c0", "#0078b0", "#00d0f0",
        "#90b058", "#80c8a8", "#a00060", "#381800",
        "#f80098", "#c82000", "#d82858", "#50f8d8",
        "#a0e0f8", "#d06820", "#008800", "#10e060",
        "#705000", "#385068", "#58a060", "#c01080",
        "#001038", "#004068", "#084800", "#f030f8",
        "#a050d0", "#b8c088", "#d8e028", "#b030b8",
        "#982830", "#001890", "#807028", "#f83000",
        "#a828f8", "#8800c8", "#b88800", "#d01030",
        "#80b020", "#1848f0", "#1888f0", "#7028f8",
        "#90f890", "#90f8d0", "#a870a8", "#20d8a0",
        "#b800b8", "#d860f8", "#d018f0", "#f8d8c8",
        "#40f890", "#3818c8", "#38c058", "#a85028",
        "#781870", "#30e018", "#d88058", "#68c070",
        "#300058", "#4858d8", "#302828", "#605850",
        "#e0f860", "#988868", "#60e028", "#f87070",
        "#e84848", "#f000c8", "#e0b878", "#00a828",
        "#502050", "#3010f8", "#30b880", "#d0a8c8",
        "#d8e098", "#f02088", "#a0c000", "#d8a098",
        "#68e088", "#30f868", "#e070d0", "#008868",
        "#00b858", "#005838", "#682800", "#f8a008",
        "#0090d0", "#705888", "#e08090", "#6870e0",
        "#28d0d8", "#800828", "#b80000", "#a0e860",
        "#28b0a8", "#28b8f8", "#8808f8", "#c0a848",
        "#68e0c8", "#b090e0", "#a0f800", "#e848a0",
        "#c88828", "#588098", "#783050", "#d0c028",
        "#78d850", "#f8d810", "#78a8e8", "#787800",
        "#909018", "#9848f8", "#f80028", "#f820d0",
        "#0000b8", "#0000f8", "#58a030", "#00d018",
        "#08a080", "#b09098", "#b84800", "#f850f0",
        "#2038d0", "#2078d0", "#080018", "#30f8f8",
        "#6020d0", "#680088", "#2890a8", "#78f8f8"
        };
  return distinct_colors[index % distinct_colors.size()];
}

std::string eps(const std::string & plot_file_name,
                const std::string & translate,
                const double scale) {
  std::ifstream plot_file(plot_file_name.c_str());
  if (!plot_file) throw Error("Could not open plot file") << plot_file_name;
  std::string line;
  std::ostringstream out;

  out << "% eps start for " << plot_file_name << "\n";
  out << "save" << "\n";
  out << translate << " translate "
      << scale << " " << scale << " scale" << "\n";
  while (getline(plot_file, line)) {
    if (line.size() && line[0] == '%') continue;
    out << line << "\n";
  }
  out << "restore" << "\n";
  out << "% eps stop for " << plot_file_name << "\n";
  return out.str();
}

template <class Value>
void plot_continuous(const std::vector<Value> & values,
                     PSGraph & graph, const std::string & color,
                     const double x_thresh = 0.0001,
                     const double y_thresh = 0.0001,
                     const double factor = 1.0) {
  PSXYSeries & series{*PSXYSeries::create(
      graph, Marker{paa::circle(), 0.3, color, 0.2, true})};
  double last_x{-100};
  Value last_y{10 * values.back()};
  for (uint64_t v{0}; v != values.size(); ++v) {
    const double x{1.0 * v / values.size()};
    const Value y{values[v]};
    if (fabs(x - last_x) > x_thresh || fabs(y - last_y) > y_thresh) {
      series.add_point(x, factor * y);
      last_x = x;
      last_y = y;
    }
  }
}

void draw_vertical_line(PSGraph & graph, const double x,
                        const std::string & color,
                        const std::string & dash = "") {
  std::ostringstream ps;
  if (dash.size()) ps << dash << " sd ";
  ps << "1 lw " << color << " c np " << x << " xc 0 yfc m "
     << x << " xc 1 yfc l sp";
  if (dash.size()) ps << " nd ";
  ps << "\n";
  graph.ps(ps.str());
}
void draw_horizontal_line(PSGraph & graph, const double y,
                          const std::string & color) {
  std::ostringstream ps;
  ps << "1 lw " << color << " c np 0 xfc " << y << " yc m 1 xfc "
     << y << " yc l sp\n";
  graph.ps(ps.str());
}

void ur_legend_box(PSGraph & graph, const std::vector<std::string> & text) {
  std::ostringstream ps;
  ps << "10 sf /xr 0.97 xfc def /yt 0.97 yfc def /ls 13 def /pd 4 def\n"
     << "/max { 2 copy lt {exch pop} {pop} ifelse } def\n"
     << "/mw 0 def\n";
  for (uint64_t l{0}; l != text.size(); ++l) {
    const std::string & line{text[l]};
    ps << "/l" << l << " (" << line << ") def"
       << " l" << l << " sw pop mw max /mw e def\n";
  }
  ps << "/yb yt 2 pd mul " << text.size() << " ls mul add sub def\n"
     << "/xl xr 2 pd mul mw add sub def 0 0 0 c\n"
     << "1 lw np xl yb m xr yb l xr yt l xl yt l cp gs 1 1 1 c fp gr sp\n"
     << "xl pd add yt 1 sub m\n";
  for (uint64_t l{0}; l != text.size(); ++l) {
    ps << "0 ls neg rm gs l" << l << " s gr\n";
  }
  graph.ps(ps.str());
}
void ur_legend_box(PSGraph & graph, const std::string & text) {
  std::istringstream text_stream{text.c_str()};
  std::vector<std::string> lines;
  std::string line;
  while (getline(text_stream, line)) lines.push_back(line);
  ur_legend_box(graph, lines);
}

PSGraph & typed_hist(PSPage & page, const std::string & title,
                     const std::string & y_axis_name,
                     const std::vector<double> & values,
                     const std::vector<std::string> & names,
                     const std::vector<std::string> & colors = {},
                     const bool show_value = false,
                     const double max_value = 100,
                     const double fraction_empty = 0.3) {
  const double gap_size{fraction_empty / (values.size() + 1)};
  const double bar_size{(1 - fraction_empty) / values.size()};
  std::ostringstream ps;
  ps << "1 lw 0 0 0 c\n";
  for (uint64_t t{0}; t != values.size(); ++t) {
    const double value{values[t]};
    const double start{gap_size + t * (gap_size + bar_size)};
    const double stop{start + bar_size};
    const double mid{(start + stop) / 2};
    const std::string & color{t < colors.size() ? colors[t] : "0.9 0 0"};
    ps << "np " << start << " 0 gc m "
         << stop << " 0 gc l "
         << stop << " " << value << " gc l "
         << start << " " << value << " gc l cp "
       << "gs " << color << " c fp gr sp\n";
    if (names.size()) ps << "0 0 0 c 10 sf " << mid << " xc 0 yfc 12 sub m ("
                         << names[t] << ") jc s\n";
    if (show_value)
      ps << "0 0 0 c 10 sf " << mid << " xc " << value << " yc 4 add m ("
         << std::setprecision(4) << value
         << std::setprecision(6) << "%) jc s\n";
  }
  PSGraph & graph{*PSGraph::create(
      page, title + ";" + (names.size() ? "." : "") + ";" + y_axis_name,
      Bounds{0, 1, 0, max_value})};
  graph.x_label().color("1 1 1");
  graph.do_x_ticks(false);
  graph.ps(ps.str());
  return graph;
}

constexpr uint64_t lll_font{14};  // Lower left legend
constexpr double lll_y{1.0 * lll_font};
constexpr double lll_off{1.9 * lll_font};

}  // namespace paa

#endif  // PAA_PSPLOT_H

