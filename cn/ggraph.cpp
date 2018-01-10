//
// ggraph.cpp
//
// CN view gui
//
// Copyright 2016 Peter Andrews @ CSHL
//

// example command
// ggraph cn chrAll.fa abspos,ratio,seg 0,1 {m,f,d,s}.txt

#include <algorithm>
#include <exception>
#include <iostream>
#include <fstream>
#include <functional>
#include <future>
#include <map>
#include <memory>
#include <mutex>
#include <sstream>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include "cn.h"
#include "error.h"
#include "files.h"
#include "genes.h"
#include "mumdex.h"
#include "strings.h"
#include "tsv.h"
#include "threads.h"
#include "x11plot.h"

using std::cin;
using std::cout;
using std::cerr;
using std::endl;
using std::exception;
using std::ifstream;
using std::future;
using std::istringstream;
using std::make_unique;
using std::map;
using std::max;
using std::min;
using std::mutex;
using std::ostringstream;
using std::pair;
using std::placeholders::_1;
using std::placeholders::_2;
using std::string;
using std::unique_lock;
using std::unique_ptr;
using std::vector;

using paa::rect;
using paa::remove_substring;
using paa::remove_including_final;
using paa::Columns;
using paa::Bin;
using paa::iBounds;
using paa::ChromosomeIndexLookup;
using paa::CN_abspos;
using paa::Error;
using paa::Event;
using paa::GeneXrefs;
using paa::KnownGene;
using paa::KnownGenes;
using paa::Point;
using paa::Radio;
using paa::Reference;
using paa::ThreadPool;
using paa::TSV;
using paa::X11App;
using paa::X11Font;
using paa::X11Graph;
using paa::X11Win;

unsigned int n_threads{std::thread::hardware_concurrency()};

const unsigned int max_gene_mb{100};
const unsigned int max_name_mb{5};

using Rect = vector<vector<int>>;
class GeneInfo {
 public:
  GeneInfo(const string & name_, const string & description_,
           const Rect & bounds_, const int low_, const int high_) :
      name{name_}, description{description_}, bounds{bounds_},
    low{low_}, high{high_} { }

  const string name;
  const string description;
  Rect bounds;
  int low;
  int high;
};

struct BestVariant {
  BestVariant(const string description_, const int start_, const int stop_) :
      description{description_}, start{start_}, stop{stop_},
    low{start}, high{stop} { }
  string description;
  int start;
  int stop;
  int low;
  int high;
  int middle() const { return (start + stop) / 2; }
  int width() const { return stop - start; }
};

bool add_genes(const Reference & ref, const CN_abspos & cn_abspos,
               X11Graph & graph, const Event & event = Event()) {
  static bool first_call{true};
  if (first_call && graph.range[0][2] > max_gene_mb * 1000000.0) return false;
  static mutex genes_mutex;
  unique_lock<mutex> lock{genes_mutex};
  static const KnownGenes genes{[&ref]() {
      try {
        const string & reference_file{ref.fasta_file()};
        const ChromosomeIndexLookup chr_lookup{ref};
        const string genes_name{reference_file + ".bin/knownGene.txt"};
        const string isoforms_name{reference_file + ".bin/knownIsoforms.txt"};
        return KnownGenes{genes_name, isoforms_name, chr_lookup, ref};
      } catch (Error & err) {
        cerr << err.what() << endl;
        return KnownGenes{ref};
      }
    }()};
  static bool missing_message_shown = false;
  if (genes.size() == 0) {
    if (!missing_message_shown) {
      cerr << "******************************************************" << endl;
      cerr << "No gene information was loaded" << endl;
      cerr << "Please download the files" << endl;
      cerr << "knownGenes.txt, knownIsoforms.txt and kgXref.txt" << endl;
      cerr << "from the UCSC website, from mumdex.com or from the" << endl;
      cerr << "G-Graph paper supplementary materials" << endl;
      cerr << "to match your reference, and place in the directory" << endl;
      cerr << ref.fasta_file() + ".bin/" << endl;
      cerr << "Only hg19 and hg38 versions of these files were tested" << endl;
      cerr << "Note chromosome names must match your reference" << endl;
      cerr << "G-Graph will continue to run without this gene info" << endl;
      cerr << "******************************************************" << endl;
      missing_message_shown = true;
    }
    return false;
  }
  static const string & reference_file{ref.fasta_file()};
  static const string kgXrefs_name{reference_file + ".bin/kgXref.txt"};
  static const GeneXrefs xref{kgXrefs_name};
  if (first_call) {
    first_call = false;
    return false;  // Enable pre-loading in separate thread
  }
  lock.unlock();

  static vector<GeneInfo> gene_info;

  if (graph.log_radios[0]) return false;

  if (event.type == Event::X) {
    switch (event.x->type) {
      case MotionNotify:
        {
          // Show gene description on hover
          const XMotionEvent & xmotion{event.x->xmotion};
          for (const GeneInfo & gene : gene_info) {
            if (xmotion.x >= max(graph.bounds[0][0], gene.bounds[0][0]) &&
                xmotion.x <= min(graph.bounds[0][1], gene.bounds[0][1]) &&
                xmotion.y >= max(graph.bounds[1][0], gene.bounds[1][0]) &&
                xmotion.y <= min(graph.bounds[1][1], gene.bounds[1][1])) {
              graph.status = gene.description;
              remove_substring(&graph.status, "Homo sapiens ");
              remove_substring(&graph.status,
                               std::string(" (") + gene.name + ")");
              remove_substring(&graph.status, std::string(" ") + gene.name);
              graph.draw_status(true);
              graph.status = "";
              return true;
            }
          }
          // Show chromosomal position
          if (graph.coord_radio) {
            if (xmotion.x >= graph.bounds[0][0] &&
                xmotion.x <= graph.bounds[0][1] &&
                xmotion.y >= graph.bounds[1][0] &&
                xmotion.y <= graph.bounds[1][1]) {
              const pair<unsigned int, unsigned int> chrpos{
                cn_abspos.chrpos(graph.icoord(0, xmotion.x))};
              std::ostringstream coordinates;
              coordinates << std::setprecision(12) << "("
                          << ref.name(chrpos.first);
              const double vals[2]{static_cast<double>(chrpos.second),
                    graph.icoord(1, xmotion.y)};
              for (const bool y : {false, true}) {
                const double val{vals[y]};
                const double res{graph.range[y][2] / graph.bounds[y][2]};
                const double pres{pow(10, floor(log10(res)))};
                const double rval{round(val / pres) * pres};
                const double nval{graph.log_radios[y] ? pow(10, rval) : rval};
                coordinates << (y ? " , " : " ") << nval;
              }
              coordinates << " )";
              graph.status = coordinates.str();
              graph.draw_status();
              return true;
            }
          }
        }
        break;
      case ButtonRelease:
        {
          // Open genome browser
          const XButtonEvent & xbutton{event.x->xbutton};
          for (const GeneInfo & gene : gene_info) {
            if (xbutton.x >= gene.bounds[0][0] &&
                xbutton.x <= gene.bounds[0][1] &&
                xbutton.y >= gene.bounds[1][0] &&
                xbutton.y <= gene.bounds[1][1]) {
              ostringstream firefox;
              const pair<unsigned int, unsigned int> chrpos_start{
                cn_abspos.chrpos(graph.icoord(0, gene.low))};
              const pair<unsigned int, unsigned int> chrpos_stop{
                cn_abspos.chrpos(graph.icoord(0, gene.high))};
              const unsigned int chr{chrpos_start.first};
              const unsigned int range{
                chrpos_stop.second - chrpos_start.second};
              const unsigned int display_start(
                  chrpos_start.second > 0.1 * range + 1 ?
                  chrpos_start.second - 0.1 * range + 1 : 1);
              const unsigned int display_stop(
                  chrpos_stop.second + 0.1 * range < ref.size(chr) ?
                  chrpos_stop.second + 0.1 * range : ref.size(chr));
#ifdef __linux__
              const bool mac{false};
#else
              const bool mac{true};
#endif
              firefox << std::string(!mac ? "firefox" : "open -a safari")
                      << " 'http://genome-mirror.cshl.edu/"
                      << "cgi-bin/hgTracks?db=hg19&"
                      << "position=" << ref.name(chr) << ":"
                      << display_start << "-" << display_stop << "' &";
              if (system(firefox.str().c_str()) == -1) {
                std::cerr << "Problem starting browser" << std::endl;
              }
              return true;
            }
          }
        }
        break;

      default:
        break;
    }
  }

  if (event.type != Event::Draw) return false;

  gene_info.clear();
  if (graph.range[0][2] > max_gene_mb * 1000000.0) return false;

  // Get chrpos from range of graph
  unsigned int chr[2];
  unsigned int pos[2];
  for (const bool high : {false, true}) {
    const unsigned int abspos(max(0.0, graph.range[0][high]));
    const pair<unsigned int, unsigned int> chrpos{cn_abspos.chrpos(abspos)};
    chr[high] = chrpos.first;
    pos[high] = chrpos.second;
  }

  // Get genes in range
  std::vector<KnownGene>::const_iterator gene_limits[2];
  for (const bool high : {false, true}) {
    const unsigned int e{500000};
    const unsigned int chrpos[2]{chr[high],
          high ? pos[high] + e : (pos[high] > e ? pos[high] - e : 0)};
    gene_limits[high] = upper_bound(
        genes.begin(), genes.end(), chrpos,
        [high](const unsigned int * cp, const KnownGene & gene) {
          if (gene.chr == cp[0]) {
            return (high ? gene.t_stop : gene.t_start) >= cp[1];
          }
          return gene.chr >= cp[0];
        });
  }
  if (gene_limits[0] > gene_limits[1]) gene_limits[0] = gene_limits[1];

  // Draw all exons and genes in range
  const int exon_height{10};
  const unsigned int gene_y(graph.bounds[1][0] * 1.2);
  const unsigned int exon_y{gene_y - exon_height / 2};
  using Item = pair<string, BestVariant>;
  map<string, BestVariant> names;
  for (std::vector<KnownGene>::const_iterator g{gene_limits[0]};
       g != gene_limits[1]; ++g) {
    // Gene line
    const KnownGene & gene{*g};
    const double min_gene{max(graph.range[0][0],
                              1.0 * cn_abspos(gene.chr, gene.t_start))};
    const double max_gene{min(graph.range[0][1],
                              1.0 * cn_abspos(gene.chr, gene.t_stop))};
    if (max_gene <= graph.range[0][0] || min_gene >= graph.range[0][1])
      continue;
    const int gene_start{graph.coord(0, min_gene)};
    const int gene_stop{max(graph.coord(0, max_gene), gene_start + 1)};
    if (graph.range[0][2] <= max_gene_mb * 1000000.0) {
      const pair<string, BestVariant> to_add{
        xref[gene.name].geneSymbol,
        {xref[gene.name].description, gene_start, gene_stop}};
      auto item = names.insert(to_add);
      if (!item.second) {
        BestVariant & inside{(*item.first).second};
        if (inside.low > to_add.second.low)
          inside.low = to_add.second.low;
        if (inside.high < to_add.second.high)
          inside.high = to_add.second.high;
        if (inside.width() < to_add.second.width()) {
          inside.start = to_add.second.start;
          inside.stop = to_add.second.stop;
          inside.description = to_add.second.description;
        }
      }
    }
    XDrawLine(graph.display(), graph.window, graph.border_gc,
              gene_start, gene_y,
              gene_stop, gene_y);

    // Exon boxes
    for (unsigned int e{0}; e != gene.exon_starts.size(); ++e) {
      const double exon_start(cn_abspos(gene.chr, gene.exon_starts[e]));
      const double exon_stop(cn_abspos(gene.chr, gene.exon_stops[e]));
      if (exon_start >= graph.range[0][1] || exon_stop <= graph.range[0][0])
        continue;
      const double mod_start{max(graph.range[0][0], exon_start)};
      const double mod_stop{min(graph.range[0][1], exon_stop)};
      const int box_start{graph.coord(0, mod_start)};
      const double frac_pixel{(mod_stop - mod_start) * graph.scale[0]};
      if (frac_pixel < 0.5) continue;
      const int box_stop{frac_pixel < 1.5 ? box_start + 1 :
            graph.coord(0, mod_stop)};
      XFillRectangle(graph.display(), graph.window, graph.border_gc,
                     box_start, exon_y, box_stop - box_start, exon_height);
    }
  }

  if (graph.range[0][2] > max_name_mb * 1000000.0) return false;

  // Get good font for gene names
  const double avail_height{(graph.bounds[1][0] - graph.border_width) * 0.6};
  static const X11Font * last_font{nullptr};
  static GC gc{[&graph]() {
      GC gc_{XCreateGC(graph.display(), graph.window, 0, nullptr)};
      XSetForeground(graph.display(), gc_, graph.app.black);
      return gc_;
    }()};
  const X11Font * fits{graph.app.fonts.fits("A", 1000, avail_height)};
  if (fits != last_font) XSetFont(graph.display(), gc, fits->id());
  last_font = fits;

  // Set clip rectangle
  static iBounds last_bounds;
  if (graph.bounds != last_bounds) {
    // if (paa::bne(graph.bounds, last_bounds)) {
    XRectangle clip_rectangle(rect(graph.bounds));
    XSetClipRectangles(graph.display(), gc, 0, 0, &clip_rectangle, 1, YXBanded);
  }

  // Draw gene names nicely
  std::vector<Item> snames(names.begin(), names.end());
  sort(snames.begin(), snames.end(), [](const Item & lhs, const Item & rhs) {
      const int l{lhs.second.low + lhs.second.high};
      const int r{rhs.second.low + rhs.second.high};
      if (l == r) return lhs.first < rhs.first;
      return l < r;
    });
  vector<int> last_right(100);
  for (const Item & entry : snames) {
    const string & name{entry.first};

    if (false && name == "ZNF462") {
      const double del_start{cn_abspos(chr[0], 109694531) + 0.0};
      const double del_stop{cn_abspos(chr[0], 109700021) + 0.0};
      const int box_start{graph.coord(0, del_start)};
      const int box_stop{graph.coord(0, del_stop)};
      const int height{graph.coord(1, 1.25)};
      XFillRectangle(graph.display(), graph.window, graph.border_gc,
                     box_start, height, box_stop - box_start, 10);
      const std::string descr{
        "SSC02971 proband 5492 base deletion chr9 109694531 - 109700021"};
      const int xpos{fits->centered_x(descr, (box_start + box_stop) / 2)};
      XDrawString(graph.display(), graph.window, gc, xpos,
                    height - 20, const_cast<char *>(descr.c_str()),
                    static_cast<unsigned int>(descr.size()));
    }

    const int tpos{entry.second.middle()};
    const int width(fits->string_width(name + "  "));
    const int left(tpos - width / 2);
    const int xpos{fits->centered_x(name, tpos)};
    for (unsigned int i{0}; i != last_right.size(); ++i) {
      if (left > last_right[i]) {
        const int ypos(exon_y + 2 * graph.border_width +
                       exon_height + i * fits->height());
        XDrawString(graph.display(), graph.window, gc, xpos,
                    fits->below_y(ypos), const_cast<char *>(name.c_str()),
                    static_cast<unsigned int>(name.size()));
        const int right{tpos + width / 2};
        last_right[i] = right;
        const Rect gene_rect{{left, right},
          {ypos, ypos + fits->height()}};
        gene_info.emplace_back(name, entry.second.description, gene_rect,
                               entry.second.low, entry.second.high);
        break;
      }
    }
  }

  // Show sequence!
  const string test_bases{"GGGGGGGGGG"};
  const double base_width{graph.bounds[0][2] / graph.range[0][2]};
  const X11Font * bfont{graph.app.fonts.fits(
      test_bases, test_bases.size() * base_width, 1000)};
  const int actual_width(bfont->string_width(test_bases) / test_bases.size());
  if (actual_width <= base_width) {
    if (bfont != last_font) XSetFont(graph.display(), gc, bfont->id());
    last_font = bfont;
    unsigned int n_bases{0};
    const int y_pos{graph.bounds[1][1] - 10 * graph.border_width};
    for (unsigned int b = max(0.0, graph.range[0][0]);
         b < graph.range[0][1]; ++b) {
      const pair<unsigned int, unsigned int> chrpos{cn_abspos.chrpos(b)};
      const string base{ref[chrpos.first][chrpos.second]};
      const int x_coord{graph.coord(0, b)};
      const int x_pos{bfont->centered_x(base, x_coord)};
      if (x_coord - actual_width / 2 > graph.bounds[0][0] &&
          x_coord + actual_width / 2 < graph.bounds[0][1]) {
        XDrawString(graph.display(), graph.window, gc,
                    x_pos, y_pos, base.c_str(), 1);
        ++n_bases;
      }
    }
    // Stupid easter egg "it's turtles all the way down"
    if (n_bases <= 1) {
      const string turtle{"turtle"};
      const X11Font * tfont{graph.app.fonts.fits(
          turtle, graph.bounds[0][2] / 10, 1000)};
      if (tfont != last_font) XSetFont(graph.display(), gc, tfont->id());
      last_font = tfont;
      const paa::Axis axis{graph.range[0][0], graph.range[0][1], 8, false};
      for (const std::pair<double, bool> tick : axis.ticks()) {
        if (!tick.second) continue;
        const int x_coord{graph.coord(0, tick.first)};
        const int x_pos{bfont->centered_x(turtle, x_coord)};
        XDrawString(graph.display(), graph.window, gc,
                    x_pos, y_pos - 1.5 * bfont->height(),
                    turtle.c_str(), static_cast<unsigned int>(turtle.size()));
      }
    }
  }
  return false;
}

bool add_chromosomes(const Reference & ref, const CN_abspos & cn_abspos,
                     X11Graph & graph, const Event & event = Event()) {
  if (event.type != Event::Draw) return false;

  const vector<unsigned int> & chromosomes{cn_abspos.chromosomes()};
  int min_w{100000000};
  vector<unsigned int> chr;
  vector<int> pos;
  vector<int> widths;
  for (unsigned int c{0}; c != chromosomes.size(); ++c) {
    double cl(cn_abspos.cn_offset(c));
    double ch(cn_abspos.cn_offset(c + 1));
    if (graph.log_radios[0]) cl = log10(cl);
    if (graph.log_radios[0]) ch = log10(ch);
    const double gl{graph.range[0][0]};
    const double gh{graph.range[0][1]};
    if (cl < gh && ch > gl) {
      const double l{max(cl, gl)};
      const double h{min(ch, gh)};
      const int tl{graph.coord(0, l)};
      const int th{graph.coord(0, h)};
      const int m{(th + tl) / 2};
      const int w{th - tl};
      min_w = min(min_w, w);
      widths.push_back(w);
      chr.push_back(chromosomes[c]);
      pos.push_back(m + 1);
    }
  }

  // Chromosome name
  const double avail_height{(graph.bounds[1][0] - graph.border_width) * 0.6};
  const std::string longest{"22"};
  static X11Font * last_font{nullptr};
  static GC gc{[&graph]() {
      GC gc_{XCreateGC(graph.display(), graph.window, 0, nullptr)};
      XSetForeground(graph.display(), gc_, graph.app.black);
      return gc_;
    }()};
  const X11Font * fits{graph.app.fonts.fits(longest, 1000, avail_height)};
  if (fits != last_font) XSetFont(graph.display(), gc, fits->id());

  static iBounds last_bounds;
  if (graph.bounds != last_bounds) {
    XRectangle clip_rectangle(rect(graph.bounds));
    XSetClipRectangles(graph.display(), gc, 0, 0, &clip_rectangle, 1, YXBanded);
  }
  for (unsigned int c{0}; c != chr.size(); ++c) {
    const std::string name{remove_substring(ref.name(chr[c]), "chr")};
    if (fits->string_width(name) < widths[c]) {
      XDrawString(graph.display(), graph.window, gc, fits->centered_x(
          name, pos[c]), fits->centered_y(
              graph.bounds[1][1] - graph.border_width - avail_height / 2),
                  const_cast<char *>(name.c_str()),
                  static_cast<unsigned int>(name.size()));
    }
  }

  // Chromosome boundary lines
  for (unsigned int c{0}; c <= chromosomes.size(); ++c) {
    double chr_start(cn_abspos.cn_offset(c));
    if (graph.log_radios[0]) chr_start = log10(chr_start);
    if (chr_start > graph.range[0][0] && chr_start < graph.range[0][1]) {
      const unsigned int x_pos(graph.coord(0, chr_start));
      XDrawLine(graph.display(), graph.window, graph.gc,
                x_pos, graph.bounds[1][0], x_pos, graph.bounds[1][1]);
    }
  }
  return false;
}

bool add_ratio_lines(const vector<double> cn_lines,
                     X11Graph & graph, const Event & event = Event()) {
  if (event.type != Event::Draw) return false;

  // Ratio lines
  for (double y : cn_lines) {
    if (graph.log_radios[1]) {
      if (y <= 0) continue;
      y = log10(y);
    }
    if (y > graph.range[1][0] && y < graph.range[1][1]) {
      const unsigned int y_pos(graph.coord(1, y));
      XDrawLine(graph.display(), graph.window, graph.gc,
                graph.bounds[0][0], y_pos, graph.bounds[0][1], y_pos);
    }
  }
  return false;
}

struct CytobandInfo {
  CytobandInfo(const string & chr,
               const unsigned int start, const unsigned int stop,
               const string & name_, const GC gc_) :
      abspos_start{start}, abspos_stop{stop},
    name{chr + " " + name_}, gc{gc_} { }
  CytobandInfo(const int start, const int stop,
               const string & name_) :
      abspos_start{static_cast<unsigned int>(start)},
    abspos_stop{static_cast<unsigned int>(stop)}, name{name_} { }
  CytobandInfo(const CytobandInfo &) = default;
  CytobandInfo & operator=(const CytobandInfo &) = default;
  unsigned int abspos_start;
  unsigned int abspos_stop;
  string name;
  GC gc{nullptr};
  bool operator<(const CytobandInfo & rhs) const {
    return abspos_start < rhs.abspos_start;
  }
};

bool add_cytobands(const Reference & ref, const CN_abspos & cn_abspos,
                   X11Graph & graph, const Event & event = Event()) {
  static vector<CytobandInfo> cytobands{[&ref, &cn_abspos, &graph]() {
      // GCs for different stain colors
      const map<string, string> stain_colors{
        {"gneg", "grey100"},
        {"gpos25", "grey90"},
        {"gpos50", "grey70"},
        {"gpos75", "grey40"},
        {"gpos100", "grey0"},
        {"gvar", "grey100"},
        {"stalk", "brown3"},
        {"acen", "brown4"}};
      map<string, GC> stain_gcs;
      for (const auto sc : stain_colors) {
        const string stain_name{sc.first};
        const string color_name{sc.second};
        auto inserted = stain_gcs.emplace(stain_name, nullptr);
        if (inserted.second) {
          GC & gc{inserted.first->second};
          XColor color;
          if (!XAllocNamedColor(graph.display(), graph.app.colormap,
                                stain_colors.at(stain_name).c_str(),
                                &color, &color))
            throw Error("Could not get color") << color_name;
          gc = graph.create_gc(color.pixel, graph.app.white, 1,
                               LineSolid, CapButt, JoinMiter);
        }
      }

      // Read bands file
      const string bands_name{ref.fasta_file() + ".bin/cytoBand.txt"};
      ifstream bands{bands_name.c_str()};
      if (!bands) {
        cerr << "Could not load cytobands file " << bands_name << endl;
        return vector<CytobandInfo>{};
      }
      // bands.ignore(1000, '\n');
      string chr_name;
      unsigned int start;
      unsigned int stop;
      string name;
      string stain;
      vector<CytobandInfo> result;
      const ChromosomeIndexLookup chr_lookup{ref};
      while (bands >> chr_name >> start >> stop) {
        bands.get();
        getline(bands, name, '\t');
        bands >> stain;
        if (!bands) throw Error("Problem parsing cytoband file") << bands_name;
        if (chr_name.find("_") != string::npos ||
            chr_name.find("M") != string::npos) {
          continue;
        }
        const unsigned int chr{chr_lookup[chr_name]};
        result.emplace_back(chr_name,
                            cn_abspos(chr, start), cn_abspos(chr, stop),
                            name, stain_gcs[stain]);
      }
      sort(result.begin(), result.end());
      return result;
    }()};

  static bool shown_missing_message{false};
  if (cytobands.empty()) {
    if (!shown_missing_message) {
      cerr << "******************************************************" << endl;
      cerr << "No cytoband information was loaded" << endl;
      cerr << "Please download the file cytoBand.txt" << endl;
      cerr << "from the UCSC website, from mumdex.com or from the" << endl;
      cerr << "G-Graph paper supplementary materials" << endl;
      cerr << "to match your reference, and place in the directory" << endl;
      cerr << ref.fasta_file() + ".bin/" << endl;
      cerr << "Only hg19 and hg38 versions of these files were tested" << endl;
      cerr << "Note chromosome names must match your reference" << endl;
      cerr << "G-Graph will continue to run without this info" << endl;
      cerr << "******************************************************" << endl;
      shown_missing_message = true;
    }
    return false;
  }

  static int band_low{0};
  static int band_high{0};
  static vector<CytobandInfo> band_info;

  if (event.type == Event::X && event.x->type == MotionNotify) {
    // Show gene description on hover
    const XMotionEvent & xmotion{event.x->xmotion};
    if (xmotion.y < band_low || xmotion.y > band_high) return false;
    if (xmotion.x < graph.bounds[0][0] ||
        xmotion.x > graph.bounds[0][1]) return false;
    auto found = upper_bound(band_info.begin(), band_info.end(), xmotion.x,
                             [](const unsigned int lhs,
                                const CytobandInfo & rhs) {
                               return lhs < rhs.abspos_start;
                             });
    if (found == band_info.begin()) return false;
    --found;
    if (xmotion.x > static_cast<int>(found->abspos_stop)) return false;
    graph.status = found->name;
    graph.draw_status(true);
    graph.status = "";
    return true;
  }
  if (event.type != Event::PreDraw) return false;

  // Draw stains
  band_info.clear();
  for (const CytobandInfo & band : cytobands) {
    if (band.abspos_start > graph.range[0][1] ||
        band.abspos_stop < graph.range[0][0]) continue;
    const double start{max(1.0 * band.abspos_start, graph.range[0][0])};
    const double stop{min(1.0 * band.abspos_stop, graph.range[0][1])};
    const double avail_height{(graph.bounds[1][0] - graph.border_width) * 0.6};
    const int x_start{graph.coord(0, start)};
    const int x_stop{graph.coord(0, stop)};
    band_low = graph.bounds[1][1] - 2.5 * avail_height;
    band_high = band_low + avail_height;
    XFillRectangle(graph.display(), graph.pixmap, band.gc, x_start, band_low,
                   x_stop - x_start, avail_height);
    band_info.emplace_back(x_start, x_stop, band.name);
  }

  return false;
}

int main(int argc, char* argv[]) try {
  // Process optional command line arguments
  --argc;

  bool set_geometry{false};
  unsigned int width{1200};
  unsigned int height{1000};
  int x_off{0};
  int y_off{0};
  char ** initial{nullptr};
  bool fullscreen{false};
  while (argc) {
    if (argv[1][0] == '-') {
      const string option{argv[1]};
      if (option == "--geometry") {
        const string geometry{argv[2]};
        istringstream geometry_stream{geometry.c_str()};
        char dummy;
        geometry_stream >> width >> dummy >> height >> x_off >> y_off;
        if (!geometry_stream) {
          throw Error("Problem parsing geometry") << geometry;
        }
        set_geometry = true;
        if (false)
          cerr << "Geometry set to width " << width << " height " << height
               << " x offset " << x_off << " y offset " << y_off << endl;
        argv += 2;
        argc -= 2;
      } else if (option == "--initial-view") {
        initial = argv + 2;
        argc -= 5;
        argv += 5;
      } else if (option == "--n-threads") {
        n_threads = atoi(argv[2]);
        argc -= 2;
        argv += 2;
      } else if (option == "--fullscreen") {
        fullscreen = true;
        argc -= 1;
        argv += 1;
      } else {
        throw Error("Unrecognized command line option") << option;
      }
    } else {
      break;
    }
  }

  if (argc < 4)
    throw Error("usage: ggraph plain|genome|cn ref xn,yn1[:l],yn2[:l]... "
                "data_file ...");

  if (fullscreen && set_geometry) {
    cerr << "Note that --fullscreen option overrides --geometry" << endl;
    sleep(1);
  }

  // Process arguments
  const string type{argv[1]};
  const string ref_name{argv[2]};
  const std::string columns{argv[3]};
  argc -= 3;
  argv += 3;

  // Names of input files
  const vector<string> names{[argc, argv] () {
      vector<string> result;
      result.reserve(argc);
      for (int a{0}; a != argc; ++a) {
        result.emplace_back(argv[a + 1]);
      }
      return result;
    }()};

  // Read in columns from data file
  const vector<Columns> results{[&names, &columns] () {
      ThreadPool pool(n_threads);
      vector<future<Columns>> futures;
      for (const string & name : names) {
        futures.push_back(
            pool.run([&columns] (const string & file_name) {
                return Columns{file_name, 1000000,
                      columns, columns.find_first_not_of(
                          "0123456789,") != string::npos};
              }, name));
      }
      vector<Columns> result;
      result.reserve(names.size());
      for (auto & fut : futures) {
        result.push_back(fut.get());
      }
      return result;
    }()};

  const vector<unsigned char> do_lines{[&results] () {
      vector<unsigned char> result;
      // Assumes all loaded files share the same column names, as is expected
      const Columns & cols{results[0]};
      for (unsigned int c{1}; c != cols.n_cols(); ++c) {
        if (cols.type(c).empty() || cols.type(c)[0] == 'p') {
          result.push_back(0);
        } else if (cols.type(c)[0] == 'l') {
          result.push_back(1);
        } else {
          throw Error("Unknown column type") << cols.type(c);
        }
      }
      // Special case for two dependent variables and no type specifications
      // make the second dependent variable lines
      if (cols.n_cols() == 3 && cols.type(1).empty() && cols.type(2).empty()) {
        result[1] = 1;
      }
      return result;
    }()};


  // App to display multiple windows
  X11App app;

  if (fullscreen) {
    width = app.display_size[0];
    height = app.display_size[1];
    x_off = 0;
    y_off = 0;
  }

  // Rearrange data in specific X11Graph format for all individuals
  const int n_sets(static_cast<unsigned int>(results.size()));
  const int n_y(results.front().n_cols() - 1);
  std::vector<std::vector<const std::vector<double> *>> data(
      n_sets * n_y, std::vector<const std::vector<double> *>(2));
  for (int r{0}; r != n_sets; ++r) {
    for (int y{0}; y != n_y; ++y) {
      data[n_sets * y + r][0] = &results[r][0];
      data[n_sets * y + r][1] = &results[r][y + 1];
    }
  }

  // All individuals graph
  X11Graph & graph{X11Graph::create_whole(app, data,
                                          width, height, x_off, y_off,
                                          "G-Graph")};
  X11Graph * graphp{&graph};

  // Adjust series positions, assign names, and make some series lines only
  for (int r{0}; r != n_sets; ++r) {
    const Columns & result{results[r]};
    const double scale{n_sets <= 4 ? 1.0 : 4.5 / n_sets};
    for (int y{0}; y != n_y; ++y) {
      graph.series_radios[n_sets * y + r].specification.y =
          2 + scale * (1.25 * n_y * (n_sets - r - 1) + (n_y - y - 1));
      const string name{remove_substring(remove_substring(
          remove_including_final(names[r], '/'),
          "_results.txt"), ".varbin.data.txt") + "  " + result.name(y + 1)};
      graph.series_radios[n_sets * y + r].description =
          string("Toggle display for ") + name;
      // + " " + std::to_string(y) + " " + std::to_string(r);
      if (do_lines[y]) graph.series_only_lines[n_sets * y + r] = true;
    }
    if (n_sets > 1) {
      Radio testradio{"Place this series on top", graphp, {-1, 2}};
      graph.extra_radios.push_back(Radio{
            string("Place series") + (n_y > 1 ? " group" : "") +
                " on top", graphp,
            {-0.7, 2 + scale *
                  (1.25 * n_y * (n_sets - r - 1) + (n_y - 0.5 - 1))},
            {[graphp, r, n_y, n_sets]() {
                // Find location of series in ordering list
                vector<unsigned int> & order{graphp->series_order};
                vector<unsigned int>::iterator riter{
                  find(order.begin(), order.end(), r)};
                const unsigned int rindex{
                  static_cast<unsigned int>(riter - order.begin())};
                for (int y{0}; y != n_y; ++y) {
                  vector<unsigned int>::iterator togo{
                    order.begin() + (y + 1) * n_sets};
                  order.insert(togo, r + y * n_sets);
                  vector<unsigned int>::iterator toremove{
                    order.begin() + y * n_sets + rindex};
                  order.erase(toremove);
                }
                // graphp->show_order("after");
                graphp->draw();
              },
                  [r, graphp, n_sets]() {
                    return (graphp->series_order[n_sets - 1]) !=
                        static_cast<unsigned int>(r);
                  }}});
      graph.extra_radios.back().radius_scale = 0.5;
      graph.extra_radios.back().id = r;
      graph.radios.push_front(&graph.extra_radios.back());
    }
  }

  if (type == "cn") {
    // Are some Ys ratios?  Then change scale for ratio lines
    const vector<double> cn_lines{[&results]() {
        bool some_ratios{false};
        for (unsigned int c{1}; c != results.front().n_cols(); ++c)
          if (results.front().name(c).find("ratio") != string::npos ||
              results.front().name(c).find("seg.mean") != string::npos)
            some_ratios = true;
        vector<double> result{0.0, 1.0, 2.0, 3.0, 4.0, 5.0};
        if (some_ratios) for (double & line : result) line /= 2;
        return result;
      }()};

    graph.add_call_back("Toggle ratio lines", X11Graph::CallBack{std::bind(
        &add_ratio_lines, cn_lines, _1, _2)});

    graph.grid_radios[1][1].toggled = false;
    graph.grid_radios[0][1].toggled = false;
  }

  future<bool> gene_future;;
  if (type == "genome" || type == "cn") {
    static const Reference ref{ref_name, true};
    static const CN_abspos cn_abspos{ref};

    // Add drawing callback to add chromosomes and ratio lines
    graph.add_call_back("Toggle chromosome boundaries and names",
                        X11Graph::CallBack{std::bind(
                            &add_chromosomes, std::cref(ref),
                            std::cref(cn_abspos), _1, _2)});

    graph.add_call_back(
        "Toggle chrpos and gene display (only shown for X axis range below " +
        std::to_string(max_gene_mb) + "MB for genes and " +
        std::to_string(max_name_mb) + " MB for names)",
        X11Graph::CallBack{std::bind(&add_genes, std::ref(ref),
                                     std::cref(cn_abspos), _1, _2)});

    // Preload gene info to avoid wait time after zoom
    static ThreadPool pool(1);
    gene_future = pool.run(add_genes, std::cref(ref), std::cref(cn_abspos),
                            std::ref(graph), Event());

    // Cytobands
    graph.add_call_back("Toggle cytobands and names",
                        X11Graph::CallBack{std::bind(
                            &add_cytobands, std::cref(ref),
                            std::cref(cn_abspos), _1, _2)}, true, false);

    // Other genome-specific tweaks
    graph.grid_radios[0][0].toggled = false;
    graph.grid_radios[1][0].toggled = false;
    graph.log_radios[0].actions.visible = [] () { return false; };
  }

  graph.n_threads(n_threads);

  // Process initial view command line arguments
  if (initial) {
    graph.get_range();
    graph.set_range(0, atof(initial[0]), atof(initial[1]));
    graph.set_range(1, atof(initial[2]), atof(initial[3]));
    graph.prepare();
  }

  // Run the app
  app.run();

  if (gene_future.valid()) gene_future.get();

  return 0;
} catch (Error & e) {
  cerr << "paa::Error:" << endl;  cerr << e.what() << endl;
  return 1;
} catch (exception & e) {
  cerr << "std::exception" << endl;
  cerr << e.what() << endl;
  return 1;
} catch (...) {
  cerr << "unknown exception was caught" << endl;
  return 1;
}
