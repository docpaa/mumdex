//
// mumdex_poster
//
// Generate a poster describing mumdex
//
// Copyright 2018 Peter Andrews @ CSHL
//

#include <exception>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "error.h"
#include "layout.h"

using std::cout;
using std::cerr;
using std::endl;
using std::exception;
using std::ifstream;
using std::istringstream;
using std::ofstream;
using std::ostringstream;
using std::string;

using paa::Error;
using paa::Geometry;
using paa::Point;

constexpr int inch{72};
constexpr int foot{12};
constexpr int points_per_foot{inch * foot};
// constexpr int page_width{4 * points_per_foot};
// constexpr int page_height{3 * points_per_foot};
constexpr int page_width{48 * inch};
constexpr int page_height{44 * inch};
constexpr int footer_height{inch * 3};
constexpr int header_height{7 * inch / 2};
constexpr int header_pos{page_height - header_height};
constexpr int ncol{4};
constexpr int colw{page_width / ncol};
constexpr int colimw{15 * colw / 16};

string include_eps(const string & file, double width,
                   const bool is_width, const bool is_top) {
  ifstream eps_file{file.c_str()};
  if (!eps_file) throw Error("Could not open eps file") << file;

  string line;
  string eps;
  int w{0};
  int h{0};
  while (getline(eps_file, line)) {
    if (line.find("%%BoundingBox") != string::npos) {
      istringstream line_stream{line.c_str()};
      string text;
      int dummy;
      line_stream >> text >> dummy >> dummy >> w >> h;
    } else {
      const size_t sp_pos{line.find("showpage")};
      if (sp_pos != string::npos) {
        line.replace(sp_pos, 8, "");
      }
    }
    eps += line + "\n";
  }

  ostringstream result;
  result << " save" << endl;
  result << "gsave" << endl;
  const double scale{1.0 * width / (is_width ? w : h)};
  double height = scale * h;
  if (!is_width) {
    height = width;
    width = scale * w;
  }
  result << -0.5 * width << " "
         << (is_top ? -1.0 : -0.5) * height << " rmoveto" << endl;
  result << "currentpoint translate " << endl;
  result << scale << " " << scale << " scale" << endl;
  result << eps << endl;
  result << "grestore" << endl;
  result << "restore" << endl;
  result << "0 " << height << " fontheight add neg rmoveto" << endl;
  return result.str();
}

string include_eps(const string & file, double width, bool is_width = true) {
  return include_eps(file, width, is_width, true);
}

string include_eps(const string & file,
                   const double x_mid, const double y_mid,
                   double width, const bool is_width) {
  ostringstream result;
  result << " " << x_mid << " " << y_mid << " moveto ";
  result << " gsave " << width
         << " 1.1 mul hh 0.65 mul wbox grestore " << endl;

  return result.str() + include_eps(file, width, is_width, false);
}


int main(int argc, char**) try {
  if (--argc != 0) throw Error("usage: mumdex_poster");

  const double imsf{0.97};

  ostringstream poster;
  poster
      << "%!PS-Adobe-3.0\n"
      << "%%BoundingBox: 0 0 " << page_width << " " << page_height << "\n"
      << R"xxx(
%%Pages: 1
%%Title: MUMdex presentation
%%Creator: Peter Andrews
%%CreationDate: Today
%%EndComments
%%BeginProlog
)xxx"
      << "/pgw " << page_width << " def\n"
      << "/pgh " << page_height << " def\n"
      << "/inch " << inch << " def\n"
      << "/ncol " << ncol << " def\n"
      << "/colw " << colw << " def\n"
      << "/colimw " << colimw << " def\n"
      << "/hh " << header_height << " def\n"
      << "/hpos " << header_pos << " def\n"
      << "/fh " << footer_height << " def\n"
      << R"xxx(
/body_height hpos fh sub def
/sf {/Helvetica findfont exch floor scalefont setfont} bind def 
/jcx {dup stringwidth pop 2 div neg} bind def
/jc {jcx 0 rmoveto} bind def
/jcs {gsave jc show grestore} def
/jrx {dup stringwidth pop neg} bind def
/jr {jrx 0 rmoveto} bind def

/cpy { currentpoint exch pop } def

% string fontsize jcsl    assumes pos at top center of desired string
/jcsl {

gsave

% get copy of fontsize and set font
/FS exch def FS sf

% move to display position and show string
dup stringwidth pop 2 div neg 1.25 FS mul neg rmoveto show

grestore

% move to desired end location
0 FS 1.5 mul neg rmoveto

} def

/jcsh {
% gsave currentpoint newpath moveto colw 2 div neg 0 rmoveto colw 0 rlineto stroke grestore
jcsl
halfline
} def

% From the bluebook examples, the Postscript language reference manual
/wordbreak ( ) def
/BreakIntoLines {
    /proc exch def
    /linewidth exch def
    /textstring exch def
    /breakwidth wordbreak stringwidth pop def
    /curwidth 0 def
    /lastwordbreak 0 def
    /startchar 0 def
    /restoftext textstring def
    {
        restoftext wordbreak search
        {
            /nextword exch def pop
            /restoftext exch def
            /wordwidth nextword stringwidth pop def
            curwidth wordwidth add linewidth gt
            {
                textstring startchar
                lastwordbreak startchar sub
                getinterval proc
                /startchar lastwordbreak def
                /curwidth wordwidth breakwidth add def
            } {
                /curwidth curwidth wordwidth add breakwidth add def
            } ifelse
            /lastwordbreak lastwordbreak
            nextword length add 1 add def
        }
        { pop exit }
        ifelse
    } loop
    /lastchar textstring length def
    textstring startchar lastchar startchar sub
    getinterval proc
} def

% /showpoint { gsave currentpoint newpath inch 3 div 0 360 arc stroke grestore} def
/showpoint { } def
/boxwidth colw 63 mul 64 div def
/boxprint {
gsave
boxwidth 2 div neg 0 rmoveto
show
grestore
} def

/titleheight 62 def
/titlefont { titleheight sf } def
/authorheight 40 def
/authorspacing authorheight 1.2 mul def
/authorfont { authorheight sf } def
/affilheight 30 def
/affilfont { affilheight sf } def
/h2height 50 def
/h3height 35 def

/urlheight 35 def
/urlspacing urlheight 1.2 mul def
/urlfont { urlheight sf } def
/fontheight 20 def
/fontspacing fontheight 1.2 mul def
/fontused { fontheight sf } def 
/backline { 0 fontheight rmoveto } def
/nextline { 0 fontheight neg rmoveto } def
/halfline { 0 fontheight 2 div neg rmoveto } def

/h1fh hh def
/h1f {h1fh sf} def
/h2fh {hh 2 div} def
/h2f {h2fh sf} def
/h3fh {hh 3 div} def
/h3f {h3fh sf} def
/h4fh {hh 4 div} def
/h4f {h4fh sf} def
/h5fh {hh 5 div} def
/h5f {h5fh sf} def
/h6fh {hh 6 div} def
/h6f {h6fh sf} def
/h7fh {hh 7 div} def
/h7f {h7fh sf} def
/h8fh {hh 8 div} def
/h8f {h8fh sf} def

/tpos colw 16 div def
/twid colw 7 mul 8 div def

/justify {
  /www exch def 
  /jlh fontheight 1.4 mul neg def

   0 fontheight neg rmoveto
   www { gsave www 2 div neg 0 rmoveto show grestore 0 jlh rmoveto }
   BreakIntoLines
  0 jlh 2 div rmoveto
} def

/cbox {
gsave 
/bxh exch def /bxw exch def
/rad 10 def
/bxhl bxh rad 2 mul sub def
/bxwl bxw rad 2 mul sub def
currentpoint translate
0 0 0 setrgbcolor
bxwl 2 div neg bxh 2 div neg rmoveto bxwl 0 rlineto
bxwl 2 div bxhl 2 div neg rad -90 0 arc 0 bxhl rlineto
bxwl 2 div bxhl 2 div rad 0 90 arc bxwl neg 0 rlineto
bxwl 2 div neg bxhl 2 div rad 90 180 arc 0 bxhl neg rlineto
bxwl 2 div neg bxhl 2 div neg rad 180 270 arc
closepath
gsave
setrgbcolor
fill
grestore
3 setlinewidth
stroke
grestore
} def

/wbox {

/bxh exch def /bxw exch def
/rad 10 def
/bxhl bxh rad 2 mul sub def
/bxwl bxw rad 2 mul sub def
gsave currentpoint translate
0 0 0 setrgbcolor
bxwl 2 div neg bxh 2 div neg rmoveto bxwl 0 rlineto
bxwl 2 div bxhl 2 div neg rad -90 0 arc 0 bxhl rlineto
bxwl 2 div bxhl 2 div rad 0 90 arc bxwl neg 0 rlineto
bxwl 2 div neg bxhl 2 div rad 90 180 arc 0 bxhl neg rlineto
bxwl 2 div neg bxhl 2 div neg rad 180 270 arc
closepath
gsave
1 1 1 setrgbcolor
fill
grestore
3 setlinewidth
stroke
grestore
} def

/gbox {
/bxh exch def /bxw exch def
/rad 10 def
/bxhl bxh rad 2 mul sub def
/bxwl bxw rad 2 mul sub def
gsave currentpoint translate
black
bxwl 2 div neg bxh 2 div neg rmoveto bxwl 0 rlineto
bxwl 2 div bxhl 2 div neg rad -90 0 arc 0 bxhl rlineto
bxwl 2 div bxhl 2 div rad 0 90 arc bxwl neg 0 rlineto
bxwl 2 div neg bxhl 2 div rad 90 180 arc 0 bxhl neg rlineto
bxwl 2 div neg bxhl 2 div neg rad 180 270 arc
closepath
gsave
gray
fill
grestore
3 setlinewidth
stroke
grestore
} def

/wbox2 {
/bxh exch def /bxw exch def
gsave
0 0 0 setrgbcolor
bxw 2 div neg bxh 2 div neg rmoveto
bxw 0 rlineto 0 bxh rlineto bxw neg 0 rlineto closepath
gsave
1 1 1 setrgbcolor
fill
grestore
3 setlinewidth
stroke
grestore
} def

/bubble1 { /bxh exch def /bxw exch def
gsave
black
bxw 2 div neg bxh neg rmoveto
bxw 0 rlineto 0 bxh rlineto bxw neg 0 rlineto closepath
gsave
gray
fill
grestore
stroke
grestore
} def

/bdefs {
/bxh exch def
/bxw exch def
currentpoint
/y exch bxh 2 div sub def
/x exch def

/cr bxw 120 div def
/xlw bxw 2 cr mul sub def
/hxlw xlw 2 div def
} def


% width bottomcap -> drawing, moving down
/bottomcap {
} def

/bubble {
showpoint
gsave
bdefs
/ylw bxh 2 cr mul sub def
/hylw ylw 2 div def

black
x y moveto
bxw 2 div neg bxh 2 div neg rmoveto
showpoint
cr 0 rmoveto
xlw 0 rlineto

x hxlw add y hylw sub cr -90 0 arc
0 ylw rlineto
x hxlw add y hylw add cr 0 90 arc
xlw neg 0 rlineto
x hxlw sub y hylw add cr 90 180 arc
0 ylw neg rlineto
x hxlw sub y hylw sub cr 180 270 arc
closepath

gsave
gray
fill
grestore
stroke
grestore
} def

/super {
gsave
currentpoint pop /xip exch def
dup sf
0 exch 0.8 mul rmoveto
show
currentpoint pop /xfp exch def
grestore
xfp xip sub 0 rmoveto
} def

/borderline {
  gsave
  currentpoint newpath moveto
  [ 1 3 2 2 1 ] 0 setdash 
  2 setlinewidth
  % 0.7 0.7 0.7 setrgbcolor
  /rad 20 def
  colw 2 div neg fontheight 2 div rad add rmoveto
  currentpoint exch rad add exch rad 180 270 arc
  colw rad 2 mul sub 0 rlineto
  currentpoint rad add rad 270 360 arc
  stroke
  % colw 2 div neg fontheight 2 div rad add neg rmoveto
  % [ ] 0 setdash
  grestore
} def

%%ENDProlog
%%Page: 1 1
%%;
)xxx";

  poster << R"xxx(
/grayv 0.95 def
/gray { grayv grayv grayv setrgbcolor} def

/blackv 0 def
/black { blackv blackv blackv setrgbcolor} def

/whitev 1 def
/white { whitev whitev whitev setrgbcolor} def

/red { 1 0 0 setrgbcolor} def

newpath

false
{
gsave
gray
0 hpos moveto pgw 0 rlineto 0 hh rlineto pgw neg 0 rlineto closepath
fill
grestore
} if

/extra -40 def

% Header line
0 0 0 setrgbcolor inch 16 div setlinewidth
% 0 hpos moveto pgw 0 rlineto

hh hpos hh 2 div add moveto pgw hh 2 mul sub 0 rlineto stroke newpath

pgw 2 div hpos hh 0.9 mul add moveto
pgw 0.76 mul hh 0.9 mul bubble

% Column separator lines
1 1 ncol 1 sub {colw mul fh extra add moveto 0 hpos fh sub extra sub rlineto} for
pgw 2 div 0 moveto 0 fh rlineto
stroke

newpath white pgw 2 div 0 moveto 0 110 rlineto currentpoint stroke black
gsave newpath moveto currentpoint -100 0 rmoveto 200 0 rlineto stroke grestore  
newpath 8 0 360 arc gsave white fill grestore stroke

newpath

% main text
titlefont pgw 2 div hpos hh 8 div 5 mul add moveto

(MUMdex: MUM-based structural variation detection and the G-Graph Copy Number GUI) jcs

% authors
/supera { exch show affilheight 0.75 mul super } def

pgw 2 div hpos hh 3 div add moveto
/authors {
authorfont
(Peter A. Andrews) (1) supera (, ) show
(Ivan Iossifov) (1,2) supera (, ) show
(Jude Kendall) (1)  supera (, ) show
(Steven Marks) (1,2) supera (, ) show
(Lakshmi Muthuswamy) (2) supera (, ) show
(Zihua Wang) (1) supera (, ) show
(Michael Wigler) (1) supera (, ) show
(and Dan Levy) (1) supera
} def
gray

gsave
pgw pgh moveto 
currentpoint pop /xia exch def
authors
currentpoint pop /xfa exch def
grestore

xfa xia sub 2 div neg 0 rmoveto
black
authors

% affiliations
affilfont pgw 2 div hpos hh 9 div add moveto
(1: Simons Center for Quantitative Biology, Cold Spring Harbor Laboratory, Cold Spring Harbor, NY 11724, USA; and 2: New York Genome Center, New York, NY 10013, USA) jcs

% column 1
showpoint
pgw 1 mul ncol 2 mul div hpos moveto nextline
(MUMs: Maximal Unique Matches) h2height jcsh
(MUMs are matches between the genome reference and sequencing read that cannot be extended and that match the reference genome in exactly one place. We use MUMs to help us find mutation and interpret the genome. ) twid fontused justify

(A toy example demonstrating MUMs) h3height jcsh
)xxx" << include_eps("mums.eps", colimw) << R"xxx(
(The short reference and read displayed above have several Maximum Unique Matches between them, each displayed as a differently colored line. The MUMs tell us where each portion of the read belongs in the reference. Structural variant detection is based on finding novel adjacencies in the sample reads relative to the reference and across populations. ) twid fontused justify

(Bridges connect two MUMs in a read) h3height jcsh
)xxx" << include_eps("bridges.eps", colimw) << R"xxx(
(Two more bases are spanned in the read than in the reference. Each MUM in the read induces its own reference coordinate system on the read. The two coordinate systems for this read are shifted by 2.  This event has an invariant of 2, indicating a possible two base insertion event. ) twid fontused justify

(Bridges bracket mutational events) h3height jcsh
)xxx" << include_eps("paper_1.eps", colimw) << R"xxx(
(A bridge is any pair of MUMs in a read, and is necessarily caused by a difference between the read and the reference genome. We count identical or similar bridges between MUMs in a read and compare counts among families and across populations.  Some bridges are easily interpretable as to what caused the bridge signature - i.e. a SNP, indel, inversion, translocation. ) twid fontused justify


(A geometric view of common event types) h3height jcsh
)xxx" << include_eps("wedges.eps", colimw * 8.5 / 11) << R"xxx(
(This is an alternate geometrical way of viewing bridge event types that may help to better show the actual event. Two toy reference chromosomes are displayed at top, with the height of the reference bases continually increasing along each chromosome as an aid to visualization when event grafting reorders pieces of the reference into a sample read. ) twid fontused justify

% column 2
pgw 3 mul ncol 2 mul div hpos moveto
nextline
(The MUMdex software package) h2height jcsh

currentrgbcolor 0 0 0.8 setrgbcolor 
(http://mumdex.com/) urlheight jcsl
(https://github.com/docpaa/mumdex/) urlheight jcsh
setrgbcolor
(MUMdex can be downloaded from either the MUMdex or github websites.  It is written entirely in C++ and requires a C++11 or later compiler such as GCC or clang. Optional MUMdex components require the gsl and Xlib libraries to be installed. MUMdex software compiles and runs under Linux, Mac OSX or Windows with Cygwin installed. ) twid fontused justify

(MUMdex uses a suffix array to align reads) h3height jcsh
 gsave colw 4 div neg 0 rmoveto
)xxx"    << include_eps("suffix.eps", colimw * 0.45) << R"xxx(
colw 4 div 0 rmoveto  currentpoint grestore colw 4 div 0 rmoveto gsave
nextline nextline
(A suffix array is a computer science data structure that can be used to quickly find all MUMs between a read and the reference genome.  The MUMdex default suffix array requires 120 GB of space, but a slower 32 GB version is also available for memory-strapped machines. The suffix array and associated data structures are built during a one-time initial procedure and then loaded using memory-mapped files for subsequent use. ) colimw 0.45 mul fontused justify
(The suffix array example at left encodes the string AGTTAGTCC as the numeric array 9408751362, which is a list of indexes into the original string sorted by string suffix that groups similar suffixes together, enabling fast lookup for querying read locations. ) colimw 0.45 mul fontused justify
grestore
moveto

nextline
(The MUMdex binary alignment format) h3height jcsh
)xxx"    << include_eps("mumdex.eps", colimw) << R"xxx(
(The MUMdex alignment format compactly encodes MUM alignment and read sequence information in a lossless binary format using the reference encoding method. Sequence that is part of a MUM can be looked up from the reference genome, so only bases not covered by MUMs need to be stored for each read.  All MUMs for each read pair are stored together, and individual MUMs or read pairs can be looked up using a genome order index. MUMdex data structures are typically loaded as needed using memory mapped files. ) twid fontused justify

(The MUMdex genome analysis pipeline) h3height jcsh
)xxx"    << include_eps("software.eps", colimw) << R"xxx(
(The MUMdex pipeline first aligns read pairs and stores the information in the MUMdex alignment format. Next, it extracts all bridges from each sample's MUMdex file independently.  Then, it analyzes bridges for all samples in a population simultaneously over a window to detect events of interest.  Family information is used to find de novo structural variants and the population is used to filter for rare events. ) twid fontused justify


(MUMdex time and space characteristics) h3height jcsh
)xxx"    << include_eps("mumdex_time_space.eps", colimw) << R"xxx(
(MUMdex sample-level processing, which includes read pair alignment, MUMdex alignment format creation and the bridge extraction step, can be performed in less than 4 hours and 90 GB of space per sample. ) twid fontused justify

% column 3
pgw 5 mul ncol 2 mul div hpos moveto
nextline
(Analysis and validation results) h2height jcsh
(We processed whole genome data from 510 quad-format Simon's Simplex Collection families, each with one autistic proband, and analyzed the results to find 6278 de novo structural variation candidates in probands and / or siblings.  4740 of these events were private to one family, while the rest were seen in up to 5 families total. MUMdex has the ability to find events of any size without significant event-size bias. MUMdex can also detect substitution (SNP) events, but those were ignored for this study. ) twid fontused justify

(6278 structural variant candidates detected) h3height jcsh
)xxx"    << include_eps("table_types.eps", colimw) << R"xxx(
(The table above lists counts for events found by event size, whether they were insertions or deletions, and whether they were found in one or more than one family.  Not tabulated above are the 13 translocation and 7 inversion events that were also detected by the method. ) twid fontused justify

(Event size and genomic position distribution) h3height jcsh
gsave colimw 4 div neg 0 rmoveto 
)xxx"    << include_eps("supplementary_figure_4.eps", imsf * colimw / 2)
         << R"xxx(
20 0 rmoveto
(Most de novo events found are indels of size 30 bases and smaller. Deletions outnumber insertions by about a factor of 5.  The events are distributed evenly across the genome, except for a marked reduction in density on the X and Y chromosomes. ) colimw 2 div 0.9 mul fontused justify
-20 0 rmoveto

grestore colimw 4 div 0 rmoveto
)xxx"    << include_eps("karyogram.eps", imsf * colimw / 2) << R"xxx(
 colimw 4 div neg 0 rmoveto

(A typical 24 base microsatellite contraction) h3height jcsh
)xxx"    << include_eps("supplementary_figure_8.eps", colimw) << R"xxx(
(This is an example of a typical 24 base deletion candidate, which appears to be a microsatellite contraction event.  Note how the MUMs in the right read overlap in the read (middle) but not in the genome (top and bottom). ) twid fontused justify

(Validation results summary) h3height jcsh
(We performed molecular validation on 106 randomly selected single-family de novo candidates from two size groups and 101 of them validated. This corresponds to a 95% overall validation rate. The events validated included 25 / 26 of the larger events tested (100 bp and more) and 76 / 80 of the shorter events. ) twid fontused justify

(The MUMdex bioRxiv preprint paper) h3height jcsh
currentrgbcolor 0 0 0.8 setrgbcolor 
(https://www.biorxiv.org/content/early/2016/09/30/078261) urlheight 0.9 mul jcsh
setrgbcolor

(The MUMdex software package and the analysis results contained in this poster are described in more detail in the MUMdex bioRxiv paper. ) twid fontused justify


% column 4
pgw 7 mul ncol 2 mul div hpos moveto
nextline
(Copy Number, G-Graph and SMASH) h2height jcsh
(The MUMdex genome analysis package includes software to perform MUM-based copy number analysis.  This includes code to select appropriate bin boundaries for your dataset and to generate profiles for individual samples at any desired resolution. The G-Graph visualization GUI allows researchers and clinicians to efficiently perform multiscale interactive data exploration, even for millions of data points. ) twid fontused justify

(SSC quad WG MUM-based 500,000 bin CN) h3height jcsh
)xxx"    << include_eps("ggraph.eps", colimw) << R"xxx(
(We checked MUMdex SSC large (over 1000 bases) insertion and deletion de novo candidates for evidence of each split-read detected event in the copy number signal. We were able to confirm 79% of deletions and 71% of duplications in this fashion. Above is a screenshot of the G-Graph GUI interface for a transmitted event, showing gene annotation from the UCSC browser. ) twid fontused justify

% A block
(G-Graph tiled, help and other views) h3height jcsh
 showpoint colimw 6 div neg 0 rmoveto
 /yip currentpoint exch pop def
)xxx"    << include_eps("tiled.eps", imsf * colimw * 2 / 3)
         << include_eps("help.eps", imsf * colimw * 2 / 3)
         << " colimw 6 div 0 rmoveto /yfp currentpoint exch pop def "
         << " colimw 3 div yfp yip sub neg rmoveto  "
         << include_eps("x11plot.eps", imsf * colimw / 3)
         << include_eps("colors.eps", imsf * colimw / 3) << R"xxx(
colimw 3 div neg 0 rmoveto
(Above we show four views from the G-Graph and related visualization programs.  These include the tiled view mode, a help screen view, a series chooser interface and a series color selector window view. ) twid fontused justify

(Lower cost SMASH-based CN method) h3height jcsh
)xxx"    << include_eps("smash.eps", colimw) << R"xxx(
(The SMASH sequencing method is a 5 times cheaper way to get data for copy number analysis. SMASH packs on average 5 unique genome reference maps per read pair so you get more information for less sequencing money.  The SMASH method uses the MUMdex aligner to decompose each read pair in-silico into its constituent mappings. We then use MUMdex copy number analysis software to generate segmented copy number profiles. ) twid fontused justify

% signature and copyright
/endh fh 1.5 mul def
colw 1.5 mul endh moveto
nextline
(MUMdex software is open source) h3height jcsh
backline
halfline
(copyright 2018 by Peter Andrews @ CSHL) h3height jcsh


/footer
{

% gsave 0 0 moveto pgw 0 rlineto 0 fh rlineto pgw neg 0 rlineto closepath gray fill grestore

% acknowledgements
/fwid colw 1.9 mul def
/fhei fh 0.9 mul def
/bwid fwid 0.95 mul def

colw fh extra add moveto
fwid fhei extra add bubble

(Acknowledgements) h2height jcsh

(We thank the Mike Wigler laboratory at Cold Spring Harbor Laboratory for supporting this work, the Simons Foundation for project funding, and the New York Genome Center for performing the sequencing and providing the datasets used.) bwid fontused justify

% references
/extra 40 def

colw 3 mul fh extra add moveto
fwid fhei extra add bubble

(Selected References) h2height jcsh

(1. Andrews PA, Iossifov I, Kendall J, Marks S, Muthuswamy L, Wang Z, Levy D, Wigler M. (2016) MUMdex: MUM-based structural variation detection. bioRxiv 078261; doi: http://dx.doi.org/10.1101/078261 ) bwid fontused justify
backline

(2. Fischbach GD, Lord C. (2010) The Simons Simplex Collection: a resource for identification of autism genetic risk factors. Neuron, 68 (2), 192-195) bwid fontused justify
backline

(3. Kent WJ, Sugnet CW, Furey TS, Roskin KM, Pringle TH, Zahler AM, Haussler D. (2002) The human genome browser at UCSC. Genome Res. 12 (6), 996-1006.) bwid fontused justify
} def

footer

false {
black
newpath moveto [ 80 5 65 5 ] 0 setdash 0 -155 rlineto currentpoint stroke
newpath moveto [ 50 5 35 10 ] 0 setdash 0 -100 rlineto currentpoint stroke
newpath moveto [ 20 15 5 20 ] 0 setdash 0 -60 rlineto stroke
    } if

)xxx";

  // CSHL and NYGC logos
  const int image_width{header_height};
  poster << include_eps("cshl.eps",
                        image_width * 13 / 16,
                        page_height - 0.55 * header_height, image_width, true);
  poster << include_eps("nygc.eps",
                        page_width - image_width * 12 / 16,
                        page_height - 0.55 * header_height, 0.8 * image_width,
                        true);

  poster << "stroke\n";
  poster << "showpage\n"
         << "%%EndPage: 1\n"
         << "%%EOF\n";

  std::cout << poster.str() << std::endl;
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
