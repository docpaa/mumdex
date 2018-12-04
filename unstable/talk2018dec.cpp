//
// talk2018dec.cpp
//
// Dec 6, 2018 at NYGC
//
// Copyright 2018 Peter Andrews
//
// ~/mumdex/talk2018dec
// cp ~/mumdex/talk2018dec ~/web/talk2-18dec/index.cgi
// firefox http://localhost/talk2018dec/
//
#include <climits>
#include <deque>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "error.h"
#include "files.h"
#include "strings.h"
#include "utility.h"

std::string message{""};

namespace paa {

using StringPair = std::pair<std::string, std::string>;
using StringSet = std::set<std::string>;
using WebRoot = StringPair;

class Title {
 public:
  Title(const std::string & title__,
        const std::string & subtitle__,
        const std::string & short_title__) :
      title_{title__},
    subtitle_{subtitle__},
    short_title_{short_title__} { }
  Title(const char * const title__) :  // NOLINT
      title_{get_part(title__, 0)},
    subtitle_{get_part(title__, 1)},
    short_title_{get_part(title__, 2)} { }
  Title(const std::string & title__) :  // NOLINT
      title_{get_part(title__, 0)},
    subtitle_{get_part(title__, 1)},
    short_title_{get_part(title__, 2)} { }

  std::string get_part(const std::string full, unsigned int part) {
    uint64_t pos{0};
    while (part--) {
      pos = full.find(";", pos + 1) + 1;
      if (pos == std::string::npos) return "";
    }
    std::string result;
    std::istringstream in{full.substr(pos)};
    getline(in, result, ';');
    return result;
  }

  std::string title() const { return title_; }
  Title & title(const std::string & title__) {
    title_ = title__;
    return *this;
  }
  std::string subtitle() const { return subtitle_; }
  Title & subtitle(const std::string & subtitle__) {
    subtitle_ = subtitle__;
    return *this;
  }
  std::string short_title() const { return short_title_; }
  Title & short_title(const std::string & short_title__) {
    short_title_ = short_title__;
    return *this;
  }

 private:
  std::string title_;
  std::string subtitle_;
  std::string short_title_;
};

class Employee {
 public:
  Employee(const std::string & first__,
         const std::string & last__,
         const std::string & lab__,
         const std::string & employer__,
         const std::string & title__) :
      first_{first__},
    last_{last__},
    lab_{lab__},
    employer_{employer__},
    title_{title__} { }

  std::string name() const { return first_ + " " + last_; }
  std::string first() const { return first_; }
  Employee & first(const std::string & first__) {
    first_ = first__;
    return *this;
  }
  std::string last() const { return last_; }
  Employee & last(const std::string & last__) {
    last_ = last__;
    return *this;
  }
  std::string lab() const { return lab_; }
  Employee & lab(const std::string & lab__) {
    lab_ = lab__;
    return *this;
  }
  std::string employer() const { return employer_; }
  Employee & employer(const std::string & employer__) {
    employer_ = employer__;
    return *this;
  }
  std::string title() const { return title_; }
  Employee & title(const std::string & title__) {
    title_ = title__;
    return *this;
  }

 private:
  std::string first_;
  std::string last_;
  std::string lab_;
  std::string employer_;
  std::string title_;
};
using Employees = std::vector<Employee>;

class Seminar : public Title {
 public:
  Seminar(const std::string & title__,
          const std::string & subtitle__,
          const std::string & short_title__,
          const std::string & series__,
          const std::string & department__,
          const std::string & room__,
          const std::string & building__,
          const std::string & date__,
          const std::string & year__) :
      Title{title__, subtitle__, short_title__},
    series_{series__},
    department_{department__},
    room_{room__},
    building_{building__},
    date_{date__},
    year_{year__} { }

  std::string series() const { return series_; }
  Seminar & series(const std::string & series__) {
    series_ = series__;
    return *this;
  }
  std::string department() const { return department_; }
  Seminar & department(const std::string & department__) {
    department_ = department__;
    return *this;
  }
  std::string room() const { return room_; }
  Seminar & room(const std::string & room__) {
    room_ = room__;
    return *this;
  }
  std::string building() const { return building_; }
  Seminar & building(const std::string & building__) {
    building_ = building__;
    return *this;
  }
  std::string date() const { return date_; }
  Seminar & date(const std::string & date__) {
    date_ = date__;
    return *this;
  }
  std::string year() const { return year_; }
  Seminar & year(const std::string & year__) {
    year_ = year__;
    return *this;
  }

 private:
  std::string series_;
  std::string department_;
  std::string room_;
  std::string building_;
  std::string date_;
  std::string year_;
};

//
// input
//

using Strings = std::vector<std::string>;
Strings strings(const int n, char ** p) {
  Strings strings_;
  for (int s{0}; s != n; ++s) strings_.push_back(*(p+s));
  return strings_;
}

class InputMethod {
 public:
  InputMethod(std::istream & in__,
              int argc__, char ** argv__,
              char ** envp__,
              const std::string dir__) :
      in_{in__},
    argc_{argc__},
    argv_{argv__},
    envc_{0},
    envp_{envp__},
    dir_{dir__} {
      if (!in_) throw Error("Expected stdin to be good");
      if (argc_ != 1) throw Error("does not expect command line args") << argc_;
      // if (!readable(dir_)) throw Error("cannot read from dir") << dir_;
      while (*(envp_ + envc_) && envc_ < 10000) { ++envc_; }
      if (envc_ == 10000) throw Error("Too much environment");
      for (int e{0}; e != n_envs(); ++e) {
        std::istringstream env_{env(e)};
        std::string name;
        std::string value;
        getline(env_, name, '=');
        getline(env_, value, static_cast<char>(0));
        env_map[name] = value;
      }
  }
  InputMethod(const InputMethod & rhs) = default;
  // in_{rhs.in_}, argc_{rhs.argc_}, argv_{rhs.argv_},
  // envc_{rhs.envc_}, envp_{rhs.envp_}, dir_{rhs.dir_} {}
  InputMethod & operator=(const InputMethod &) = default;

  std::istream & in() const { return in_; }
  int n_args() const { return argc_; }
  std::string arg(const int arg_) const { return *(argv_ + arg_); }
  char ** args_begin() const { return argv_; }
  char ** args_end() const { return argv_ + argc_; }
  Strings args() const { return strings(argc_, argv_); }
  int n_envs() const { return envc_; }
  std::string env(const int env_) const { return *(envp_ + env_); }
  char ** envs_begin() const { return envp_; }
  char ** envs_end() const { return envp_ + envc_; }
  Strings envs() const { return strings(envc_, envp_); }
  std::string env(const std::string & name) const try {
    return env_map.at(name);
  } catch(...) {
    if (name == "REQUEST_URI") return "/talk2018dec/";
    if (name == "SERVER_NAME") return "wigclust19.cshl.edu";
    if (name == "QUERY_STRING") return "page=0";
    if (name == "SCRIPT_NAME") return "/talk2018dec/talk.html";
    throw Error("Problem looking up env in InputMethod") << name;
    return "";
  }
  std::string dir() const { return dir_; }

  std::string env() const {
    std::ostringstream env_;
    for (int e{0}; e != n_envs(); ++e)
      env_ << replace_all(env(e), "<", "&lt;") << "\n";
    return env_.str();
  }

  std::string server_name() const { return env("SERVER_NAME"); }
  std::string request_uri() const { return env("REQUEST_URI"); }
  std::string query_string() const { return env("QUERY_STRING"); }
  std::string script_name() const { return env("SCRIPT_NAME"); }

 private:
  std::istream & in_;
  int argc_;
  char ** argv_;
  int envc_;
  char ** envp_;
  std::map<std::string, std::string> env_map{};
  std::string dir_;
};

std::string fake_message{""};
using BoolString = std::pair<bool, std::string>;
class Allowed {
 public:
  using DefaultSet = std::pair<std::string, StringSet>;
  using AllowedLookup = std::map<std::string, DefaultSet>;
  explicit Allowed(const AllowedLookup & allowed__) :
      allowed_{allowed__} { }

  BoolString get_value(const std::string & name,
                       const std::string & value) const {
    if (1) {  // VERY WEIRD if 0, page lookup fails!
      fake_message += "(look " + name + ", " + value + " in";
      for (const std::string & choice : allowed_.at(name).second) {
        fake_message += " " + choice;
      }
      fake_message += ")\n";
    }
    auto found_name = allowed().find(name);
    if (found_name == allowed().end()) return BoolString{false, ""};
    const StringSet & values{found_name->second.second};
    auto found_value = values.find(value);
    if (found_value == values.end()) {
      // message += "end";
      return BoolString{true, found_name->second.first};  // def
    }
    return BoolString{value == found_name->second.first, value};
  }
  Strings names() const {
    Strings names_;
    for (const auto & name_def : allowed_) names_.push_back(name_def.first);
    return names_;
  }

 private:
  AllowedLookup allowed() const { return allowed_; }
  Allowed & allowed(const AllowedLookup & allowed__) {
    allowed_ = allowed__;
    return *this;
  }
  AllowedLookup allowed_;
};

class QueryVals {
 public:
  explicit QueryVals(const InputMethod & input, const Allowed & allowed) {
    std::map<std::string, std::string> user_values;
    std::istringstream query_string_{input.query_string().c_str()};
    std::string name_vals;
    while (getline(query_string_, name_vals, '&')) {
      std::istringstream name_vals_{name_vals.c_str()};
      std::string name;
      std::string value_;
      getline(name_vals_, name, '=');
      getline(name_vals_, value_);
      user_values[name] = value_;
    }
    for (const auto name : allowed.names()) {
      const std::string user_value{user_values[name]};
      const BoolString def_value_{allowed.get_value(name, user_value)};
      if (0) {
        message += name + ", " + user_value + ", " +
            (def_value_.first ? "T" : "F") + ", " + def_value_.second + "\n";
      }
      lookup[name] = names_.size();
      names_.push_back(name);
      values_.push_back(def_value_.second);
      defaults_.push_back(def_value_.first);
    }
    if (0) {
      message += "page is " + user_values["page"] + ", " + value("page") + "\n";
      for (unsigned int n{0}; n != allowed.names().size(); ++n) {
        message += allowed.names()[n] +
            " " + names_[n] +
            " " + values_[n] +
            " " + (defaults_[n] ? "true" : "false") + "\n";
      }
    }
    output_suffix_ = input.script_name().substr(
        input.script_name().find_last_of('.') + 1);
    server_name_ = input.env("SERVER_NAME");
    request_uri_ = "http://" + input.env("SERVER_NAME") +
        input.env("REQUEST_URI");
    script_name_ = replace_substring(input.env("SCRIPT_NAME"), "index.cgi", "");
  }

  std::string value(const std::string & name) const try {
    return values_[lookup.at(name)];
  } catch(...) {
    throw Error("Problem looking up value in QueryVals::value") << name;
    return "";
  }
  uint64_t ivalue(const std::string & name) const try {
    std::istringstream value_{values_[lookup.at(name)].c_str()};
    uint64_t ivalue_{0};
    value_ >> ivalue_;
    return ivalue_;
  } catch(...) {
    throw Error("Problem looking up value in QueryVals::ivalue") << name;
    return 0;
  }
  BoolString def_value(const std::string & name) const try {
    const uint64_t index{lookup.at(name)};
    return BoolString{defaults_[index], values_[index]};
  } catch(...) {
    throw Error("Problem looking up value in QueryVals::def_value") << name;
    return {false, ""};
  }
  std::string value(const uint64_t index) const {
    return values_[index];
  }
  uint64_t size() const { return names_.size(); }
  const Strings & names() const { return names_; }
  const Strings & values() const { return values_; }
  template <class Value>
  std::string query_string(const std::string name, const Value & value_) const {
    return query_string(name, std::to_string(value_));
  }
  std::string query_string(const std::string name,
                           const char * const value_) const {
    return query_string(name, std::string(value_));
  }
  std::string query_string(const std::string name,
                           const std::string & value_) const {
    Strings name_value_pairs;
    for (unsigned int p{0}; p != names_.size(); ++p) {
      const std::string & pname{names_[p]};
      if (name != pname && defaults_[p]) continue;
      const std::string & pvalue{name == pname ? value_ : values_[p]};
      name_value_pairs.emplace_back(pname + "=" + pvalue);
    }
    if (name_value_pairs.empty()) return "";
    std::ostringstream query;
    query << script_name() + "?";
    for (uint64_t nv{0}; nv != name_value_pairs.size() ; ++nv) {
      const std::string & name_value{name_value_pairs[nv]};
      if (nv) query << "&";
      query << name_value;
    }
    return query.str();
  }

  std::string output_suffix() const { return output_suffix_; }
  std::string request_uri() const { return request_uri_; }
  std::string server_name() const { return server_name_; }
  std::string script_name() const { return script_name_; }

 private:
  std::map<std::string, uint64_t> lookup{};
  Strings names_{};
  Strings values_{};
  std::vector<bool> defaults_{};
  std::string output_suffix_{"html"};
  std::string request_uri_{""};
  std::string server_name_{""};
  std::string script_name_{""};
};

//
// content building blocks
//

class PageElement {
 public:
  using PEP = std::shared_ptr<PageElement>;
  using PEPS = std::vector<PEP>;
  template <class ... Input>
  explicit PageElement(Input && ... input) {
    process(std::forward<Input>(input)...);
  }
  virtual ~PageElement() = default;
  virtual std::string text() const {
    std::ostringstream text_;
    for (const PEP & elem : elements) text_ << elem->text() << "\n";
    return text_.str();
  }
  virtual std::string html() const {
    std::ostringstream html_;
    for (const PEP & elem : elements) html_ << elem->html() << "\n";
    return html_.str();
  }

  template <class ... Input>
  void process(const std::string & string, Input && ... input);

  template <class ... Input>
  void process(const char * const string, Input && ... input);

  template <class Element, class ... Input>
  void process(const Element & element, Input && ... input) {
    elements.push_back(std::make_shared<Element>(element));
    process(std::forward<Input>(input)...);
  }

  void process() { }

  std::string type_html() const {
    return (type_.size() ? " class=\"" + type_ + "\"" : "");
  }
  std::string type() const { return type_; }
  PageElement & type(const std::string & type__) {
    type_ = type__;
    return *this;
  }

 protected:
  std::string type_{""};
  PEPS elements{};
};

template <class Element>
Element & type(const Element & element, const std::string & type_) {
  element.type(type_);
  return element;
}
template <class Element>
Element && type(Element && element, const std::string & type_) {
  element.type(type_);
  return std::move(element);
}

class Text : public PageElement {
 public:
  explicit Text(const std::string & text__) : text_{text__} { }
  virtual ~Text() = default;
  virtual std::string text() const { return text_; }
  virtual std::string html() const { return text_; }

 private:
  std::string text_;
};

template <class ... Input>
void PageElement::process(const std::string & string, Input && ... input) {
  elements.push_back(std::make_shared<Text>(string));
  process(std::forward<Input>(input)...);
}
template <class ... Input>
void PageElement::process(const char * const string, Input && ... input) {
  elements.push_back(std::make_shared<Text>(string));
  process(std::forward<Input>(input)...);
}

template<std::string (*TAG)()>
class Tag : public PageElement {
 public:
  template <class ... Input>
  explicit Tag(Input && ... input) :
      PageElement{std::forward<Input>(input)...} { }

  virtual ~Tag() = default;
  virtual std::string html() const {
    return "<" + TAG() + type_html() + ">" +
        PageElement::html() + "</" + TAG() + ">";
  }
};

std::string bold_fun() { return "b"; }
using bold = Tag<&bold_fun>;
std::string heading_fun() { return "h1"; }
using heading = Tag<&heading_fun>;
std::string heading2_fun() { return "h2"; }
using heading2 = Tag<&heading2_fun>;
std::string heading3_fun() { return "h3"; }
using heading3 = Tag<&heading3_fun>;
std::string span_fun() { return "span"; }
using span = Tag<&span_fun>;
std::string pre_fun() { return "pre"; }
using pre = Tag<&pre_fun>;
std::string div_fun() { return "div"; }
using division = Tag<&div_fun>;
std::string para_fun() { return "p"; }
using para = Tag<&para_fun>;
std::string italic_fun() { return "i"; }
using italic = Tag<&italic_fun>;

class png : public PageElement {
 public:
  png(const std::string & file__,
      const std::string & alt__ = "",
      const std::string & width__ = "",
      const std::string & link__ = "") :    // NOLINT
      file_{file__},
    alt_{alt__.size() ? alt__ : file__},
    width_{width__}, link_{link__} { }
  png(const char * const file__) :    // NOLINT
      file_{file__}, alt_{file__}, width_{""}, link_{""} { }
  virtual ~png() = default;

  virtual std::string text() const { return file(); }
  virtual std::string html() const {
    std::ostringstream html__;
    if (link_.size()) html__ << "<a href=\"" << link_ << "\">";
    html__ <<  "<img src=\"" << file() << "\""
           << " alt=\"" << alt() << "\""
           << " title=\"" << alt() << "\"";
    if (width().size()) html__ << " style=\"width:" << width() << ";\"";
    html__ << " />";
    if (link_.size()) html__ << "</a>";
    return html__.str();
  }

  std::string file() const { return file_; }
  std::string alt() const { return alt_; }
  std::string width() const { return width_; }
  png & width(const std::string & width__) {
    width_ = width__;
    return *this;
  }
  std::string link() const { return link_; }
  png & link(const std::string & link__) {
    link_ = link__;
    return *this;
  }

 private:
  std::string file_;
  std::string alt_;
  std::string width_;
  std::string link_;
};

class PDF {};
class PS {};
class Movie {};

class list : public PageElement {
 public:
  template <class ... Input>
  explicit list(Input && ... input) :
      PageElement{std::forward<Input>(input)...} { }
  virtual ~list() = default;

  virtual std::string html() const {
    std::ostringstream html__;
    html__ << "<ul" << type_html() << ">\n";
    for (const auto & element : elements)
      html__ << "<li>" << element->html() << "</li>\n";
    html__ << "</ul>\n";
    return html__.str();
  }
};


template<std::string (*TYPE)()>
class Type : public PageElement {
 public:
  template <class ... Input>
  explicit Type(Input && ... input) :
      PageElement{std::forward<Input>(input)...} {
    type_ = TYPE();
  }
  virtual ~Type() = default;

  virtual std::string html() const {
    std::ostringstream html__;
    html__ << "<span" << type_html() << ">\n";
    for (const auto & element : elements)
      html__ << element->html() << "\n";
    html__ << "</span>\n";
    return html__.str();
  }
};

std::string bigger_fun() { return "bigger"; }
using bigger = Type<&bigger_fun>;
std::string huge_fun() { return "huge"; }
using huge = Type<&huge_fun>;
std::string huger_fun() { return "huger"; }
using huger = Type<&huger_fun>;
std::string smaller_fun() { return "smaller"; }
using smaller = Type<&smaller_fun>;
std::string tiny_fun() { return "tiny"; }
using tiny = Type<&tiny_fun>;
std::string red_fun() { return "red"; }
using red = Type<&red_fun>;
std::string lower_fun() { return "lower"; }
using lower = Type<&lower_fun>;

class numbered : public PageElement {
 public:
  template <class ... Input>
  explicit numbered(Input && ... input) :
      PageElement{std::forward<Input>(input)...} { }
  virtual ~numbered() = default;

  virtual std::string html() const {
    std::ostringstream html__;
    html__ << "<ol" << type_html() << ">\n";
    for (const auto & element : elements)
      html__ << "<li>" << element->html() << "</li>\n";
    html__ << "</ol>\n";
    return html__.str();
  }
};

class hlink : public PageElement {
 public:
  template <class ... Input>
  hlink(const std::string & href__, Input && ... input) :
      PageElement{std::forward<Input>(input)...},
    href_{href__} {
      if (elements.empty()) elements.push_back(std::make_shared<Text>(href_));
    }

  virtual std::string text() const {
    return PageElement::text() + "(" + href() + ")";
  }
  virtual std::string html() const {
    std::ostringstream html__;
    html__ <<  "<a href=\"" << href() << "\"" + type_html() + ">"
           << PageElement::html() << "</a>";
    return html__.str();
  }

  std::string href() const { return href_; }
  hlink & href(const std::string & href__) {
    href_ = href__;
    return *this;
  }

 private:
  std::string href_;
};

class super : public PageElement {
 public:
  template <class ... Input>
  super(const std::string & base__, const std::string & exponent__,
        Input && ... input) :
      PageElement{std::forward<Input>(input)...},
    base_{base__},
    exponent_{exponent__} { }

  virtual std::string text() const {
    return base() + " ^ " + exponent() + PageElement::text();
  }
  virtual std::string html() const {
    return base() + "<sup>" + exponent() + PageElement::html() + "</sup>";
  }

  std::string base() const { return base_; }
  super & base(const std::string & base__) {
    base_ = base__;
    return *this;
  }
  std::string exponent() const { return exponent_; }
  super & exponent(const std::string & exponent__) {
    exponent_ = exponent__;
    return *this;
  }

 private:
  std::string base_;
  std::string exponent_;
};


class iframe : public PageElement {
 public:
  template <class ... Input>
  iframe(const std::string & href__, Input && ... input) :
      PageElement{std::forward<Input>(input)...},
    href_{href__} {}

  virtual std::string text() const {
    return PageElement::text();
  }
  virtual std::string html() const {
    std::ostringstream html__;
    html__ <<  "<iframe src=\"" << href() << "\"" << type_html() << ">"
           << PageElement::html() << "</iframe>";
    return html__.str();
  }

  std::string href() const { return href_; }
  iframe & href(const std::string & href__) {
    href_ = href__;
    return *this;
  }

 private:
  std::string href_;
};

class Table : public PageElement {
  using Data = const std::vector<Strings>;

 public:
  Table(const std::string & title__, const Data & data__) :
      title_{title__}, data_{data__} {}
  virtual ~Table() = default;

  virtual std::string text() const {
    std::ostringstream text__;
    text__ << title() << "\n";
    for (uint64_t r{0}; r != data_.size(); ++r) {
      for (uint64_t c{0}; c != data_[r].size(); ++c) {
        if (c) text__ << "\t";
        text__ << data_[r][c];
      }
      text__ << "\n";
    }
    return text__.str();
  }
  virtual std::string html() const {
    std::ostringstream html__;
    html__ << "<h2>" << title() << "</h2>\n<table>\n";
    for (uint64_t r{0}; r != data_.size(); ++r) {
      html__ << "<tr>";
      for (uint64_t c{0}; c != data_[r].size(); ++c) {
        html__ << "<td>" << data_[r][c] << "</td>";
      }
      html__ << "</tr>\n";
    }
    html__ << "</table>\n";
    return html__.str();
  }
  const std::string title() const { return title_; }

 private:
  std::string title_;
  Data & data_;
};
//        font-size:calc(1vw + 3pt); }

std::string default_style() {
    return R"xxx(
body { margin:0em; padding:0em; padding-top:1.5em; padding-bottom:5em;
       font-size:calc(0.9vw + 12pt); }
h1 { text-align:center; margin: 0em 0em; font-size:150%; }
h2 { margin: 0.2em 0; font-size:120%; }
h3 { margin: 0.2em 0; font-size:110%; }
h4 { margin: 0.2em 0; font-size:100%; }
ol, ul { margin-top:0em; padding-top:0.5em; }
li { padding-bottom:0.3em; }
.clear { clear:both; }
.bold { font-weight:bold; }
.tiny { font-size:35%; }
.smaller { font-size:70%; }
.bigger { font-size:130%; }
.huge { font-size:150%; }
.huger { font-size:180%; }
.hugelist { font-size:200%; padding-left:3em; }
.lower { position:relative; top:0.25em; }
.third { display:inline; float:left; width:33%; }
.thirdplus { display:inline; float:left; width:41.66%; }
.thirdminus { display:inline; float:left; width:29.16%; }
.quarter { display:inline; float:left; width:25%; }
.half { display:inline; float:left; width:50%; }
.halfminus { display:inline; float:left; width:58.33%; }
table.header { width:100%; z-index:100;
               background-color:#DDD; font-weight:bold; font-size:70%; }
table#header { position:fixed; top:0px;
               border-bottom:4px solid #B00; 
               padding-bottom:0px; margin-bottom:0px; }
table#footer { position:fixed; bottom:0px;
               border-top:4px solid #B00; 
               padding-top:0px; margin-top:0px; }
table.header tr { padding:0em; margin:0em; }
table.header td { vertical-align:center; padding:0.25em; margin:0em; }
.leftdiv { display:block; text-align:left; font-size:220%;
           padding-left:15%; padding-right:0%; }
.left { text-align:left; }
.center { text-align:center; }
.right { text-align:right; }
.lefthead { text-align:left; width:18em; }
.righthead { text-align:right; width:18em; }
td.left td.right { white-space:nowrap; }
div.thumb { float:left; display:inline; margin:0px; padding:0px;
            width:25%; min-width:300px; max-width:100vw;
            font-size:30%; height:340px; }
div.thumb div.thumbslide { border:1px solid black; margin:10px; padding:10px;
                           height:300px; overflow:auto; }
div#body { margin:0em; padding:10px 10px 10px 10px; }
div#body.slide { margin:0em; padding:10px 2vw 10px 2vw; }
div#body p { margin:0em; padding:0.5em; }
a { text-decoration:none; color:#000; }
pre { font-size:120%; }
img { }
iframe { border:none; }
iframe.indels { border:4px solid #000; width:90vw; height:67vh; }
div.thumb iframe.indels { width:90vw; height:80vh; min-height:250px; }
iframe.spatialfull { width:90vw; height:80vh; }
iframe.spatial { width:73vh; max-width:65vw; height:70vh; }
iframe.spatialentry { width:50vh; max-width:33vw; height:70vh; }
div.thumb iframe.spatialfull { width:90%; height:70%; min-height:250px; }
div.thumb iframe.spatial { width:70%; height:80%; min-height:220px; }
div.thumb iframe.spatialentry { width:20%; height:80%; min-height:220px; 
                                max-width:150px; max-height:220px; }
.red { color:#B00; }
.bigred { font-size:150%; color:#B00; }
a.red { color:#B00; }
)xxx";
}

std::string style_extras{R"xxx(
img { }
b.chosen { color:#BB0000; }
a { text-decoration:none; color:#55F; font-weight:bold; }

tr:nth-child(even) { background-color:#f2f2f2; }
tr:hover { background-color:#ddd; }
tr#top { position:-webkit-sticky; position:sticky; top:5em; background-color:#FFF; }
td, th { padding:3px; text-align:justified; }
th { text-align:left; background-color:#4c50af; color:white; }
)xxx"};


class HTML : public PageElement {
 public:
  HTML(const std::string & content__ = "",
       const std::string & title__ = "",
       const std::string & style__ = "",
       const std::string & script__ = "") :
      title_{title__},
    style_{style__},
    script_{script__},
    content_{content__} { }
  virtual ~HTML() = default;

  virtual std::string html() const {
    std::ostringstream html_;
    html_ << "<!DOCTYPE html>\n"
          << "<html xmlns=\"http://www.w3.org/1999/xhtml\" "
          << "lang=\"en\" xml:lang=\"en\">\n"
          << "<head>\n"
          << "<meta http-equiv=\"Content-Type\""
          << " content=\"text/html; charset=utf-8\"/>\n";

    if (title().size()) html_ << "<title>" << title() << "</title>\n";
    if (style().size()) html_ << "<style>\n" << style() << "</style>\n";
    if (style().size()) html_ << "<script>\n" << script() << "</script>\n";

    html_ << "</head>\n";
    html_ << "<body onkeypress=\"navigate()\">\n" << content() << "</body>\n";
    html_ << "</html>\n";
    return html_.str();
  }
  virtual std::string text() const {
    return html();
  }
  std::string title() const { return title_; }
  std::string style() const { return style_; }
  std::string script() const { return script_; }
  std::string content() const { return content_; }

 private:
  const std::string title_;
  const std::string style_;
  const std::string script_;
  const std::string content_;
};

//
// content
//

class Slide : public PageElement {
 public:
  explicit Slide(const Title & title__,
                 const std::string & part__,
                 const Employees & people__ = Employees()) :
      title_{title__},
    part_{part__},
    people_{people__} { }
  virtual ~Slide() = default;

  virtual std::string text() const {
    std::ostringstream slide;
    slide << heading(title_.title()).text() << "\n";
    for (const PEP & pe : elements)
      slide << (*pe).text() << "\n";
    return slide.str();
  }
  virtual std::string html() const {
    std::ostringstream slide;
    if (center()) slide << "<div class=\"center\">";
    slide << heading(title_.title()).html() << "\n";
    for (const PEP & pe : elements) slide << (*pe).html() << "\n";
    if (center()) slide << "</div>";
    return slide.str();
  }
  Title title() const { return title_; }
  std::string part() const { return part_; }
  Employees people() const { return people_; }

  template <class Element>
  Slide & operator<<(const Element & element__) {
    elements.push_back(std::make_shared<Element>(element__));
    return *this;
  }
  Slide & operator<<(const char * const  text__) {
    elements.push_back(std::make_shared<Text>(std::string(text__)));
    return *this;
  }

  bool center() const { return center_; }
  Slide & center(const bool & center__) {
    center_ = center__;
    return *this;
  }

 private:
  Title title_;
  std::string part_;
  Employees people_;
  bool center_{true};
};

template<>
Slide & Slide::operator<<(const std::string & text__) {
  elements.push_back(std::make_shared<Text>(text__));
  return *this;
}

const WebRoot web_root{
  "/data/safe/paa/web/talk2018dec/", "http://wigclust19.cshl.edu/talk2018dec/"};

class Header : public PageElement {
 public:
  Header(const std::string & left__,
         const std::string & center__,
         const std::string & right__,
         const bool top__ = true) :
      left{left__},
    center{center__},
    right{right__},
    top{top__} { }
  virtual ~Header() = default;
  virtual std::string text() const {
    return left + "\n" + center + "\n" + right = "\n";
  }
  virtual std::string html() const {
    std::ostringstream header;
    header << "<table id=\"" << (top ? "header" : "footer")
           << "\" class=\"header\"><tr>\n"
           << "<td class=\"lefthead\">" << left << "</td>\n"
           << "<td class=\"center\">" << center << "</td>\n"
           << "<td class=\"righthead\">" << right
           << "</td></tr></table>\n";
    return header.str();
  }

 private:
  const std::string left;
  const std::string center;
  const std::string right;
  const bool top;
};

class Talk : public PageElement {
 public:
  Talk(const QueryVals query__,
       const std::string & name__,
       const Employee & person__,
       const Seminar & seminar__,
       const png & png__) :
      query{query__},
    name{name__},
    person{person__},
    seminar{seminar__},
    png_{png__} {
    }
  virtual ~Talk() = default;
  virtual std::string text() const {
    std::ostringstream talk;
    if (page() == "all") {
      for (const Slide & slide_ : slides) talk << slide_.text() << "\n";
    } else {
      talk << slides[int_page()].text() << "\n";
    }
    if (message.size()) talk << pre(message).text();
    return talk.str();
  }
  virtual std::string html() const {
    std::ostringstream talk;
    std::ostringstream style;
    style << default_style();
    // page decorations

    // Page links
    std::string left_link{query.query_string("page", "all")};
    std::string right_link{query.query_string("page", "all")};
    std::string tiles_link{query.query_string("page", "all")};

    // header
    std::string valid{page() == "all" || int_page() == 0 ?
          "Peter Andrews" :
          (page() + " - " + slides[int_page()].part() + " - " +
           slides[int_page()].title().short_title())};
    if (query.server_name().find("wigclust") == std::string::npos)
      valid = hlink("https://validator.w3.org/check?uri=" +
                   query.request_uri(), valid).html();
    if (page() == "all") {
      talk << Header("NYGC",
                     valid,
                     "CSHL", true).html();
      left_link = query.query_string("page", slides.size() - 1);
      right_link = query.query_string("page", 0);
    } else if (int_page() == 0) {
      talk << Header("NYGC",
                     valid,
                     ">> " +
                     hlink(right_link = query.query_string("page", 1),
                           slides[1].title().short_title()).html(),
                     true).html();
    } else {
      talk << Header(
          hlink(left_link = query.query_string("page", int_page() - 1),
                slides[int_page() - 1].title().short_title()).html() +
          " &lt;&lt;",
          valid,
          "&gt;&gt; " +
          (int_page() + 1 == slides.size() ?
           std::string("End") :
           hlink(right_link =
                 (int_page() + 1 == slides.size() ?
                  query.query_string("page", "all") :
                  query.query_string("page", int_page() + 1)),
                 slides[int_page() + 1].title().short_title()).html()),
          true).html();
    }

    // main content
    if (page() == "all") {
      talk << "<div id=\"body\">\n";
    } else {
      talk << "<div id=\"body\" class=\"slide\">\n";
    }
    if (message.size()) talk << pre(message).html();
    if (page() == "all") {
      for (uint64_t s{0}; s != slides.size(); ++s) {
        const Slide & slide_{slides[s]};
        talk << hlink(query.query_string("page", s),
                      "<div class=\"thumb\"><div class=\"thumbslide\">\n" +
                     slide_.html() + "</div></div>\n").html();
      }
    } else {
      talk << slides[int_page()].html();
    }
    // talk << "<pre style=\"clear:both;\">\n\n\n\n\n\n</pre>";
    talk << "</div>\n";


    std::ostringstream script;
    script << "document.addEventListener('keydown',\n"
           << "  function(event) {\n";
    if (left_link.size()) {
      script << "    if (event.key == \"PageUp\") {\n"
             << "      window.location = \"" << left_link << "\"\n"
             << "    }\n";
    }
    if (right_link.size()) {
      script << "    if (event.key == \"PageDown\") {\n"
             << "      window.location = \"" << right_link << "\"\n"
             << "    }\n";
    }
    script << "    if (event.key == \"0\") {\n"
           << "      window.location = \"" << tiles_link << "\"\n"
           << "    }\n";

    script << "  })\n";

    return HTML(talk.str(), seminar.title(), style.str(), script.str()).html();
  }

  Slide & slide(const Title & title__,
                const std::string & part__,
                const Employees & people__ = Employees()) {
    slides.emplace_back(title__, part__, people__);
    return slides.back();
  }
  static Allowed allowed() {
    std::set<std::string> allowed_pages;
    allowed_pages.insert("all");
    for (unsigned int p{0}; p != 100; ++p)
      allowed_pages.insert(std::to_string(p));
    const Allowed::AllowedLookup allowed_{
      {"page", {"all", allowed_pages}}
    };
    return Allowed{allowed_};
  }

  std::string page() const { return query.value("page"); }
  uint64_t int_page() const { return query.ivalue("page"); }


 private:
  const QueryVals query;
  const std::string name;
  const Employee person;
  const Seminar seminar;
  const png png_;
  std::deque<Slide> slides{};
};

//
// output
//

class OutputMethodImp {
 public:
  virtual ~OutputMethodImp() = default;
  virtual OutputMethodImp & operator<<(const std::string & text) = 0;
  virtual OutputMethodImp & operator<<(const PageElement & element) = 0;
};

enum class OutputMethodID { HTTP, FS, COUT, UI, PS };

class HTTP : public OutputMethodImp {
 public:
  static constexpr OutputMethodID omid_{OutputMethodID::HTTP};
  OutputMethodID omid() const { return omid_; }

  HTTP(std::ostream & out_,
       const std::string & format_) : out{out_}, format{format_} { }
  virtual ~HTTP() = default;
  void output_header(const std::string & format_ = "html") {
    std::map<std::string, std::string> format2content{
      {"txt", "text/plain"},
      {"cgi", "text/html"},
      {"html", "text/html"}
    };
    if (!sent_header) {
      out << "Content-Type: " << format2content.at(format_)
          << "; charset=utf-8\n";
      // out << "Cache-Control: no-cache, no-store, must-revalidate\n";
      out << "\n";
      sent_header = true;
    }
  }
  virtual OutputMethodImp & operator<<(const std::string & text) {
    output_header();
    out << text;
    return *this;
  }
  virtual OutputMethodImp & operator<<(const HTML & html__) {
    output_header();
    out << HTML(html__.html()).html();
    return (*this);
  }
  virtual OutputMethodImp & operator<<(const PageElement & element) {
    if (0) std::cerr << "format is " << format << std::endl;
    output_header(format);
    if (format == "txt") {
      out << element.text();
    } else {
      out << element.html();
    }
    return (*this);
  }

#if 0
  const std::string output_suffix{query.output_suffix()};
  if (output_suffix == "") 1;
#endif

  HTTP & operator<<(const QueryVals &) {
    output_header();
    // HTML html{out};
    return *this;
  }

 private:
  std::ostream & out;
  std::string format;
  bool sent_header{false};
};

class OutputMethod {
 public:
  OutputMethod(const QueryVals & query, std::ostream & out) {
    output = std::make_unique<HTTP>(out, query.output_suffix());
  }

  template <class Something>
  OutputMethod & operator<<(const Something & something) {
    (*output) << something;
    return *this;
  }

 private:
  OutputMethodID get_output_method() const {
    return OutputMethodID::HTTP;
  }
  std::unique_ptr<OutputMethodImp> output{};
};

template <>
OutputMethod & OutputMethod::operator<<(const HTML & something) {
  (*output) << something;
  return *this;
}

}  // namespace paa

using namespace paa;  // NOLINT

// xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

int main(int argc, char ** argv, char ** envp) {
  try {
    paa::exit_on_pipe_close();

    const InputMethod input{
      std::cin, argc, argv, envp, "/data/safe/paa/web/talk2018dec/"};

    const Allowed allowed{Talk::allowed()};
    const QueryVals query{input, allowed};

    OutputMethod output{query, std::cout};

    const Employee peter{"Peter", "Andrews", "Wigler Lab",
          "Cold Spring Harbor Laboratory", "Senior Computer Scientist"};
    const Employee mike{"Mike", "Wigler", "Wigler Lab",
          "Cold Spring Harbor Laboratory", "Professor"};
    const Employee dan{"Dan", "Levy", "Levy Lab",
          "Cold Spring Harbor Laboratory", "Associate Professor"};

    const Seminar seminar{
      "SV and CNV in Autism, CHD and Cancer",
          "Analysis and Visualization with MUMdex",
          "SNV, CNV, MUMdex",
          "Weekly&nbsp;Lab&nbsp;Meeting",
          "New&nbsp;York&nbsp;Genome&nbsp;Center",
          "Training&nbsp;Room", "NYGC",
          "December 6", "2018"};

    Talk talk{query, "talk2018dec", peter, seminar,
          png("atan.png", "atan Profile", "100%", "atan_full.png")};

    talk.slide(seminar, "Introduction", Employees{peter})
          << heading2(seminar.subtitle())
          << png("./sv_denovo.png", "SV de novo", "49%", "./sv_denovo.pdf")
          << png("./cn_denovo.png", "CN de novo", "49%")
          << heading3(seminar.date() + ", " + seminar.year())
          << type(hlink("http://mumdex.com/talk2018dec/"), "red");


    talk.slide("Abstract", "Introduction").center(false)
        << R"xxx(
<div style="text-align:justify; font-size:105%">
<p>
Structural Variation (SV) and Copy Number Variation (CNV) are responsible for many diseases and developmental syndromes, but many existing analysis pipelines have trouble identifying and combining them.
</p>
<p>
Our work examines SV and CNV using data from the autism Simons Simplex Collection (9077 people) and in-house studies on congenital heart disease (537 people) and metastatic breast cancer (202 samples from 67 people).
</p>
<p>
We developed a set of data processing and visualization tools bundled in the MUMdex software package to analyze multiple types of sequence data. The data processing tools include a maximal unique match (MUM) based aligner, a split-read based SV detector, <i>de novo</i> and transmitted mutation finders and a classical CNV pipeline. The visualization tools include a genome browser for CNV [GGraph], a MUM alignment viewer [mview], and websites to investigate copy-number profiles, <i>de novo</i> indel candidates and detected CNV events.
</p>
<p>
Our findings include a moderate increase in sensitivity and equal specificity for MUMdex <i>de novo</i> indel candidates over a standard method, the near elimination of event-size bias and a good correspondence between large MUMdex indels and CNV-based event calls.
</p>
</div>
)xxx";

    talk.slide("Previous Seminars", "Introduction").center(false)
        << type(division(
            division(
                type(division(
                    heading2("December 2015"),
                    list("MUMdex", "suffix arrays", "bridges",
                         "<i>de novo</i> search", "<i>de novo</i> pseudogenes",
                         "SMASH sequencing")), "half"),
                png("./suffix_array.png", "A Suffix Array", "40%")),
            type(heading2("February 2018"), "clear"),
            list("CNV analysis", "G-Graph",
                 "recurrence in AML: 5 people * 3 samples",
                 "indels + CN for SSC: 510 families * 4 people",
                 "pediatric SMASH CN: 126 syndromes")), "thirdplus")
        << type(division(
            png("./mumdex_format.png", "", "49%"),
            png("./bridges_read.png", "", "49%"),
            png("./pseudogene.png", "", "49%"),
            png("./smash_protocol.png", "", "49%"),
            png("./cn_primer.png", "", "49%"),
            png("./pediatric.png", "", "49%")), "halfminus")
        << division(type(heading(
            "This seminar builds on my prior work"), "center"));

    talk.slide("This Talk", "Introduction")
        << division(
            type(division(
                heading2("1. MUMdex concepts"),
                png("./bridges_read.png", "", "100%")), "half"),
            type(division(
                heading2("2. MUMdex software tools"),
                png("./mview.png", "mview application", "60%")), "half"))
        << type(division(
            heading2("3. SV and CNV discovery and visualization using MUMdex"),
            png("./indel_deletion.png", "Indel", "32%"),
            png("./CZ30.png", "", "45%")), "clear");

    talk.slide("What are MUMs and MUMdex?", "Concepts").center(false)
        << heading2("A Maximal Unique Match (MUM) "
                    "between two sequences:")
        << list("is a subsequence in common, matching exactly (match)",
                "cannot be extenced to the left or right (maximal)",
                "exists only once in each sequence (unique)")
        << png("./mum.png", "An example of a MUM", "90%")
        << heading2("MUMdex software uses MUMs for genomic analysis "
                    "and visualization:")
        << list("fast suffix array aligner outputs MUMdex alignment format",
                "<i>de novo</i> and transmitted structural variation "
                "mutation finder",
                "MUM-based traditional copy number analysis",
                "SV and CN visualization software");

    talk.slide("Bridges", "Concepts")
        << heading2("Bridges are pairs of MUMs in a read")
        << "moo";


    talk.slide("New Variant Detection", "Main")
        << heading2("I specialize in Structural and Copy Number "
                    "Variant Detection")
        << png("./indel_deletion.png", "A small deletion", "48%",
                 "./indel_deletion.pdf")
        << png("./cn_deletion.png", "A large deletion", "50%",
                 "http://mumdex.com/chd/500000/5105P/"
                 "denovo_loss.chr12.124661924-125013032/")
        << heading("SV and CNV are responsible")
        << heading("for many diseases and developmental syndromes");

    talk.slide("Structural Variant Detection", "Main")
        << png("./indel_deletion.png", "A small deletion", "60%",
                 "./indel_deletion.pdf")
        << heading2("MUMdex finds structural variants by split-read mapping");

    talk.slide({"Types of Structural Variation", "", "SV Types"},
               "Main").center(false)
        << type(division(
            type(list(
                "small insertions",
                "small deletions",
                "microsatellite instability",
                "tandem duplications",
                "large insertions (usually CNV)",
                "large deletions (usually CNV)",
                "inversions",
                "translocations"), "huge")), "thirdplus")
        << png("./wedges.png", "Types of SV events", "44%", "./wedges.png");

    talk.slide({"Split-Read Structural Variation in the SSC", "",
            "SSC SV Detection"}, "Main")
        << heading2("9079 30x coverage whole genome samples "
                    "in 2380 SSC autism families")
        << type(iframe("http://mumdex.com/indels/"), "indels")
        << heading2("33376 <i>de novo</i> SV candidates were found, "
                    "with a 95% successful validation rate");

    talk.slide({"Standard (BWA, GATK, etc) Indel Method Comparison ", "",
            "SV Detection Comparison"}, "Main")
        << png("bridge_count_compare.png", "Bridge Count Comparison", "30%")
        << png("invariant_compare.png", "Event Size Comparison", "30%")
        << png("repeat_compare.png", "Repeat Content Comparison", "30%")
        << png("parents_compare.png", "Parent Coverage Comparison", "30%")
        << png("offset_compare.png", "Anchor Offset Comparison", "30%")
        << png("max_other_compare.png", "Max Population Count Comparison",
               "30%")
        << heading2("MUMdex finds 15.5% more <i>de novo</i> indel candidates, "
                    "and they validate equally well");

    talk.slide({"Other Uses for Structural Variants", "",
            "SV Variant Utility"}, "Main")
        << heading2("Luria-Delbr&uuml;ck was used to show "
                    "that indel sequencing errors are very rare!")
        << png("./luria.png", "luria-delbruck", "50%")
        << heading2("Most indel errors are due to "
                    "machine error, not PCR error")
        << heading2("We now use indels as very sensitive probes, "
                    "many with error rates &lt; ", super("10", "-6"));

    talk.slide({"Structural Variation Visualization with MView", "",
            "SV Variant Visualization"}, "Main")
        << png("./mview.png", "mview application", "80%")
        << heading2("This MUMdex GUI is designed to visualize MUM alignments");

    talk.slide("Copy Number Variant Detection", "Main")
        << png("./cn_deletion.png", "A large deletion", "60%",
                 "./cn_deletion.png")
        << heading2("MUMdex finds CN variants by bin counting read pairs "
                    "and then segmenting");

     talk.slide({"Classical Copy Number Variation in Metastatic Breast Cancer",
             "", "CN in Breast Cancer"}, "Main")
         << heading2("Website for project uses new CN scale "
                     "to visualize highly aneuploid profiles")
         << type(iframe("http://mumdex.com/cn/?n_x=3"), "indels")
         << heading2("Features include scale choices, zoom to chromosome, "
                     "profile ordering and grouping");

     talk.slide({"SMASH Copy Number Variation in Congenital Heart Disease",
             "", "CN in CHD"}, "Main")
         << heading2("To find <i>de novo</i> and transmitted CN events in CHD")
         << type(iframe("http://mumdex.com/chd/100000/denovo.html"), "indels")
         << heading2("Most known events were found, plus many more");

     talk.slide({"A New Copy Number Display Scale",
             "", "New CN Scale"}, "Main")
         << png("./atan_fun.png", "CN function definition", "45%")
         << png("./atan.png", "New CN Scale", "80%")
         << heading2(
             "Smoothly projects CN range of 0 - ",
             lower(huge("&#x221e;")), " to graph Y axis range of 0 - 1");

     talk.slide({"G-Graph Copy Number Visualization", "", "G-Graph"}, "Main")
         << heading2("Used to easily explore copy number profiles")
         << png("./ggraph.png", "G-Graph Application", "70%")
         << heading2("G-Graph is a GUI that runs on Linux, Mac and Windows");

     talk.slide({"SV and CNV Complementarity in SSC Data", "", "SV and CNV"},
                "Main")
         << heading2("<i>De novo</i> SV events often also show up "
                     "as <i>de novo</i> CNV events, and vice-versa")
         << png("./sv_denovo.png", "SV de novo", "49%", "./sv_denovo.pdf")
         << png("./cn_denovo.png", "CN de novo", "49%")
         << heading2("~70% of large SV are found as CNV "
                     "and ~30% of CNV are found as SV")
         << heading2("Exploiting this will be one focus for me in the future");

      talk.slide({"MUMdex and Related Tools From this Seminar",
              "", "MUMdex Tools"}, "Summary").center(false)
          << type(division(png("./CZ30.png", "", "95%"),
                           png("./indel_deletion.png", "", "95%"),
                           png("./chd_view.png", "CHD website", "95%")),
                  "quarter")
          << type(division(
              heading2("MUMdex: SV, CN analysis"),
              list(hlink("http://mumdex.com/",
                         "http://mumdex.com/ for downloads and tutorial"),
                   hlink("http://mumdex.com/indels/",
                         "http://mumdex.com/indels/ for SSC analysis results"),
                   hlink("http://mumdex.com/cn/",
                         "http://mumdex.com/cn/ for html-based CN viewer")),
              heading2("SMASH: CN analysis for SMASH data"),
              list(hlink("http://mumdex.com/chd/",
                         "http://mumdex.com/chd/ for CHD CN analysis results")),
              heading2("G-Graph: CN exploratory analysis"),
              list(hlink("http://mumdex.com/ggraph/",
                         "http://mumdex.com/ggraph/ "
                         "for download and tutorial")),
              heading2("MView: MUM alignment viewer"),
              list(hlink("http://mumdex.com/",
                         "http://mumdex.com/ included with MUMdex package"))),
                  "half")
          << type(division(png("./indels_web.png", "", "100%"),
                           png("./mview.png", "mview application", "100%")),
                  "quarter");

    talk.slide("Other Current Projects", "Conclusion").center(false)
        << type(list("quantitative sensitive detection (QSD) with "
                     "single nucleotide polymorphisms (SNPs) and indels",
                     "<i>De novo</i> and transmitted indel sharing rates "
                     "in autism",
                     "sibling sharing of <i>de novo</i> and transmitted events",
                     "maximally distant barcode creation for multiplexing",
                     "machine learning for determining candidate thresholds",
                     "machine learning for judging sequence alignment validity",
                     "C++ unified display architecture for multiple output "
                     " formats (html, text, pdf, png, svg, ps, ...) and "
                     "modalities (file, http, X11, printer, app, stdout, ...) "
                     "- developed in part to present this seminar"),
                "huger");

    talk.slide("Thanks!", "Conclusion").center(false)
        << type(division(
            type(division(
                list("Mike Wigler",
                     "Dan Levy",
                     "Ivan Iossifov",
                     "Zihua Wang",
                     "Jude Kendall",
                     "Andrea Moffitt",
                     "Chris Yoon",
                     "Michael Ronemus",
                     "Matthew Wroten",
                     "Yoon-Ha Lee",
                     "Alexander Krasnitz")), "thirdminus"),
            type(division(
                list("Patty Bird",
                     "David Trimboli",
                     "Mike Riggs",
                     "Linda Rodgers",
                     "Beicong Ma",
                     "Mona Spector",
                     "Joan Alexander",
                     "Asya Stepansky",
                     "David McCandlish",
                     "Jason Sheltzer",
                     "Anders Zetterberg")), "thirdminus"),
            type(division(
                list("The Simons Foundation",
                     "The CSHL IT Department",
                     "Cold Spring Harbor Laboratory",
                     "The New York Genome Center", "Lots More!"),
                png("thanks.png", "Thanks!", "90%")),
                 "thirdplus")), "huge");
    talk.slide({"Behind the Talk Curtains", "", "Talk Secrets"},
               "Conclusion")
        << png("./emacs.png", "Talk source code", "32%", "./emacs.png")
        << png("./gemacs.png", "Graph source code", "32%", "./emacs.png")
        << png("./iemacs.png", "Graph source code", "32%", "./emacs.png")
        << heading2("Seven slides worth of source code plus the graph demo");

#if 0
    talk.slide("Talk program debug details", "Technical").center(false)
        << heading2("page = " + query.value("page"))
        << heading2("suffix = " + query.output_suffix())
        << heading2("environment:")
        << pre(input.env());
#endif

    output << talk;

    return 0;
  }

  catch(std::exception & e) {
    std::cerr << e.what() << std::endl;
    return 1;
  }
  catch(...) {
    std::cerr << "Some exception was caught." << std::endl;
    return 1;
  }
  return 0;
}
