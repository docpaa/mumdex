//
// html.h
//
// useful html functions
//
// Copyright Peter Andrews 2021 @ CSHL
//

#ifndef UTILITY_HTML_H_
#define UTILITY_HTML_H_

#include <map>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "paastrings.h"
#include "utility.h"

namespace paa {

const std::string nb{"&nbsp;"};

// Warnings and errors
void warning(const std::string & message) {
  std::cout << "<p class=\"warning\">" << message << "</p>";
}
void success(const std::string & message) {
  std::cout << "<p class=\"success\">" << message << "</p>";
}

// Encode/decode text for attribute, html, url display
std::string encode(std::string text) {
  replace_all_inplace(text, "&", "&amp;");
  replace_all_inplace(text, "<", "&lt;");
  replace_all_inplace(text, ">", "&gt;");
  replace_all_inplace(text, "\"", "&quot;");
  replace_all_inplace(text, "'", "&#x27;");
  replace_all_inplace(text, "/", "&#x2F;");
  return text;
}
std::string attribute_encode(std::string attribute) {
  attribute = encode(attribute);
  replace_all_inplace(attribute, nl, "&#10;");
  return attribute;
}

// Decode strings with hex coding like %A7 to the correct characters
std::string url_decode(const std::string & text) {
  std::istringstream text_stream{text.c_str()};
  char c;
  char h;
  char l;
  std::string result;
  std::string hex;
  while (text_stream.get(c)) {
    if (c == '%') {
      if (text_stream.get(h)) {
        if (isxdigit(h) && text_stream.get(l)) {
          if (isxdigit(l)) {
            hex = {h, l};
            c = stoi(hex, nullptr, 16);
            result += c;
          } else {
            result += c;
            result += h;
            result += l;
          }
        } else {
          result += c;
          result += h;
        }
      } else {
        result += c;
      }
    } else {
      result += c;
    }
  }
  return result;
}
void test_url_decode() {
  const std::string chars{"`~!@#$%^&*()-_=+[{]};:'\",<.>/?"};
  std::cout << "<h2>Character encoding</h2>" << nl << "<ul>";
  for (const char c : chars) {
    std::ostringstream hex_stream;
    hex_stream << '%' << std::setfill('0') << std::setw(2) << std::hex
               << uint64_t(c);
    const std::string hex{hex_stream.str()};
    const std::string C(1, c);
    std::cout << "<li>" << encode(C) << " -> " << hex
         << " -> " << encode(url_decode(hex)) << "</li>" << nl;
  }
  std::cout << "</ul>" << nl;
  const std::string test{"abcde%28fghi%29jlk%2dlm%no"};
  std::cout << "<p>" << url_decode(test) << "</p>" << nl;
}

// Server environment variables
const std::string server_name{get_environment("SERVER_NAME", "localhost")};
const std::string remote_address{get_environment("REMOTE_ADDR", "command")};
const std::string request_uri{get_environment("REQUEST_URI", "/covid/")};
const bool served{remote_address != "command"};

// Parse a POST request
using StringSet = std::set<std::string>;
using StringMap = std::map<std::string, std::string>;
using StringsMap = std::map<std::string, std::vector<std::string>>;
using StringMapValue = std::pair<const std::string, std::string>;
StringMap parse_post() {
  StringMap result;
  if (!served) return result;
  std::string boundary;
  getline(std::cin, boundary, '\r');
  if (!std::cin) return result;
  std::cin.get();
  std::string line;
  std::string data;
  char c;
  while (std::cin.get(c)) data += c;
  size_t pos{0};
  size_t next_pos;
  while ((next_pos = data.find(boundary, pos)) != std::string::npos) {
    const std::string block{data.substr(pos, next_pos - pos - 2)};
    std::istringstream block_stream{block.c_str(), std::ios::binary};
    std::string name;
    while (getline(block_stream, line)) {
      if (line == "\r") break;
      const std::string name_start{" name=\""};
      const size_t name_pos{line.find(name_start)};
      if (name_pos != std::string::npos) {
        std::istringstream name_stream{
          line.substr(name_pos + name_start.size())};
        getline(name_stream, name, '"');
      }
      const std::string file_start{"filename=\""};
      const size_t file_pos{line.find(file_start)};
      if (file_pos != std::string::npos) {
        std::istringstream file_stream{
          line.substr(file_pos + file_start.size())};
        std::string file_name;
        getline(file_stream, file_name, '"');
        result[name + "_name"] = file_name;
      }
    }
    if (name.empty())
      throw Error("Could not find name in parse_post in ") << block;
    const size_t dpos{block.find("\r\n\r\n")};
    std::string & name_data{result[name]};
    name_data = block.substr(dpos + 4);
    pos = next_pos + boundary.size() + 2;
  }
  return result;
}

StringMap parse_cookies() {
  StringMap result;
  const std::string cookies_string{get_environment("HTTP_COOKIE")};
  std::istringstream cookies_stream{cookies_string.c_str()};
  std::string cookie;
  while (getline(cookies_stream, cookie, ';')) {
    std::istringstream cookie_stream{cookie.c_str()};
    if (cookie_stream.peek() == ' ') cookie_stream.get();
    std::string value;
    getline(cookie_stream, cookie, '=');
    getline(cookie_stream, value);
    if (cookie_stream) {
      result[cookie] = value;
      // result["message"] += cookie + "=" + value + ";";
    }
  }
  return result;
}

// Get input from a GET query string and POST
StringMap query() {
  const std::string query_string{get_environment("QUERY_STRING")};
  std::istringstream query_stream{query_string.c_str()};
  std::string token;
  std::vector<std::string> query_parts;
  while (getline(query_stream, token, '&'))
    if (token.size())
      query_parts.push_back(token);
  // if (query_parts.empty()) throw Error("No query string available");
  if (false) {
    std::cout << "<h2>Query parts</h2>" << nl;
    std::cout << "<ul>" << nl;
    for (const std::string & part : query_parts)
      std::cout << "<li>" << part << "</li>" << nl;
    std::cout << "</ul>" << nl;
  }
  StringMap query_keys;
  StringMap cookies{parse_cookies()};
  for (StringMapValue & item : cookies)
    query_keys[item.first] = move(item.second);
  for (const std::string & part : query_parts) {
    std::istringstream part_stream{part.c_str()};
    std::string key;
    std::string value;
    getline(part_stream, key, '=');
    if (part_stream) {
      getline(part_stream, value);
      if (part_stream) {
        replace_all_inplace(value, "+", " ");
        value = url_decode(value);
        if (value.size()) query_keys[key] = value;
      }
    }
  }
  StringMap post{parse_post()};
  for (StringMapValue & item : post)
    query_keys[item.first] = move(item.second);
  if (query_keys.count("action") && query_keys["action"] == "logout") {
    query_keys.erase("action");
    query_keys.erase("user");
    query_keys.erase("password");
  }
  if (false) {
    std::cout << "<h2>Key values</h2>" << nl;
    std::cout << "<dl>" << nl;
    for (const auto & kv : query_keys)
      std::cout << "<dt>" << kv.first << "</dt><dd>" << kv.second << "</dd>"
                << nl;
    std::cout << "</dl>" << nl;
  }
  return query_keys;
}

// HTML attributes
inline std::string attrs(const StringPairs & attributes) {
  static const StringSet boolean_attr{
    "required", "hidden", "selected", "checked", "disabled", "multiple"};
  std::string result;
  for (const StringPair & attr : attributes) {
    result += " " + attr.first;
    if (!boolean_attr.count(attr.first)) result += "=\"" + attr.second + "\"";
  }
  return result;
}
template <class ... Attrs>
std::string attrs(const StringPairs & first, Attrs && ... more) {
  return attrs(first) + attrs(std::forward<Attrs>(more) ...);
}
std::string attrs() {
  return "";
}

// Common attributes
const StringPair nodisplay{"style", "display: none"};
const StringPair required{"required", ""};
const StringPair disabled{"disabled", ""};
const StringPair multiple{"multiple", ""};
const StringPair hidden{"hidden", ""};
const StringPair selected{"selected", ""};
const StringPair checked{"checked", ""};
const StringPair small{"style", "width: 3em"};
const StringPair checkbox_t{"type", "checkbox"};
const StringPair file_t{"type", "file"};
const StringPair tel_t{"type", "tel"};
const StringPair tel_pattern{"pattern", "\\([0-9]{3}\\) [0-9]{3} - [0-9]{4}"};
const StringPair tel_format{"title", "Required format: (NNN) NNN-NNNN"};
const StringPairs tel_attrs{required, tel_t, tel_pattern, tel_format};
const StringPair uint_pattern{"pattern", "[1-9][0-9]*"};
const StringPair uint_format{"title", "Enter a positive integer"};
const StringPairs uint_attrs{small, required, uint_pattern, uint_format};
const std::string date_string{"2[0-9]{3}-[0-1][0-9]-[0-3][0-9]"};
const StringPair date_pattern{"pattern", date_string};
const StringPair date_format{"title", "Required format: YYYY-MM-DD"};
const StringPairs date_attrs{required, date_pattern, date_format};
const StringPair time_pattern{"pattern", date_string +
  " [0-2][0-9]:[0-6][0-9]:[0-6][0-9]"};
const StringPair time_format{"title", "Required format: YYYY-MM-DD HH:MM:SS"};
const StringPairs time_attrs{required, time_pattern, time_format};
const StringPair word_pattern{"pattern", "[A-Za-z]+"};
const StringPair word_format{"title", "Only letters are allowed"};
const StringPairs word_attrs{required, word_pattern, word_format};
const StringPair name_pattern{"pattern", "[A-Za-z ]+"};
const StringPair name_format{"title", "Only letters and spaces are allowed"};
const StringPairs name_attrs{required, name_pattern, name_format};
const StringPair name_num_pattern{"pattern", "[A-Za-z0-9 ]+"};
const StringPair name_num_format{"title",
  "Only letters, numbers and spaces are allowed"};
const StringPairs name_num_attrs{required, name_num_pattern, name_num_format};
const StringPair places_pattern{"pattern", "[A-Za-z0-9, ]+"};
const StringPair places_format{"title",
  "Only letters, numbers, spaces, and commas are allowed"};
const StringPairs places_attrs{required, places_pattern, places_format};
const StringPair password_t{"type", "password"};
const StringPairs password_attrs{required, password_t, {"id", "password"}};
const StringPair email_t{"type", "email"};
const StringPairs email_attrs{required, email_t};
const std::string font_pattern{"[0-9]{1,4}(.[0-9]{1,4})?(pt|%|em|in|px|pc)"};
const StringPairs font_attrs{
  required, {"pattern", font_pattern + "(/" + font_pattern + ")?(NU)?"},
  {"title", "Like: " + font_pattern}};

// HTML tag
class Tag {
 public:
  explicit Tag(const std::string & name_,
               const std::string & other_ = "",
               const StringPairs & default_attributes_ = {}) :
      name{name_}, other{other_}, default_attributes{default_attributes_} {}
  std::string operator()(const std::string & content,
                         const std::string & other_,
                         const StringPairs & attrs_ = {},
                         const bool end = true) const {
    if (other_.size() && other.empty())
      throw Error("Bad other tag usage") << other << other_;
    return "<" + name +
        (other.size() && other_.size() ?
         " " + other + "=\"" + other_ + "\"" : "") +
        attrs(default_attributes) + attrs(attrs_) + ">" +
        content + (end ? "</" + name + ">" : "");
  }
  template<class String>
  std::string operator()(const String & content,
                         const StringPairs & attrs_ = {},
                         const bool end = true) const {
      return (*this)(content, "", attrs_, end);
  }
  std::string operator()(const StringPairs & attrs_,
                         const bool end) const {
      return (*this)("", "", attrs_, end);
  }
  template<class String>
  std::string operator()(const String & content, const bool end) const {
    return (*this)(content, "", {}, end);
  }
  std::string operator()(const bool end) const {
    return (*this)("moo", "", {}, end);
  }
  template<class String1, class String2>
  std::string operator()(const String1 & content,
                         const String2 & other_,
                         const bool end = true) const {
    return (*this)(content, other_, {}, end);
  }

 private:
  std::string name;
  std::string other;
  StringPairs default_attributes;
};

// HTML tags
const Tag a{"a", "href"};
const Tag button("button", "", {{"type", "button"}});
const Tag div("div");
const Tag label("label");
const Tag h1("h1");
const Tag h2("h2");
const Tag h3("h3");
const Tag Input("input");
const Tag li("li");
const Tag option("option");
const Tag p("p");
const Tag pre("pre");
const Tag script("script");
const Tag span("span");
const Tag style("style");
const Tag table_tag("table");
const Tag td("td");
const Tag textarea("textarea");
const Tag th("th");
const Tag title("title");
const Tag tr("tr");
const Tag Select("select");
const Tag submit{"button", "name", {{"type", "submit"}}};

// Special html tags
inline std::string input(const std::string & name, StringPairs attrs = {}) {
  return Tag("input", "", {{"name", name}})(attrs, false);
}
inline std::string input(const std::string & name, const std::string & value,
                         const StringPairs & attrs = {}) {
  return Tag("input", "", {{"name", name}, {"value", value}})(attrs, false);
}
inline std::string input(const StringPairs & attrs = {}) {
  return Tag("input")(attrs, false);
}

// Other
inline std::string add(const std::string & type,
                       const std::string & more = "") {
  return nj(label(capitalize(type) + " name:"),
            input("add_" + type, name_attrs) + nlmaybe(more),
            submit("Add " + type));
}
inline std::string change(const std::string & name, const std::string & type) {
  return nj(label(capitalize(type) + " name:"),
            input("change_name", name,
                  (type == "place") ? name_num_attrs : name_attrs),
            submit("Change name"));
}
StringPairs anchor(const std::string & name) {
  return {{"id", name + "_a"}, {"class", "a"}};
}

// To remove
const std::string a_base{"<a href=\"?"};

const StringMap query_keys{query()};
std::string query(const std::string & key,
                  const std::string & default_value = "") {
  auto value_iter = query_keys.find(key);
  if (value_iter == query_keys.end()) return default_value;
  return value_iter->second;
}

// Page view and action
const std::string overview{"Overview"};
const std::string view{query("view", overview)};
const std::string action{query("action")};
const bool defaults{action == "defaults"};
const bool reset_{action == "reset"};
const bool reset{action == "reset!"};
bool reload{action == "reload"};

// Construct a new query string
std::string query(const StringPairs & key_vals = {}) {
  std::string result{"?"};
  StringSet seen;
  std::string anchor;
  uint64_t n{0};
  for (uint64_t i{0}; i != key_vals.size(); ++i) {
    const StringPair & key_val{key_vals[i]};
    if (key_val.first == "anchor") {
      anchor = key_val.second + "_a";
      continue;
    }
    if (n++) result += "&";
    if (key_val.first != "view" || key_val.second != overview)
      result += key_val.first + "=" + key_val.second;
    seen.insert(key_val.first);
  }
  if (!seen.count("view") && view.size() && view != overview)
    result += std::string(n ? "&" : "") + "view=" + view;
  if (anchor.size()) result += "#" + anchor;
  if (result.size() == 1) result = "./";
  return result;
}


template <class Type>
Type dc(const std::string & name, const Type & def) {
  return defaults ? def : query(name, def);
}
std::string dc(const std::string & name, const char * const def) {
  return defaults ? def : query(name, def);
}
bool dc(const std::string & name, const bool def) {
  return defaults ? def :
      query(name, std::string(def ? "true" : "false")) == "true";
}

// The start and end of a form
std::string sform(const StringPairs & pieces = {},
                  const StringPairs & attributes = {}) {
  std::string result{"<form method=\"post\" enctype=\"multipart/form-data\""
                " action=\"" + query(pieces) + "\""};
  for (const auto & item : attributes)
    result += " " + item.first + "=\"" + item.second + "\"";
  result += ">";
  return result;
}
const std::string eform{"</form>"};

// A complete form
// Join strings
inline std::string form_pieces() {
  return eform;
}
template <class ... Strings>
std::string form_pieces(const std::string & first, Strings && ... rest) {
  return first + nl + form_pieces(std::forward<Strings>(rest) ...);
}
template <class ... Strings>
void form(const StringPairs & pieces, const StringPairs & attributes,
          const std::string & first, Strings && ... rest) {
  std::string attr;
  for (const auto & item : attributes)
    attr += " " + item.first + "=\"" + item.second + "\"";
  std::cout << nj("<form method=\"post\" enctype=\"multipart/form-data\""
                  " action=\"" + query(pieces) + "\"" + attr + ">", first,
                  form_pieces(std::forward<Strings>(rest) ...), "");
}
template <class ... Strings>
void form(const StringPairs & pieces, Strings && ... rest) {
  form(pieces, {}, std::forward<Strings>(rest) ...);
}
template <class ... Strings>
void form(Strings && ... rest) {
  form({}, {}, std::forward<Strings>(rest) ...);
}

// Delete button
std::string del(const bool permission, const std::string & delete_type,
                const StringPairs & others = {}) {
  if (!permission) return "";
  const std::string did{delete_type + "_id"};
  std::string result{
    "," + nl + "concat('" + sform(
        {{"view", view}, {"delete_" + delete_type, "', " + did + ", '"}}) +
    submit("Delete", {{"class", "table_button"}})};
  for (const StringPair & other : others)
    result += input(other.first, other.second, {hidden});
  result += eform + "') as 'Delete" + br + delete_type + "'";
  return result;
}

// Site hyperlink
std::string link(const StringPairs & pieces, const std::string & text,
                 const StringPairs & attributes = {}) {
  return a(text, query(pieces), attributes);
}

// Grant access to a function
template <class Fun>
void control(Fun & fun, const bool permission) {
  if (permission) {
    fun();
  } else {
    warning("You require additional privileges to " +
            (action.size() ? "perform the action '" + action + "'" :
             "view the contents of the " + view + " page"));
  }
}

}  // namespace paa

#endif  // UTILITY_HTML_H_
