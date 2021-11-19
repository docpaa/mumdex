//
// paastrings.h
//
// useful string functions
//
// Copyright Peter Andrews 2015 @ CSHL
//

#ifndef LONGMEM_STRINGS_H_
#define LONGMEM_STRINGS_H_

#include <set>
#include <string>
#include <utility>
#include <vector>

namespace paa {

// Some simple strings
const std::string nl{"\n"};
const std::string indent{"  "};
const std::string nli{nl + indent};
const std::string cn{"," + nl};
const std::string br{"<br />"};
const std::string brnl{br + nl};
const std::string md{" &mdash; "};
const std::string mdnl{" &mdash;" + nl};

using StringPair = std::pair<std::string, std::string>;
using StringPairs = std::vector<StringPair>;
using Strings = std::vector<std::string>;
using StringSet = std::set<std::string>;

template<class Stream>
inline std::string getline(Stream & input) {
  std::string line;
  // using std::getline;
  getline(input, line);
  return line;
}

inline std::string first_word(const std::string & text) {
  auto found = text.find(" ");
  if (found == std::string::npos) {
    return text;
  } else {
    return text.substr(0, found);
  }
}

inline Strings split(const std::string & text, const std::string & splitter) {
  size_t start{0};
  Strings result;
  while (true) {
    const size_t found{text.find(splitter, start)};
    if (found == std::string::npos) {
      result.push_back(text.substr(start));
      break;
    }
    result.push_back(text.substr(start, found - start));
    start = found + splitter.size();
  }
  return result;
}

inline std::string capitalize(std::string name) {
  name[0] = toupper(name[0]);
  return name;
}

inline std::string nlmaybe(const std::string & text) {
  return (text.size() ? nl + text : "");
}
inline std::string choose(const bool condition,
                          const std::string & text1,
                          const std::string & text2 = "") {
  return condition ? text1 : text2;
}
inline std::string maybe(const bool condition, const std::string & text) {
  return choose(condition, text);
}

// Join strings
template <class String>
inline std::string join(const std::string &, String && first) {
  return std::forward<String>(first);
}
template <class String, class ... Strings>
inline std::string join(const std::string & joiner, String && first,
                        Strings && ... rest) {
  return std::forward<String>(first) + joiner +
      join(joiner, std::forward<Strings>(rest) ...);
}

template <class ... Strings>
std::string nj(Strings && ... values) {
  return join(nl, std::forward<Strings>(values) ...);
}
template <class ... Strings>
std::string nij(Strings && ... values) {
  return join(nli, std::forward<Strings>(values) ...);
}
template <class ... Strings>
std::string cnj(Strings && ... values) {
  return join(cn, std::forward<Strings>(values) ...);
}
template <class ... Strings>
std::string cnij(Strings && ... values) {
  return join(cn + indent, std::forward<Strings>(values) ...);
}
template <class ... Strings>
std::string mdj(Strings && ... values) {
  return join(mdnl, std::forward<Strings>(values) ...);
}
template <class ... Strings>
std::string csj(Strings && ... values) {
  return join(", ", std::forward<Strings>(values) ...);
}

template <typename V>
inline std::string dash(const V & val) {
  std::ostringstream out;
  out << val;
  return out.str();
}

template <typename V, typename... Vals>
inline std::string dash(const V & val, Vals... vals) {
  return dash(val) + "-" + dash(vals...);
}

template <class Type>
inline std::string nice_string(const Type & type) {
  std::ostringstream out;
  out << type;
  return out.str();
}

// Convert to string
template <class T>
inline std::string as_string(T && val) {
  std::ostringstream out;
  out << std::forward<T>(val);
  return out.str();
}
template <class Val>
inline std::string str(Val && val) {
  return as_string(std::forward<Val>(val));
}

template <class Val>
inline uint64_t as_uint(Val && text) {
  std::istringstream text_stream{as_string(std::forward<Val>(text)).c_str()};
  uint64_t result;
  text_stream >> result;
  return result;
}

inline std::string trim(std::string input) {
  while (isspace(input.front())) input = input.substr(1);
  while (isspace(input.back())) input.pop_back();
  return input;
}

std::string to_lower(const std::string & in);
std::string to_upper(const std::string & in);
// void remove(std::string & input, const std::string & search);
inline void replace_inplace(std::string & input, const char a, const char b) {
  size_t pos = 0;
  while ((pos = input.find(a, pos)) != std::string::npos) {
    input[pos] = b;
    ++pos;
  }
}
inline std::string replace(const std::string & input,
                           const char a, const char b) {
  size_t pos = 0;
  std::string copy{input};
  while ((pos = copy.find(a, pos)) != std::string::npos) {
    copy[pos] = b;
    ++pos;
  }
  return copy;
}

inline std::string replace_substring(const std::string & str,
                                     const std::string & oldstr,
                                     const std::string & newstr) {
  std::string result{str};
  const size_t pos = result.find(oldstr);
  if (pos != std::string::npos) {
    result.replace(pos, oldstr.size(), newstr);
  }
  return result;
}

inline std::string replace_one(const std::string & str,
                               const std::string & oldstr,
                               const std::string & newstr) {
  std::string result{str};
  const size_t pos = result.find(oldstr);
  if (pos != std::string::npos) {
    result.replace(pos, oldstr.size(), newstr);
  }
  return result;
}

inline std::string replace_all(std::string result,
                               const std::string & oldstr,
                               const std::string & newstr) {
  size_t pos = 0;
  while ((pos = result.find(oldstr, pos)) != std::string::npos) {
    result.replace(pos, oldstr.size(), newstr);
    pos += newstr.size();
  }
  return result;
}

inline void replace_all_inplace(std::string & result,
                                const std::string & oldstr,
                                const std::string & newstr) {
  size_t pos = 0;
  while ((pos = result.find(oldstr, pos)) != std::string::npos) {
    result.replace(pos, oldstr.size(), newstr);
    pos += newstr.size();
  }
}

inline void replace_substring_inplace(std::string & str,
                                      const std::string & oldstr,
                                      const std::string & newstr) {
  const size_t pos = str.find(oldstr);
  if (pos != std::string::npos) {
    str.replace(pos, oldstr.size(), newstr);
  }
}

inline std::string remove_substring(const std::string & str_,
                                           const std::string & substr) {
  std::string str{str_};
  const size_t pos = str.find(substr);
  if (pos != std::string::npos) {
    str.replace(pos, substr.size(), "");
  }
  return str;
}
inline void remove_substring_inplace(std::string * str_,
                                     const std::string & substr) {
  const size_t pos = str_->find(substr);
  if (pos != std::string::npos) {
    str_->replace(pos, substr.size(), "");
  }
}
inline std::string remove_including_final(const std::string & str_,
                                          const char c) {
  std::string str{str_};
  const size_t pos = str.find_last_of(c);
  if (pos != std::string::npos) {
    str = str.substr(pos + 1);
  }
  return str;
}

inline std::string remove_including_initial(const std::string & str_,
                                            const char c) {
  std::string str{str_};
  const size_t pos = str.find_first_of(c);
  if (pos != std::string::npos) {
    str.resize(pos);
  }
  return str;
}

inline void replace_substring_inplace_c(std::string & str,
                                        const char * const oldstr,
                                        const char * const newstr) {
  return replace_substring_inplace(
      str, std::string(oldstr), std::string(newstr));
}

}  // namespace paa

#endif  // LONGMEM_STRINGS_H_
