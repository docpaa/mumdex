//
// lpe.h
//
// lexer, parser, evaluator
//
// Copyright 2020 Peter Andrews @ CSHL
//

#ifndef PAA_UTILITY_LPE_H_
#define PAA_UTILITY_LPE_H_

#include <algorithm>
#include <cmath>
#include <exception>
#include <functional>
#include <iostream>
#include <map>
#include <regex>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace paa {

// want:
//
// unary, binary, ternary operators
// functions of 0-N arguments
// literal numbers
// symbols representing numbers
// constants
//
// simple tokenization
//
// left-to-right parsing of standard infix notation
// representation as infix, rpn, or parse tree
// reduction of expressions

// Type of token
enum class TokenType {
  end,
  open,
  close,
  comma,
  number,
  name,
  Operator,
  function
};

// Token text name
const char * name(const TokenType type) {
  switch (type) {
    case TokenType::end:
      return "end";
    case TokenType::open:
      return "open";
    case TokenType::close:
      return "close";
    case TokenType::comma:
      return "comma";
    case TokenType::number:
      return "number";
    case TokenType::name:
      return "name";
    case TokenType::Operator:
      return "operator";
    case TokenType::function:
      return "function";
    default:
      throw Error("Unknown TokenType in name");
  }
}

// Type of function
enum class FunctionType {
  none,
  add,
  subtract,
  multiply,
  divide,
  mod,
  pow,
  max,
  max3,
  max4,
  max5,
  min,
  min3,
  min4,
  min5,
  sqrt
};

// FunctionType name lookup
const std::map<std::string, FunctionType> function_types{
  {"max", FunctionType::max},
  {"max3", FunctionType::max3},
  {"max4", FunctionType::max4},
  {"max5", FunctionType::max5},
  {"min", FunctionType::min},
  {"min3", FunctionType::min3},
  {"min4", FunctionType::min4},
  {"min5", FunctionType::min5},
  {"sqrt", FunctionType::sqrt},
  {"pow", FunctionType::pow},
};

// Function text name
const char * name(const FunctionType type) {
  switch (type) {
    case FunctionType::none:
      return "none";
    case FunctionType::add:
      return "add";
    case FunctionType::subtract:
      return "subtract";
    case FunctionType::multiply:
      return "multiply";
    case FunctionType::divide:
      return "divide";
    case FunctionType::mod:
      return "mod";
    case FunctionType::pow:
      return "pow";
    case FunctionType::max:
      return "max";
    case FunctionType::max3:
      return "max3";
    case FunctionType::max4:
      return "max4";
    case FunctionType::max5:
      return "max5";
    case FunctionType::min:
      return "min";
    case FunctionType::min3:
      return "min3";
    case FunctionType::min4:
      return "min4";
    case FunctionType::min5:
      return "min5";
    case FunctionType::sqrt:
      return "sqrt";
    default:
      throw Error("Unknown FunctionType in name");
  }
}

// Function operator precedence
unsigned int precedence(const FunctionType type) {
  switch (type) {
    case FunctionType::none:
      throw Error("Disallowed none FunctionType in precedence");
    case FunctionType::add:
    case FunctionType::subtract:
      return 1;
    case FunctionType::multiply:
    case FunctionType::divide:
    case FunctionType::mod:
      return 2;
    case FunctionType::pow:
    case FunctionType::max:
    case FunctionType::max3:
    case FunctionType::max4:
    case FunctionType::max5:
    case FunctionType::min:
    case FunctionType::min3:
    case FunctionType::min4:
    case FunctionType::min5:
    case FunctionType::sqrt:
      return 100;
    default:
      throw Error("Unknown FunctionType in precedence");
  }
}

// Token regular expressions
class TokenRegex {
 public:
  TokenRegex(const std::string & regex__, const TokenType type__) :
      regex_{regex__}, type_{type__} {}
  const std::regex & regex() const { return regex_; }
  TokenType type() const { return type_; }

 private:
  std::regex regex_;
  TokenType type_;
};
using TokenRegexes = std::vector<TokenRegex>;
const TokenRegexes token_regexes{
  {R"xxx(^[+-]?\d+\.?\d*(?:[eE][+-]?\d+)?)xxx", TokenType::number},
  {R"xxx(^(?:m(?:ax|in)\d?|sqrt|pow)\b)xxx", TokenType::function},
  {R"xxx(^[[:alpha:]]\w*)xxx", TokenType::name}
};

// Any token
class Token {
 public:
  Token(const std::string & symbol__,
        const TokenType token_type__,
        const FunctionType function_type__) :
      symbol_{symbol__},
      token_type_{token_type__},
      function_type_{function_type__} { }
  const std::string & symbol() const { return symbol_; }
  TokenType token_type() const { return token_type_; }
  FunctionType function_type() const { return function_type_; }

 private:
  std::string symbol_;
  TokenType token_type_;
  FunctionType function_type_;
};
std::ostream & operator<<(std::ostream & out, const Token & token) {
  return out << token.symbol();
}

// Simple tokenizer
using CI = std::string::const_iterator;
Token get_token(CI & ci, const CI end, const TokenType last) {
  // Remove leading space
  while (ci != end && isspace(*ci)) ++ci;

  // Return if empty
  if (ci == end) return Token{"", TokenType::end, FunctionType::none};

  // Return if comment
  if (*ci == '#') {
    ci = end;
    return Token{"", TokenType::end, FunctionType::none};
  }

  // Look for single character operators
  switch (*ci) {
    case '(':
      ++ci;
      return Token{"(", TokenType::open, FunctionType::none};
    case ')':
      ++ci;
      return Token{")", TokenType::close, FunctionType::none};
    case ',':
      ++ci;
      return Token{",", TokenType::comma, FunctionType::none};
    case '+':
    case '-':
      // Allow for unary +/- for all not following name, number, close
      switch (last) {
        case TokenType::close:
        case TokenType::number:
        case TokenType::name:
          if (*ci == '+') {
            ++ci;
            return Token{"+", TokenType::Operator, FunctionType::add};
          } else {
            ++ci;
            return Token{"-", TokenType::Operator, FunctionType::subtract};
          }
        case TokenType::end:
        case TokenType::open:
        case TokenType::comma:
        case TokenType::Operator:
        case TokenType::function:
        default:
          break;
      }
      break;
    case '*':
      ++ci;
      return Token{"*", TokenType::Operator, FunctionType::multiply};
    case '/':
      ++ci;
      return Token{"/", TokenType::Operator, FunctionType::divide};
    case '%':
      ++ci;
      return Token{"%", TokenType::Operator, FunctionType::mod};
    case '^':
      throw Error("Power operator ^ not implemented");
      ++ci;
      return Token{"^", TokenType::Operator, FunctionType::pow};
    default:
      break;
  }

  // Look for a token defined by a regex
  for (const TokenRegex & regex : token_regexes) {
    std::smatch matches;
    if (regex_search(ci, end, matches, regex.regex())) {
      const FunctionType function_type{
        regex.type() == TokenType::function ?
            function_types.at(matches.str()) : FunctionType::none};
      ci += matches.str().size();
      return Token{matches.str(), regex.type(), function_type};
    }
  }

  // An error if string not interpreted by above rules
  throw Error("Parse error beginning at:") << std::string(ci, end);
}

// Get token vector from expression string
using Tokens = std::vector<Token>;
Tokens get_tokens(const std::string & expression) {
  std::string::const_iterator pos{expression.begin()};
  Tokens tokens;
  TokenType last{TokenType::end};
  while (pos != expression.end()) {
    Token token{get_token(pos, expression.end(), last)};
    if ((last = token.token_type()) == TokenType::end) break;
    tokens.push_back(std::move(token));
  }
  return tokens;
}

using Words = std::vector<std::string>;

// Compare tokens with words
bool operator !=(const Tokens & tokens, const Words & words) {
  if (tokens.size() != words.size()) return true;
  for (uint64_t t{0}; t != tokens.size(); ++t)
    if (tokens[t].symbol() != words[t]) return true;
  return false;
}

// Test effect of spacing on tokenization
void test_tokenization(const std::string & expression,
                       const bool must_succeed = true) try {
  // Get space separated tokens
  const Words tokens{[&expression]() {
      std::istringstream stream{expression.c_str()};
      std::string token;
      Words result;
      while (stream >> token) result.push_back(token);
      return result;
    }()};
  if (tokens.size() < 1)
    throw Error("Unexpected too few tokens") << tokens.size();

  // Construct all possible expressions with spaces or not between tokens
  // And test tokenization on each version
  const uint64_t n_spacings{1ul << (tokens.size() + 1)};
  for (uint64_t spacing{0}; spacing != n_spacings; ++spacing) {
    uint64_t reduced{spacing};
    std::string rewrite{""};
    for (uint64_t token{0}; token != tokens.size(); ++token) {
      if (reduced % 2 == 1) rewrite += ' ';
      reduced /= 2;
      rewrite += tokens[token];
    }
    if (reduced % 2 == 1) rewrite += ' ';
    // std::cerr << spacing << ": " << rewrite << std::endl;
    const Tokens & tokens2{get_tokens(rewrite)};
    if (tokens2 != tokens)
      throw Error("Tokens mismatch for expression\n")
          << "\n" << tokens << "\n" << tokens2;
  }
  if (must_succeed) {
    std::cout << "Good tokenization for: " << expression << std::endl;
  } else {
    std::cerr << "BAD tokenization success for: " << expression << std::endl;
  }
} catch (Error & err) {
  if (must_succeed) throw Error(err.what()) << "in" << expression;
  std::cout << "Good tokenization failure for: " << expression << std::endl;
}

// Holds a number, either constant or not
class Number {
 public:
  Number(const double number__, const bool is_constant__ = false) :  // NOLINT
      number_{number__}, is_constant_{is_constant__} {}  // NOLINT
  operator double() const { return number_; }
  double * pointer() { return &number_; }
  const double * pointer() const { return &number_; }
  bool is_constant() const { return is_constant_; }
  void unchecked_set(const double number__) { number_ = number__; }
  void set(const double number__) {
    if (is_constant_) throw Error("Attempt to set a constant");
    unchecked_set(number__);
  }

 private:
  double number_;
  bool is_constant_;
};
using Symbols = std::map<std::string, Number>;

class Value {
 public:
  Value() {}
  explicit Value(const double value) {
    std::ostringstream stream;
    stream << value;
    auto emplaced = literals_.emplace(stream.str(), Number{value, true});
    number_ = &emplaced.first->second;
  }
  explicit Value(const std::string & literal) {
    std::istringstream stream{literal.c_str()};
    double value;
    stream >> value;
    if (!stream)
      throw Error("Problem converting literal") << literal << "to double";
    std::string rest;
    getline(stream, rest);
    if (rest.size())
      throw Error("Extra text found at end of literal") << literal;
    auto emplaced = literals_.emplace(literal, Number{value, true});
    number_ = &emplaced.first->second;
  }
  Value(Symbols & symbols__, const std::string & symbol,
        const bool ok_undefined) {
    // std::cerr << "Symbol in expression" << std::endl;
    if (ok_undefined) {
      auto emplaced = symbols__.emplace(symbol, 0.0);
      number_ = &emplaced.first->second;
    } else {
      try {
        number_ = &symbols__.at(symbol);
      } catch (...) {
        throw Error("Undefined symbol") << symbol;
      }
    }
  }
  Value(Symbols & symbols__,
        const std::string & symbol,
        const double value,
        const bool is_constant_ = false) {
    auto emplaced = symbols__.emplace(symbol, Number{value, is_constant_});
    if (!emplaced.second) throw Error("Redefinition of symbol") << symbol;
    number_ = &emplaced.first->second;
  }
  Value(const std::string & symbol,
        const double value,
        const bool is_constant_ = false) :
      Value{symbols_, symbol, value, is_constant_} {}
  operator double() const { return *number_; }
  bool is_constant() const { return number_->is_constant(); }

 private:
  static Symbols literals_;
  static Symbols symbols_;
  Number * number_{nullptr};
};
Symbols Value::literals_;
Symbols Value::symbols_;

class NumFun;

using FlexArgs = std::vector<double>;
using FlexFunction = std::function<double (FlexArgs &)>;

// Function also handles operators
inline double dummy(FlexArgs &) { return 0.0; }
class Function {
 public:
  Function() {}
  explicit Function(const FunctionType function_type) {
    switch (function_type) {
      case FunctionType::add:
        n_args_ = 2;
        function = [](FlexArgs & args) {
          const double second{args.back()};
          args.pop_back();
          const double first{args.back()};
          args.pop_back();
          return first + second;
        };
        break;
      case FunctionType::subtract:
        n_args_ = 2;
        function = [](FlexArgs & args) {
          const double second{args.back()};
          args.pop_back();
          const double first{args.back()};
          args.pop_back();
          return first - second;
        };
        break;
      case FunctionType::multiply:
        n_args_ = 2;
        function = [](FlexArgs & args) {
          const double second{args.back()};
          args.pop_back();
          const double first{args.back()};
          args.pop_back();
          return first * second;
        };
        break;
      case FunctionType::divide:
        n_args_ = 2;
        function = [](FlexArgs & args) {
          const double second{args.back()};
          args.pop_back();
          const double first{args.back()};
          args.pop_back();
          return first / second;
        };
        break;
    case FunctionType::mod:
      n_args_ = 2;
      function = [](FlexArgs & args) {
        const double second{args.back()};
        args.pop_back();
        const double first{args.back()};
        args.pop_back();
        return static_cast<uint64_t>(first) % static_cast<uint64_t>(second);
      };
      break;
    case FunctionType::pow:
      n_args_ = 2;
      function = [](FlexArgs & args) {
        const double second{args.back()};
        args.pop_back();
        const double first{args.back()};
        args.pop_back();
        return pow(first, second);
      };
      break;
    case FunctionType::max:
      n_args_ = 2;
      function = [](FlexArgs & args) {
        const double second{args.back()};
        args.pop_back();
        const double first{args.back()};
        args.pop_back();
        return std::max(first, second);
      };
      break;
    case FunctionType::max3:
      n_args_ = 3;
      function = [](FlexArgs & args) {
        const double third{args.back()};
        args.pop_back();
        const double second{args.back()};
        args.pop_back();
        const double first{args.back()};
        args.pop_back();
        return std::max(std::max(first, second), third);
      };
      break;
    case FunctionType::max4:
      n_args_ = 4;
      function = [](FlexArgs & args) {
        const double fourth{args.back()};
        args.pop_back();
        const double third{args.back()};
        args.pop_back();
        const double second{args.back()};
        args.pop_back();
        const double first{args.back()};
        args.pop_back();
        return std::max(std::max(std::max(first, second), third), fourth);
      };
      break;
    case FunctionType::max5:
      n_args_ = 5;
      function = [](FlexArgs & args) {
        const double fifth{args.back()};
        args.pop_back();
        const double fourth{args.back()};
        args.pop_back();
        const double third{args.back()};
        args.pop_back();
        const double second{args.back()};
        args.pop_back();
        const double first{args.back()};
        args.pop_back();
        using std::max;
        return max(max(max(max(first, second), third), fourth), fifth);
      };
      break;
    case FunctionType::min:
      n_args_ = 2;
      function = [](FlexArgs & args) {
        const double second{args.back()};
        args.pop_back();
        const double first{args.back()};
        args.pop_back();
        return std::min(first, second);
      };
      break;
    case FunctionType::min3:
      n_args_ = 3;
      function = [](FlexArgs & args) {
        const double third{args.back()};
        args.pop_back();
        const double second{args.back()};
        args.pop_back();
        const double first{args.back()};
        args.pop_back();
        return std::min(std::min(first, second), third);
      };
      break;
    case FunctionType::min4:
      n_args_ = 4;
      function = [](FlexArgs & args) {
        const double fourth{args.back()};
        args.pop_back();
        const double third{args.back()};
        args.pop_back();
        const double second{args.back()};
        args.pop_back();
        const double first{args.back()};
        args.pop_back();
        return std::min(std::min(std::min(first, second), third), fourth);
      };
      break;
    case FunctionType::min5:
      n_args_ = 5;
      function = [](FlexArgs & args) {
        const double fifth{args.back()};
        args.pop_back();
        const double fourth{args.back()};
        args.pop_back();
        const double third{args.back()};
        args.pop_back();
        const double second{args.back()};
        args.pop_back();
        const double first{args.back()};
        args.pop_back();
        using std::min;
        return min(min(min(min(first, second), third), fourth), fifth);
      };
      break;
    case FunctionType::sqrt:
      n_args_ = 1;
      function = [](FlexArgs & args) {
        const double first{args.back()};
        args.pop_back();
        return sqrt(first);
      };
      break;
    case FunctionType::none:
      throw Error("Bad none FunctionType in Function");
    default:
      throw Error("Unknown FunctionType in Function");
    }
  }

  double apply(FlexArgs & args) const {
    return function(args);
  }
  uint64_t n_args() const { return n_args_; }

 private:
  FlexFunction function{dummy};
  uint64_t n_args_{0};
};

// Holds a number or a function
class NumFun {
 public:
  explicit NumFun(Value && number__) :
      number_{std::move(number__)}, is_number_{true} {}
  explicit NumFun(Function && function__) :
      function_{std::move(function__)}, is_number_{false} {}

  bool is_number() const { return is_number_; }
  bool is_function() const { return !is_number_; }
  const Value & number() const { return number_; }
  const Function & function() const { return function_; }

 private:
  Value number_{};
  Function function_{};
  bool is_number_;
};

// A expression to be evaluated later,
// after quantities might change
// with symbols resolved to functions or numbers
class Expression {
 public:
  Expression(const std::string & expression__, Symbols & symbols) :
      Expression{expression__.begin(), expression__.end(), symbols} {}
  Expression(CI ci, const CI end, Symbols & symbols) :
      expression_{ci, end} {
    // std::cerr << "Expression" << std::endl;
    std::vector<Token> token_stack;
    uint64_t n_constant_at_top{0};
    auto op2ops = [this, &token_stack, &n_constant_at_top]() {
      Function function{token_stack.back().function_type()};
      token_stack.pop_back();
      if (n_constant_at_top >= function.n_args()) {
        std::vector<double> numbers;
        for (uint64_t a{0}; a != function.n_args(); ++a) {
          numbers.push_back(
              num_funs[num_funs.size() - function.n_args() + a].number());
        }
        for (uint64_t a{0}; a != function.n_args(); ++a) num_funs.pop_back();
        num_funs.emplace_back(Value{function.apply(numbers)});
        n_constant_at_top += 1;
        n_constant_at_top -= function.n_args();
      } else {
        num_funs.emplace_back(std::move(function));
        n_constant_at_top = 0;
      }
    };
    // Shunting yard algorithm, with simple constant reduction, infix -> rpn
    TokenType type{TokenType::end};
    uint64_t n_args{0};
    while (ci != end) {
      Token token{get_token(ci, end, type)};
      type = token.token_type();
      // std::cerr << "Token " << name(type) << std::endl;
      if (type == TokenType::end) break;
      switch (type) {
        case TokenType::end:
          throw Error("Unexpected end TokenType in Expression");
        case TokenType::open:
          n_args = 0;
          token_stack.push_back(token);
          break;
        case TokenType::close:
          while (token_stack.size() &&
                 token_stack.back().token_type() != TokenType::open) {
            op2ops();
            ++n_args;
          }
          if (token_stack.size() &&
              token_stack.back().token_type() == TokenType::open)
            token_stack.pop_back();
          break;
        case TokenType::comma:
          break;
        case TokenType::number:
          num_funs.emplace_back(Value{token.symbol()});
          if (0) std::cerr << "Number " << token.symbol() << " "
                           << num_funs.back().number() << std::endl;
          ++n_constant_at_top;
          break;
        case TokenType::name:
          num_funs.emplace_back(Value{symbols, token.symbol(), true});
          if (0) std::cerr << "Number " << token.symbol() << " "
                           << num_funs.back().number() << std::endl;
          n_constant_at_top = 0;
          break;
        case TokenType::Operator:
          // only works with left associative operators!
          while (token_stack.size() &&
                 token_stack.back().token_type() != TokenType::open &&
                 (token_stack.back().token_type() == TokenType::function ||
                  precedence(token_stack.back().function_type()) >=
                  precedence(token.function_type()))) op2ops();
          if (0) std::cerr << "Function "
                           << name(token.function_type()) << std::endl;
          token_stack.push_back(token);
          break;
        case TokenType::function:
          token_stack.push_back(token);
          break;
        default:
          throw Error("Unexpected case in Expression");
      }
    }
    while (token_stack.size()) op2ops();
    // std::cerr << "num_funs " << num_funs.size() << std::endl;
  }
  const std::string & expression() const { return expression_; }
  // Perform the operation specified by the compiled expression
  double apply() const {
    // std::cerr << "apply" << std::endl;
    std::vector<double> numbers;
    for (const NumFun & num_fun : num_funs) {
      if (num_fun.is_number()) {
        numbers.push_back(num_fun.number());
      } else {
        const Function & function{num_fun.function()};
        if (numbers.size() < function.n_args())
          throw Error("Too few numbers in operation") << numbers.size();
        numbers.push_back(function.apply(numbers));
      }
      // std::cerr << "Numbers " << numbers << std::endl;
    }
    if (numbers.size() != 1)
      throw Error("Number stack nonempty") << numbers.size();
    return numbers.back();
  }
  uint64_t size() const { return num_funs.size(); }

 private:
  std::vector<NumFun> num_funs{};
  std::string expression_{};
};

bool test_expression(const std::string & expression_string, const double answer,
                     Symbols & symbols) {
  const Expression expression{expression_string, symbols};
  const double result{expression.apply()};
  if (result < answer || result > answer)
    throw Error("result : answer : expression problem.")
        << result << ":" << answer << ":" << expression_string;
  std::cout << "GOOD result for expression of size "
            << expression.size() << ": "
            << expression_string << " = " << result << std::endl;
  return true;
}

}  // namespace paa

#endif  // PAA_UTILITY_LPE_H_
