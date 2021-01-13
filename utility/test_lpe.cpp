//
// test_lpe.h.cpp
//
// test the lpe.h lexer, parser, evaluator code
//
// Copyright 2020 Peter Andrews @ CSHL
//

#include <iostream>
#include <string>
#include <vector>

#include "error.h"
#include "lpe.h"

using std::cerr;
using std::cin;
using std::cout;
using std::endl;
using std::exception;
using std::string;
using std::vector;

using paa::test_expression;
using paa::test_tokenization;
using paa::Error;
using paa::Symbols;


int main(int argc, char ** argv) try {
  const string usage{"test_lpe"};
  if (--argc != 0) throw Error(usage);
  if (0) cerr << argv[0];

  // NOTE: spaces in expressions here are deliberate as a test of success
  // the actual tokenizer does not require spaces to work
  // test_tokenization tests tokenization 1lpwsuccess for all possible spacing

  // These should succeed
  test_tokenization("abc + -325 / sin ( 3.2e-2 - z )");
  test_tokenization("f ( 2 - 3.7 , x * 2 )");
  test_tokenization("a % 4");
  test_tokenization("a + +3.2");

  // These should fail successfully
  test_tokenization("f ( 2 - 3.7 , x );", false);
  test_tokenization("-.2", false);
  test_tokenization("-3.2.3", false);
  test_tokenization("+", false);
  test_tokenization("bad$varname", false);
  test_tokenization("", false);

  // These should fail, but do not fail
  test_tokenization("* /", false);  // must catch problem later
  test_tokenization("-3.2 moo", false);  // must catch problem later

  // These should produce the required answers
  Symbols symbols{{"a", 1}, {"b", -3}, {"abc", 25}};
  test_expression("3+4-a", 6, symbols);
  test_expression("a*b", -3, symbols);
  test_expression("abc/a+a", 26, symbols);
  test_expression("a+sqrt(a-b)+a", 4, symbols);
  test_expression("pow(a,abc)", 1, symbols);
  test_expression("max3(2, 3, 1)", 3, symbols);

  return 0;
} catch (Error & e) {
  cerr << "paa::Error:" << endl;
  cerr << e.what() << endl;
  return 1;
} catch(exception & e) {
  cerr << e.what() << endl;
  return 1;
} catch(...) {
  cerr << "Some exception was caught." << endl;
  return 1;
}
