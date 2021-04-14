//
// sha1.h
//
// useful sha1 functions
//
// Copyright Peter Andrews 2021 @ CSHL
//

#ifndef UTILITY_SHA1_H_
#define UTILITY_SHA1_H_

#include <openssl/sha.h>
#include <string>

#include "private.h"

namespace paa {

// Use sha1 hash for (slightly) better password security
// Client sends only hashed password, then we store only a salted hash of that
std::string sha1(const std::string & data) {
  unsigned char obuf[20];
  SHA1(reinterpret_cast<const unsigned char *>(&data[0]), data.size(), obuf);
  int i;
  std::string hash(40, 'x');
  for (i = 0; i < 20; i++) sprintf(&hash[2 * i], "%02x", obuf[i]);
  return hash;
}
std::string salted_sha1(const std::string & input,
                        const std::string & user_salt = "") {
  return sha1(password_salt + input + user_salt);
}

}  // namespace paa

#endif  // UTILITY_SHA1_H_
