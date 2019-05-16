//
// ivan_metadata.h
//
// interface to access Ivan's metadata format
//
// Copyright 2019 Peter Andrews CSHL
//

#ifndef PAA_IVAN_METADATA_H
#define PAA_IVAN_METADATA_H

#include <iostream>
#include <map>
#include <string>

namespace paa {

struct IvanSampleInfo {
  struct eof {};
  explicit IvanSampleInfo(std::istream & in) {
    getline(in, atts, '\t');
    char line_end;
    getline(in, bam, '\t');
    getline(in, sample_id, '\t');
    getline(in, family_id, '\t');
    getline(in, ind_id, '\t');
    getline(in, mom, '\t');
    if (mom.empty()) mom = "0";
    getline(in, dad, '\t');
    if (dad.empty()) dad = "0";
    getline(in, sex, '\t');
    in >> role;
    line_end = in.get();
    if (!in) {
      if (in.eof()) throw eof();
      throw Error("IvanSampleInfo parse error");
    }
    if (role != "prb" && role != "sib" && role != "mom" && role != "dad")
      throw Error("IvanSampleInfo role parse error")
          << bam << sample_id << role;
    if (line_end != '\n') throw Error("IvanSampleInfo line end parse error");
  }

  bool is_proband() const { return role == "prb"; }
  bool is_male() const { return sex == "M"; }
  bool is_parent() const { return role == "mom" || role == "dad"; }

  std::string atts{};
  std::string bam{};
  std::string sample_id{};
  std::string family_id{};
  std::string ind_id{};
  std::string mom{};
  std::string dad{};
  std::string sex{};
  std::string role{};
};

class IvanMetadata {
 public:
  explicit IvanMetadata(const std::string & file_name) {
    std::ifstream file{file_name.c_str()};
    if (!file) {
      throw Error("Problem opening lookup file") << file_name;
    }
    file.ignore(100000, '\n');
    while (file) {
      try {
        IvanSampleInfo info{file};
        const auto found = lookup.emplace(info.ind_id, info);
        if (found.second == false)
          throw Error("Duplicate entry for individual") << info.ind_id;
      } catch (const IvanSampleInfo::eof &) {
        break;
      }
    }
  }
  const IvanSampleInfo & operator[](const std::string & individual) const {
    try {
      return lookup.at(individual);
    } catch (...) {
      std::cerr << "Lookup problem for " << individual << std::endl;
      throw;
    }
  }
  std::string role(const std::string & individual) const {
    try {
      return lookup.at(individual).role;
    } catch (...) {
      std::cerr << "Role lookup problem for " << individual << std::endl;
      throw;
    }
  }
  std::string sex(const std::string & individual) const {
    try {
      return lookup.at(individual).sex;
    } catch (...) {
      std::cerr << "Sex lookup problem for " << individual << std::endl;
      throw;
    }
  }
  std::string mom(const std::string & individual) const {
    try {
      return lookup.at(individual).mom;
    } catch (...) {
      std::cerr << "Mom lookup problem for " << individual << std::endl;
      throw;
    }
  }
  std::string dad(const std::string & individual) const {
    try {
      return lookup.at(individual).dad;
    } catch (...) {
      std::cerr << "Dad lookup problem for " << individual << std::endl;
      throw;
    }
  }
  std::string family(const std::string & individual) const {
    try {
      return lookup.at(individual).family_id;
    } catch (...) {
      std::cerr << "Family lookup problem for " << individual << std::endl;
      throw;
    }
  }
  std::string sample(const std::string & individual) const {
    try {
      return lookup.at(individual).sample_id;
    } catch (...) {
      std::cerr << "Sample lookup problem for " << individual << std::endl;
      throw;
    }
  }
  std::string bam(const std::string & individual) const {
    try {
      return lookup.at(individual).bam;
    } catch (...) {
      std::cerr << "Bam lookup problem for " << individual << std::endl;
      throw;
    }
  }
  std::string first_family() const {
    if (lookup.empty()) throw Error("Lookup on empty map in IvanMetadata");
    return lookup.begin()->second.family_id;
  }

 private:
  std::map<std::string, IvanSampleInfo> lookup{};
};

}  // namespace paa


#endif  // PAA_IVAN_METADATA_H

