//
// sssa.h
//
// species-specific suffix array
//
// Copyright 2016 Peter Andrews @ CSHL
//

#ifndef PAA_SSSA_H
#define PAA_SSSA_H

#include <fstream>
#include <limits>
#include <map>
#include <string>

#include "error.h"
#include "files.h"
#include "longSA.h"
#include "utility.h"

namespace paa {

constexpr unsigned int high_bit(const uint64_t value) {
  return value ? high_bit(value >> 1) + 1 : 0;
}

template <unsigned int n_bits>
struct BitInfo {
  static_assert(n_bits, "No bits in BitInfo");
  using WordType = typename BitInfo<n_bits - 1>::WordType;
};

template <> struct BitInfo<1> { using WordType = uint8_t; };
template <> struct BitInfo<9> { using WordType = uint16_t; };
template <> struct BitInfo<17> { using WordType = uint32_t; };
template <> struct BitInfo<33> { using WordType = uint64_t; };
template <> struct BitInfo<64> { using WordType = uint64_t; };
template <> struct BitInfo<128> { using WordType = uint64_t; };
template <> struct BitInfo<256> { using WordType = uint64_t; };
template <> struct BitInfo<512> { using WordType = uint64_t; };
template <> struct BitInfo<1024> { using WordType = uint64_t; };
template <> struct BitInfo<2048> { using WordType = uint64_t; };

template <unsigned int species_bits>
class Species_t {
 public:
  static constexpr unsigned int n_bits() { return species_bits; }
  static constexpr unsigned int max_n() { return (1ul << n_bits()) - 1; }
  static constexpr unsigned int n_values() { return max_n() + 1; }
  using SPECIES = typename BitInfo<n_bits()>::WordType;
  static constexpr SPECIES multiple() { return max_n(); }
  static constexpr SPECIES mask() { return max_n(); }

  explicit Species_t(const uint64_t species__) :
      species_{static_cast<SPECIES>(species__)} { }
  SPECIES species() const { return species_; }
  SPECIES operator()() const { return species_; }
  void combine(const Species_t species__) {
    if (species_ != species__.species_) {
      species_ = multiple();
    }
  }
  bool operator==(const SPECIES rhs) const { return species_ == rhs; }
  bool operator<(const Species_t rhs) const { return species_ < rhs.species_; }

 private:
  SPECIES species_{0};

  static_assert(n_bits() <= 32, "Too many bits in Species");
  static_assert(max_n(), "max_n() == 0 in Species");
};

template <unsigned int n_bases>
class Kmer_t {
 public:
  static_assert(n_bases, "No bases in KmerInfo");
  using WordType = typename BitInfo<n_bases * 2>::WordType;
  static constexpr unsigned int word_size() {
    return sizeof(WordType);
  }
  static constexpr unsigned int word_bits() {
    return 8 * word_size();
  }
  static constexpr unsigned int word_bases() {
    return 4 * word_size();
  }
  static constexpr unsigned int n_words() {
    return n_bases / word_bases() +
        static_cast<bool>(n_bases % word_bases());
  }
  static constexpr unsigned int N_WORDS{n_words()};
  static constexpr unsigned int bits_available() {
    return n_words() * word_bits();
  }
  static constexpr unsigned int bits_used() {
    return n_bases * 2;
  }
  static constexpr unsigned int extra_bits() {
    return bits_available() - bits_used();
  }
  static constexpr unsigned int first_word_bases() {
    return word_bases() - 2 * extra_bits();
  }
  static constexpr double n_possible() {
    return pow(4, n_bases);
  }
  explicit Kmer_t(const std::string & sequence_) {
    if (sequence_.size() < n_bases) {
      throw Error("Kmer constructed with too small a sequence");
    }
    unsigned int b{0};
    data[n_words() - 1] = 0;
    for (unsigned int w{0}; w != n_words(); ++w) {
      for (unsigned int wb{0}; b != n_bases && wb != word_bases(); ++wb) {
        (data[w] <<= 2) |= base2int(sequence_[b++]);
      }
    }
  }
  explicit Kmer_t(const char * b, const char * const e) {
    if (e - b < n_bases) {
      throw Error("Kmer constructed with too small a sequence");
    }
    data[n_words() - 1] = 0;
    for (unsigned int w{0}; w != n_words(); ++w) {
      for (unsigned int wb{0}; b != e && wb != word_bases(); ++wb) {
        (data[w] <<= 2) |= base2int(*(b++));
      }
    }
  }
  std::string sequence() const {
    std::string result;
    unsigned int b{0};
    for (unsigned int w{0}; w != n_words(); ++w) {
      for (unsigned int wb{0}; b != n_bases && wb != word_bases(); ++wb) {
        result += int2base(data[w] >> 2 * (word_bases() - wb - 1) & 3);
      }
    }
    return result;
  }
  std::string operator()() const { return sequence(); }

  struct bad_base { };
 private:
  static WordType base2int(const char c) {
    switch (c) {
      case 'A': return 0;
      case 'C': return 1;
      case 'G': return 2;
      case 'T': return 3;
      default:
        throw bad_base{};
    }
  }
  static constexpr char int2base(const WordType d) {
    return (d == 0) ? 'A' : (d == 1 ? 'C': (d == 2 ? 'G' : 'T'));
  }
#if 0
  KMER operator>>(const unsigned int shift) const { return data >> shift; }
  KMER operator<<(const unsigned int shift) const { return data << shift; }
  KMER operator&(const KMER mask) const { return data & mask; }
  KMER operator|(const KMER species) const { return data & species; }
#endif

  WordType data[N_WORDS];
};

#if 0
template <unsigned int n_bases, unsigned int max_n_species>
class KmerSpecies {
 public:
  using Species = Species_t<max_n_species>;
  using KMER = uint64_t;


 public:
  KmerSpecies(const std::string & sequence, const Species & species_) :
      data{(Kmer{sequence} << kmer_shift()) | (species+ & species_mask())} {
    static_assert(sizeof(Kmer) > sizeof(Species),
                  "KMER size is not bigger than SPECIES size");
  }

  Kmer kmer() const { return data >> kmer_shift(); }
  Species species() const { return Species(data & species_mask()); }

 private:
  static constexpr unsigned int n_bits() { return sizeof(KMER) * 8; }
  static constexpr unsigned int kmer_bits() {
    return n_bits() - Species::n_bits();
  }
  static constexpr unsigned int max_kmer_size() { return kmer_bits / 2; }
  static constexpr Kmer species_mask() { return Species::mask(); }
  static constexpr unsigned int kmer_shift() { return Species::n_bits(); }
  static_assert(n_bases * 2 + kmer_shift() <= n_bits(),
                "Too many bases and species for Kmer size");

  Kmer data;
};
#endif

template <unsigned int n_bases, unsigned int max_n_species>
class SSSA_t {
 public:
  using Kmer = Kmer_t<n_bases>;
  static constexpr unsigned int species_bits() {
    return high_bit(max_n_species);
  }
  using Species = Species_t<species_bits()>;

  // Empty SSSA
  SSSA_t() { }

  // Load from fasta or find SSSA dir and load from there
  explicit SSSA_t(const std::string & fasta,
                  std::string name__ = "") {
    if (!readable(fasta)) throw Error("fasta file not found") << fasta;

    // See if needs to be loaded or created
    const std::string sssa_dir{[&fasta] {
        std::ostringstream result;
        result << fasta << ".sssa.kmer." << n_bases
               << ".species." << species_bits() << ".bin/";
        return result.str();
      }()};

    if (readable(sssa_dir)) {
    } else {
      // Load or create suffix array for species
      const longSA sa{fasta, true, true, true};
      if (n_bases > sa.N) throw Error("kmer size greater than sequence size");


      // Write species name to file
      mkdir(sssa_dir);
      std::ofstream names{sssa_dir + "names.txt"};
      if (!names) throw Error("Could not open names file for writing");
      names << (name__.size() ? name__ : fasta) << " " << 0 << "\n";

      // Transform SA into SSSA
      BinWriter kmers_file{sssa_dir + "kmers.bin"};
      const uint64_t n{sa.N};
      Progress progress{n, 0.1, "SSSA creation"};
      for (uint64_t i{0}; i != n; ++i) {
        progress();
        if (sa.LCP[i] >= n_bases) continue;
        const uint64_t r{sa.SA[i]};
        if (r + n_bases > n) continue;
        try {
          kmers_file << Kmer{sa.ref.seq + r, sa.ref.seq + r + n_bases};
        } catch (typename Kmer::bad_base &) {
          continue;
        }
      }

      // Empty species file for one species SSSA
      BinWriter species_file{sssa_dir + "species.bin"};
    }
    new (this) SSSA_t(sssa_dir, true);
  }

  // Load existing SSSA dir
  SSSA_t(const std::string & sssa_dir, bool) :
      kmers_{sssa_dir + "kmers.bin"},
    species_{sssa_dir + "species.bin", false} {
      std::ifstream names{sssa_dir + "names.txt"};
      if (!names) throw Error("Could not open names file for reading");
      std::string name__;
      unsigned int species__;
      while (names >> name__ >> species__) {
        add(name__, species__);
      }
    }

  void merge(const SSSA_t &) {
  }
  uint64_t size() const { return kmers_.size(); }
  unsigned int n_species() const { return name_species_.size(); }
  Species species(const uint64_t i) const { return species_[i]; }
  Kmer kmer(const uint64_t i) const { return kmers_[i]; }
  static std::string info() {
    std::ostringstream out;
    out << "kmer size " << n_bases << " "
        << "max_n_species " << max_n_species << " "
        << "kmer bytes " << sizeof(Kmer) << " "
        << "kmer bits " << Kmer::bits_used() << " "
        << "n kmers " << Kmer::n_possible() << " "
        << "species bits " << Species::n_bits() << " "
        << "extra bits " << Kmer::extra_bits() - Species::n_bits();
    return out.str();
  }
  double kmer_density() const { return size() / Kmer::n_possible(); }

 private:
  void add(const std::string & name__, const unsigned int species__) {
    const Species species___{species__};
    name_species_.emplace(name__, species___);
    species_names_.emplace(species___, name__);
  }

  std::map<std::string, Species> name_species_{};
  std::map<Species, std::string> species_names_{};
  UnMappedVector<Kmer> kmers_{};
  UnMappedVector<Species> species_{};
};


template <unsigned int n_bases>
using SSSA = SSSA_t<n_bases, Species_t<Kmer_t<n_bases>::extra_bits()>::max_n()>;

}  // namespace paa

#endif  // PAA_SSSA_H

