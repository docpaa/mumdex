//
// gene_info.cpp
//
// Output gene information
//
// Copyright 2020 Peter Andrews @ CSHL
//

#include <algorithm>
#include <exception>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "error.h"
#include "genes.h"
#include "mumdex.h"

using std::cout;
using std::cerr;
using std::endl;
using std::exception;
using std::ostringstream;
using std::pair;
using std::string;
using std::vector;

using paa::ChromosomeIndexLookup;
using paa::Error;
using paa::GeneLookup;
using paa::GeneXrefs;
using paa::KnownGene;
using paa::KnownGenes;
using paa::Reference;

int main(int argc, char * argv[]) try {
  // Check initial command line arguments
  const string usage{"usage: gene_info reference gene ..."};
  if (--argc < 2) throw Error(usage);
  const string reference_file{argv[1]};
  const Reference reference{reference_file};
  const ChromosomeIndexLookup chr_lookup{reference};
  const KnownGenes genes{chr_lookup, reference};
  const GeneXrefs xref{reference};
  const GeneLookup gene_lookup{genes, xref};
  --argc;
  argv += 2;

  // Process all gene names
  for (int a{0}; a != argc; ++a) {
    const string gene{argv[a]};
    if (a) cout << "\n";
    cout << "Gene " << gene << "\n";
    vector<pair<double, string>> ordered;
    const GeneLookup::ER isoforms{gene_lookup.isoforms(gene)};
    for (GeneLookup::Iter isoform_iter{isoforms.first};
         isoform_iter != isoforms.second; ++isoform_iter) {
      ostringstream out;
      const KnownGene & isoform{*isoform_iter->second};
      uint64_t total_exon_length{0};
      for (uint64_t exon{0}; exon != isoform.n_exons; ++exon)
        total_exon_length += isoform.exon_stops[exon] -
            isoform.exon_starts[exon];
      out << isoform.name << ' ' << isoform.chr_name
           << ' ' << isoform.t_start + 1 << ' ' << isoform.t_stop + 1
           << ' ' << isoform.c_start + 1 << ' ' << isoform.c_stop + 1
           << ' ' << isoform.n_exons << ' ' << total_exon_length << "\n";
        for (uint64_t exon{0}; exon != isoform.n_exons; ++exon)
        out << exon << ' ' << isoform.exon_starts[exon] + 1
            << ' ' << isoform.exon_stops[exon] + 1 << ' '
            << reference.subseq(isoform.chr, isoform.exon_starts[exon],
                                isoform.exon_stops[exon]) << '\n';
        ordered.emplace_back(-1.0 * total_exon_length, out.str());
    }
    sort(ordered.begin(), ordered.end());
    for (auto info : ordered) cout << info.second;
  }

  return 0;
} catch (Error & e) {
  cerr << "paa::Error:" << endl;
  cerr << e.what() << endl;
  return 1;
} catch (exception & e) {
  cerr << "std::exception" << endl;
  cerr << e.what() << endl;
  return 1;
} catch (...) {
  cerr << "unknown exception was caught" << endl;
  return 1;
}
