// ******************************************************
// vcfCTools (c) 2011 Alistair Ward
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ------------------------------------------------------
// Last modified: 18 February 2011
// ------------------------------------------------------
// Generate the statistics class.
// ******************************************************

#ifndef STATS_H
#define STATS_H

#include "stats.h"
#include "vcf.h"

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <string>
#include <getopt.h>
#include <stdlib.h>
#include <map>
#include <vector>
#include <algorithm>

using namespace std;

namespace vcfCTools {

struct variantStruct {
  unsigned int novelTransitions;
  unsigned int knownTransitions;
  unsigned int novelTransversions;
  unsigned int knownTransversions;
  unsigned int multiAllelic;

  unsigned int hapmap;

  map<unsigned int, unsigned int> insertions;
  map<unsigned int, unsigned int> deletions;

// Overload the + operator for structures
  variantStruct operator+(variantStruct& vs) {
    variantStruct result;
    result.novelTransitions   = this->novelTransitions   + vs.novelTransitions;
    result.knownTransitions   = this->knownTransitions   + vs.knownTransitions;
    result.novelTransversions = this->novelTransversions + vs.novelTransversions;
    result.knownTransversions = this->knownTransversions + vs.knownTransversions;
    result.multiAllelic       = this->multiAllelic       + vs.multiAllelic;

    result.hapmap             = this->hapmap + vs.hapmap;

    result.insertions = this->insertions;
    result.deletions  = this->deletions;
    for (map<unsigned int, unsigned int>::iterator iter = vs.insertions.begin(); iter != vs.insertions.end(); iter++) {
      result.insertions[(*iter).first] = result.insertions[(*iter).first] + vs.insertions[(*iter).first];
    }
    for (map<unsigned int, unsigned int>::iterator iter = vs.deletions.begin(); iter != vs.deletions.end(); iter++) {
      result.deletions[(*iter).first] = result.deletions[(*iter).first] + vs.deletions[(*iter).first];
    }

    return result;
  }
};

class statistics {
  public:
    statistics(void);
    ~statistics(void);
    void generateStatistics(vcf&);
    void printStatistics(ostream*);
    void countByFilter();
    void printVariantStruct(ostream*, string&, variantStruct&);

  public:
    bool isTransition;
    bool isTransversion;
    bool inDbsnp;
    bool inHapmap;
    string currentReferenceSequence;
    unsigned int lastSnpPosition;
    map<string, string> referenceSequences;
    map<string, map<string, variantStruct> > variants;
    map<string, map<string, variantStruct> > totalVariants;
    map<string, map<string, unsigned int> > distributions;
    map<unsigned int, unsigned int> snpDistribution;
};

} // namespace vcfCTools

#endif
