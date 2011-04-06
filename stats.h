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
  unsigned int diffKnownTransitions;
  unsigned int novelTransversions;
  unsigned int knownTransversions;
  unsigned int diffKnownTransversions;
  unsigned int multiAllelic;

  map<unsigned int, unsigned int> mnps;

  map<unsigned int, unsigned int> insertions;
  map<unsigned int, unsigned int> deletions;
  map<string, unsigned int> annotationsTs;
  map<string, unsigned int> annotationsTv;
  map<string, unsigned int> annotationsIns;
  map<string, unsigned int> annotationsDel;

// Overload the + operator for structures.
  variantStruct operator+(variantStruct& vs) {
    variantStruct result;
    map<unsigned int, unsigned int>::iterator iter;
    map<string, unsigned int>::iterator sIter;

    result.novelTransitions       = this->novelTransitions       + vs.novelTransitions;
    result.knownTransitions       = this->knownTransitions       + vs.knownTransitions;
    result.diffKnownTransitions   = this->diffKnownTransitions   + vs.diffKnownTransitions;
    result.novelTransversions     = this->novelTransversions     + vs.novelTransversions;
    result.knownTransversions     = this->knownTransversions     + vs.knownTransversions;
    result.diffKnownTransversions = this->diffKnownTransversions + vs.diffKnownTransversions;
    result.multiAllelic           = this->multiAllelic           + vs.multiAllelic;

    // MNPs
    for (iter = vs.mnps.begin(); iter != vs.mnps.end(); iter++) {
      result.mnps[iter->first] = this->mnps[iter->first] + vs.mnps[iter->first];
    }

    // Insertions.
    for (iter = vs.insertions.begin(); iter != vs.insertions.end(); iter++) {
      result.insertions[iter->first] = this->insertions[iter->first] + vs.insertions[iter->first];
    }

    // Deletions.
    for (iter = vs.deletions.begin(); iter != vs.deletions.end(); iter++) {
      result.deletions[iter->first] = this->deletions[iter->first] + vs.deletions[iter->first];
    }

    // Annotations (transitions).
    for (sIter = vs.annotationsTs.begin(); sIter != vs.annotationsTs.end(); sIter++) {
      result.annotationsTs[sIter->first] = this->annotationsTs[sIter->first] + vs.annotationsTs[sIter->first];
    }

    // Annotations (insertions).
    for (sIter = vs.annotationsIns.begin(); sIter != vs.annotationsIns.end(); sIter++) {
      result.annotationsIns[sIter->first] = this->annotationsIns[sIter->first] + vs.annotationsIns[sIter->first];
    }

    // Annotations (transversions).
    for (sIter = vs.annotationsDel.begin(); sIter != vs.annotationsDel.end(); sIter++) {
      result.annotationsDel[sIter->first] = this->annotationsDel[sIter->first] + vs.annotationsDel[sIter->first];
    }

    // Annotations (transversions).
    for (sIter = vs.annotationsTv.begin(); sIter != vs.annotationsTv.end(); sIter++) {
      result.annotationsTv[sIter->first] = this->annotationsTv[sIter->first] + vs.annotationsTv[sIter->first];
    }

    return result;
  }
};

class statistics {
  public:
    statistics(void);
    ~statistics(void);
    void generateStatistics(vcf&);
    void printSnpStatistics(ostream*);
    void printSnpAnnotations(ostream*);
    void printMnpStatistics(ostream*);
    void printIndelStatistics(ostream*);
    void countByFilter();
    void printVariantStruct(ostream*, string&, variantStruct&);
    void printSnpAnnotationStruct(ostream*, string&, variantStruct&, string&);
    void printMnpFilter(string&, ostream*);

  public:
    bool isTransition;
    bool isTransversion;
    bool inDbsnp;
    bool inHapmap;
    bool hasSnp;
    bool hasMnp;
    bool hasIndel;
    bool hasAnnotations;
    string currentReferenceSequence;
    unsigned int lastSnpPosition;
    map<string, string> referenceSequences;
    map<string, map<string, variantStruct> > variants;
    map<string, map<string, variantStruct> > totalVariants;
    map<string, map<string, unsigned int> > distributions;
    map<unsigned int, unsigned int> snpDistribution;
    map<string, unsigned int> annotationNames;
};

} // namespace vcfCTools

#endif
