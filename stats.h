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

#include "info.h"
#include "stats.h"
#include "variant.h"

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

struct indel {
  unsigned int insertions;
  unsigned int deletions;
};

struct snpTypes {
  unsigned int novelTransitions;
  unsigned int knownTransitions;
  unsigned int novelTransversions;
  unsigned int knownTransversions;
};

struct variantStruct {
  unsigned int novelTransitions;
  unsigned int knownTransitions;
  unsigned int diffKnownTransitions;
  unsigned int novelTransversions;
  unsigned int knownTransversions;
  unsigned int diffKnownTransversions;
  unsigned int multiAllelic;

  map<unsigned int, unsigned int> mnps;

  map<unsigned int, indel> indels;
  map<unsigned int, snpTypes> acs;
  map<double, snpTypes> afs;
  map<string, unsigned int> annotationsTs;
  map<string, unsigned int> annotationsTv;
  map<string, unsigned int> annotationsTriallelicSnp;
  map<string, unsigned int> annotationsQuadallelicSnp;
  map<string, unsigned int> annotationsMnp;
  map<string, unsigned int> annotationsIns;
  map<string, unsigned int> annotationsDel;

// Overload the + operator for structures.
  variantStruct operator+(variantStruct& vs) {
    variantStruct result;
    map<unsigned int, indel>::iterator indelIter;
    map<unsigned int, snpTypes>::iterator acsIter;
    map<double, snpTypes>::iterator afsIter;
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

    // Indels.
    for (indelIter = vs.indels.begin(); indelIter != vs.indels.end(); indelIter++) {
      result.indels[indelIter->first].insertions = this->indels[indelIter->first].insertions + vs.indels[indelIter->first].insertions;
      result.indels[indelIter->first].deletions = this->indels[indelIter->first].deletions + vs.indels[indelIter->first].deletions;
    }

    // Allele frequency spectrum.
    for (afsIter = vs.afs.begin(); afsIter != vs.afs.end(); afsIter++) {
      result.afs[afsIter->first].novelTransitions = this->afs[afsIter->first].novelTransitions + vs.afs[afsIter->first].novelTransitions;
      result.afs[afsIter->first].knownTransitions = this->afs[afsIter->first].knownTransitions + vs.afs[afsIter->first].knownTransitions;
      result.afs[afsIter->first].novelTransversions = this->afs[afsIter->first].novelTransversions + vs.afs[afsIter->first].novelTransversions;
      result.afs[afsIter->first].knownTransversions = this->afs[afsIter->first].knownTransversions + vs.afs[afsIter->first].knownTransversions;
    }

    // Allele count spectrum.
    for (acsIter = vs.acs.begin(); acsIter != vs.acs.end(); acsIter++) {
      result.acs[acsIter->first].novelTransitions = this->acs[acsIter->first].novelTransitions + vs.acs[acsIter->first].novelTransitions;
      result.acs[acsIter->first].knownTransitions = this->acs[acsIter->first].knownTransitions + vs.acs[acsIter->first].knownTransitions;
      result.acs[acsIter->first].novelTransversions = this->acs[acsIter->first].novelTransversions + vs.acs[acsIter->first].novelTransversions;
      result.acs[acsIter->first].knownTransversions = this->acs[acsIter->first].knownTransversions + vs.acs[acsIter->first].knownTransversions;
    }

    // Annotations (transitions).
    for (sIter = vs.annotationsTs.begin(); sIter != vs.annotationsTs.end(); sIter++) {
      result.annotationsTs[sIter->first] = this->annotationsTs[sIter->first] + vs.annotationsTs[sIter->first];
    }

    // Annotations (transversions).
    for (sIter = vs.annotationsTv.begin(); sIter != vs.annotationsTv.end(); sIter++) {
      result.annotationsTv[sIter->first] = this->annotationsTv[sIter->first] + vs.annotationsTv[sIter->first];
    }

    // Annotations (tri-allelic SNPs).
    for (sIter = vs.annotationsTriallelicSnp.begin(); sIter != vs.annotationsTriallelicSnp.end(); sIter++) {
      result.annotationsTriallelicSnp[sIter->first] = this->annotationsTriallelicSnp[sIter->first] + vs.annotationsTriallelicSnp[sIter->first];
    }

    // Annotations (quad-allelic SNPs).
    for (sIter = vs.annotationsQuadallelicSnp.begin(); sIter != vs.annotationsQuadallelicSnp.end(); sIter++) {
      result.annotationsQuadallelicSnp[sIter->first] = this->annotationsQuadallelicSnp[sIter->first] + vs.annotationsQuadallelicSnp[sIter->first];
    }

    // Annotations (MNPs).
    for (sIter = vs.annotationsMnp.begin(); sIter != vs.annotationsMnp.end(); sIter++) {
      result.annotationsMnp[sIter->first] = this->annotationsMnp[sIter->first] + vs.annotationsMnp[sIter->first];
    }

    // Annotations (insertions).
    for (sIter = vs.annotationsIns.begin(); sIter != vs.annotationsIns.end(); sIter++) {
      result.annotationsIns[sIter->first] = this->annotationsIns[sIter->first] + vs.annotationsIns[sIter->first];
    }

    // Annotations (transversions).
    for (sIter = vs.annotationsDel.begin(); sIter != vs.annotationsDel.end(); sIter++) {
      result.annotationsDel[sIter->first] = this->annotationsDel[sIter->first] + vs.annotationsDel[sIter->first];
    }

    return result;
  }
};

class statistics {
  public:
    statistics(void);
    ~statistics(void);
    void generateStatistics(variant&, vcf&, int, bool, vector<string>&, bool);
    void getAnnotations(vector<string>&, variantInfo&, map<string, unsigned int>&);
    void printSnpStatistics(ostream*);
    void printSnpAnnotations(ostream*);
    void printAcs(ostream*);
    void printAfs(ostream*);
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
    int lastSnpPosition;
    int lastMnpPosition;
    int lastIndelPosition;
    map<string, string> referenceSequences;
    map<string, map<string, variantStruct> > variants;
    map<string, map<string, variantStruct> > totalVariants;
    map<string, map<string, unsigned int> > distributions;
    map<int, unsigned int> snpDistribution;
    map<string, unsigned int> annotationNames;
};

} // namespace vcfCTools

#endif
