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

#include "header.h"
#include "info.h"
#include "output.h"
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

// Define structures for use with the statistics.

// Define a structure to hold counts of all the different variant types
// and annotations etc.  This also includes overloaded '+' operation.
struct variantStruct {
  unsigned int novelTransitions;
  unsigned int knownTransitions;
  unsigned int diffKnownTransitions;
  unsigned int novelTransversions;
  unsigned int knownTransversions;
  unsigned int diffKnownTransversions;
  unsigned int multiAllelic;
  unsigned int novelAminations;
  unsigned int novelDeaminations;
  unsigned int knownAminations;
  unsigned int knownDeaminations;

  map<unsigned int, unsigned int> mnps;

  map<unsigned int, indel> indels;
  //map<unsigned int, snpTypes> acs;
  //map<double, snpTypes> afs;
  map<string, unsigned int> annotationsTs;
  map<string, unsigned int> annotationsTv;
  map<string, unsigned int> annotationsTriallelicSnp;
  map<string, unsigned int> annotationsQuadallelicSnp;
  map<string, unsigned int> annotationsMnp;
  map<string, unsigned int> annotationsIns;
  map<string, unsigned int> annotationsDel;

// Overload the + operator for structures.
  //variantStruct operator+(variantStruct& vs) {
  variantStruct operator+(variantStruct& vs) {
    variantStruct result = *this;
    map<unsigned int, indel>::iterator indelIter;
    //map<unsigned int, snpTypes>::iterator acsIter;
    //map<double, snpTypes>::iterator afsIter;
    map<unsigned int, unsigned int>::iterator iter;
    map<string, unsigned int>::iterator sIter;

    result.novelTransitions       = this->novelTransitions       + vs.novelTransitions;
    result.novelTransversions     = this->novelTransversions     + vs.novelTransversions;
    result.novelAminations        = this->novelAminations        + vs.novelAminations;
    result.novelDeaminations      = this->novelDeaminations      + vs.novelDeaminations;
    result.knownTransitions       = this->knownTransitions       + vs.knownTransitions;
    result.diffKnownTransitions   = this->diffKnownTransitions   + vs.diffKnownTransitions;
    result.knownTransversions     = this->knownTransversions     + vs.knownTransversions;
    result.diffKnownTransversions = this->diffKnownTransversions + vs.diffKnownTransversions;
    result.knownAminations        = this->knownAminations        + vs.knownAminations;
    result.knownDeaminations      = this->knownDeaminations      + vs.knownAminations;
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
//    for (afsIter = vs.afs.begin(); afsIter != vs.afs.end(); afsIter++) {
//      result.afs[afsIter->first].novelTransitions = this->afs[afsIter->first].novelTransitions + vs.afs[afsIter->first].novelTransitions;
//      result.afs[afsIter->first].knownTransitions = this->afs[afsIter->first].knownTransitions + vs.afs[afsIter->first].knownTransitions;
//      result.afs[afsIter->first].novelTransversions = this->afs[afsIter->first].novelTransversions + vs.afs[afsIter->first].novelTransversions;
//      result.afs[afsIter->first].knownTransversions = this->afs[afsIter->first].knownTransversions + vs.afs[afsIter->first].knownTransversions;
//    }

    // Allele count spectrum.
//    for (acsIter = vs.acs.begin(); acsIter != vs.acs.end(); acsIter++) {
//      result.acs[acsIter->first].novelTransitions = this->acs[acsIter->first].novelTransitions + vs.acs[acsIter->first].novelTransitions;
//      result.acs[acsIter->first].knownTransitions = this->acs[acsIter->first].knownTransitions + vs.acs[acsIter->first].knownTransitions;
//      result.acs[acsIter->first].novelTransversions = this->acs[acsIter->first].novelTransversions + vs.acs[acsIter->first].novelTransversions;
//      result.acs[acsIter->first].knownTransversions = this->acs[acsIter->first].knownTransversions + vs.acs[acsIter->first].knownTransversions;
//    }

    // Annotations (transitions).
    for (sIter = vs.annotationsTs.begin(); sIter != vs.annotationsTs.end(); sIter++) {
      result.annotationsTs[sIter->first] = this->annotationsTs[sIter->first] + sIter->second;
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

// Define a structure to hold statistics for individual samples.
struct sampleStats {

  // SNPs.
  unsigned int knownAminations;
  unsigned int knownDeaminations;
  unsigned int knownHetTransitions;
  unsigned int knownHetTransversions;
  unsigned int knownHomTransitions;
  unsigned int knownHomTransversions;
  unsigned int novelAminations;
  unsigned int novelDeaminations;
  unsigned int novelHetTransitions;
  unsigned int novelHetTransversions;
  unsigned int novelHomTransitions;
  unsigned int novelHomTransversions;

  // MNPs.
  unsigned int hetMnps;

  // Indels.
  unsigned int hetDeletions;
  unsigned int hetInsertions;
  unsigned int homDeletions;
  unsigned int homInsertions;

  // Complex.
  unsigned int hetComplex;
  unsigned int homComplex;
  unsigned int homMnps;
  unsigned int homRef;

  // Other.
  unsigned unknown;
};

// Define a structure for holding boolean flags about a particular
// variant.
struct statsFlags {
  bool het;
  bool inDbsnp;
  bool isAmination;
  bool isDeamination;
  bool isTransition;
  bool isTransversion;
};

// Define the statistics class.
class statistics {
  public:
    statistics(void);
    ~statistics(void);
    void countByFilter();
    void determineSnpType(variant&, string&, unsigned int);
    void generateStatistics(vcfHeader&, variant&, bool, vector<string>&, bool, output&);
    void getAnnotations(vector<string>&, variantInfo&, map<string, unsigned int>&);
    void parseGenotypes(vcfHeader&, variant&, vector<unsigned int>);
    void printAcs(output&);
    void printAfs(output&);
    void printDetailedHeader(output&);
    void printHeader(output&, string, bool, bool, bool);
    void printIndelStatistics(output&);
    void printMnpFilter(string&, output&);
    void printMnpStatistics(output&);
    void printSampleSnps(vcfHeader&, vcf&, output&);
    void printSnpAnnotations(output&);
    void printSnpAnnotationStruct(output&, string&, variantStruct&, string&);
    void printSnpStatistics(output&);
    void printVariantStruct(output&, string, variantStruct&);
    void updateSampleLevelStats(statsFlags&, unsigned int, string&);
    void updateDetailedSnps(variant&, vcf&, unsigned int, output&);

  public:
    bool isAmination;
    bool isDeamination;
    bool isTransition;
    bool isTransversion;
    bool inDbsnp;
    bool inHapmap;
    bool hasSnp;
    bool hasMultiSnp;
    bool hasMnp;
    bool hasInsertion;
    bool hasDeletion;
    bool hasAnnotations;
    bool splitMnps;
    string currentReferenceSequence;
    int deletionSize;
    int insertionSize;
    int lastSnpPosition;
    int lastMnpPosition;
    int lastIndelPosition;
    map<string, string> referenceSequences;
    map<string, map<string, variantStruct> > variants;
    map<string, map<string, variantStruct> > totalVariants;
    map<string, map<string, unsigned int> > distributions;
    map<int, unsigned int> snpDistribution;
    map<string, unsigned int> annotationNames;

    // Sample level statistics.
    bool generateDetailed;
    bool generateSampleStats;
    double minDetailedGenotypeQuality;
    double minGenotypeQuality;
    map<string, sampleStats> sampleLevelStats;
};

} // namespace vcfCTools

#endif
