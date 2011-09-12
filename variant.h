// ******************************************************
// vcfCTools (c) 2011 Alistair Ward
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ------------------------------------------------------
// Last modified: 18 February 2011
// ------------------------------------------------------
// Define the variant class.
// ******************************************************

#ifndef VARIANT_H
#define VARIANT_H

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <map>
#include <vector>

#include "bed.h"
#include "bedStructure.h"
#include "genotype_info.h"
#include "info.h"
#include "output.h"
#include "tools.h"
#include "trim_alleles.h"
#include "vcf.h"

using namespace std;

namespace vcfCTools {

// Create astructure to hold all of the flags required to determine the intersection
// operations to be performed.
struct intFlags {
  bool annotate;
  bool findCommon;
  bool findUnion;
  bool findUnique;
  bool writeFromFirst;
  bool sitesOnly;
  bool whollyWithin;
};

// Define a structure containing the different variant
// types.
struct variantType {
  bool isBiallelicSnp;
  bool isTriallelicSnp;
  bool isQuadallelicSnp;
  bool isMnp;
  bool isInsertion;
  bool isDeletion;
};

// Define a structure that contains information about a
// particular locus.  This structure is used for variants
// still in their original form and at the position they
// appear in the vcf file.  It is possible that these
// variants can be reduced (e.g. CC -> CG could be rewritten
// as a C -> G SNP at the next position), but these are kept
// in a separate structure for use in intersections etc.
struct originalVariants {

  // General information.
  bool hasMultipleRecords;
  unsigned int numberOfRecordsAtLocus;
  int position;
  int maxPosition;
  double quality;
  string info;
  string filters;
  string referenceSequence;
  vector<int> reducedPosition;
 
  // Ref and alt allele information.
  unsigned int numberAlts;
  string ref;
  string rsid;
  vector<bool> filtered;
  vector<string> reducedRef;
  vector<string> alts;
  vector<string> reducedAlts;
  vector<variantType> type;

  // Genotype information.
  bool hasGenotypes;
  string genotypeFormat;
  string genotypes;
};

// Define a structure that stores the information about the
// modified variants.  For example CG -> CA at position 10 will
// be stored as this in the originalVariants structure, but will
// be stored as a G -> A SNP at position 11 in this structure.
// This structure will point back to the originalVariants
// structure for genotypes and info etc.
struct reducedVariants {

  // If there are multiple records at this locus, the recordNumber
  // indicates which record this variant comes from.
  unsigned int recordNumber;

  int originalPosition;
  int altID; // If alts were A,G, this will indicate 1 for A and 2 for G etc.
  string ref;
  string alt;
};

// At each locus, many different variants can exist.  The
// variantsAtLocus structure holds all of the variants
// present at this locus.
struct variantsAtLocus {
  string referenceSequence;
  vector<reducedVariants> snps;
  vector<reducedVariants> mnps;
  vector<reducedVariants> indels;
};

// Define a structure that contains information about the different reference
// sequences encountered.  This includes information such as the number of
// records parsed for each reference sequence and if an intersection is
// performed, information about whether the reference sequences have been
// used in comparisons.
struct refSeqInfo {
  unsigned int numberRecords;
  bool contiguous;
  bool usedInComparison;
};

class variant {
  public:
    variant(void);
    ~variant(void);
    void addVariantToStructure(int, variantDescription&);
    void annotateRecordBed(bedRecord&);
    void annotateRecordVcf(bool, int, unsigned int, string&, bool, bool);
    void buildOutputRecord(output&);
    bool buildVariantStructure(vcf&);
    void clearOriginalVariants(intFlags&, output&, bool);
    void clearReferenceSequence(vcf&, intFlags, string, output&, bool);
    void clearReferenceSequenceBed(vcf&, intFlags, string, output&);
    void clearType(variantType&);
    void compareVariantsSameLocus(variant&, intFlags);
    void compareAlleles(vector<reducedVariants>&, vector<reducedVariants>&, intFlags, variant&);
    void determineVariantsToProcess(bool, bool, bool, bool, bool, bool);
    void determineVariantType(string, int, string, string, variantType&, int, originalVariants&);
    vector<string> extractGenotypeField(string);
    void filterUnique();
    void updateVariantMaps(string, variantType, string, string, int, string, originalVariants&);

  public:
    unsigned int recordsInMemory;
    string referenceSequence;

    // Structure containing variant information at a particular locus
    // after the variants have been deconstructed.  For example a variant
    // included in the vcf file as CC -> TT,CG at position x, is broken
    // into CC -> TT at position x and C -> G at x+1.  The original alleles
    // and positions are stored in the originalVariantsMap.
    map<int, variantsAtLocus> variantMap;
    map<int, variantsAtLocus>::iterator vmIter;

    // Structure containing variant information in the order that it
    // appeared in the vcf file.
    map<int, vector<originalVariants> > originalVariantsMap;
    map<int, vector<originalVariants> >::iterator ovmIter;
    vector<originalVariants>::iterator ovIter;

    // Header information.
    map<string, headerInfoStruct> headerInfoFields;
    map<string, headerInfoStruct> headerFormatFields;

    // Information about encountered reference sequences.
    map<string, refSeqInfo> referenceSequenceInfo;
    map<string, refSeqInfo>::iterator refSeqIter;

    // Boolean flags.
    bool assessAlts;
    bool isDbsnp;
    bool locusHasQuadSnp;
    bool locusHasSnp;
    bool locusHasTriSnp;
    bool processSnps;
    bool processMnps;
    bool processIndels;
    bool processAll;
    bool splitMnps;
    bool removeGenotypes;
    bool storeReducedAlts;
};

} // namespace vcfCTools

#endif // VARIANT_H
