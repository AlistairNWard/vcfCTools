// ******************************************************
// vcfCTools (c) 2011 Alistair Ward
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ------------------------------------------------------
// Last modified: 18 February 2011
// ------------------------------------------------------
// Define structures to be used by all routines.
// ******************************************************

#ifndef STRUCTURES_H
#define STRUCTURES_H

using namespace std;

namespace vcfCTools {

// Create astructure to hold all of the flags required to determine the intersection
// operations to be performed.
//struct intFlags {
//  bool annotate;
//  bool findCommon;
//  bool findUnion;
//  bool findUnique;
//  bool writeFromFirst;
//  bool sitesOnly;
//  bool whollyWithin;
//};

// Define a structure containing the different variant
// types.
struct variantType {
  bool isBiallelicSnp;
  bool isTriallelicSnp;
  bool isQuadallelicSnp;
  bool isMnp;
  bool isInsertion;
  bool isDeletion;
  bool isComplex;
};

// Define a structure that contains information about a
// particular locus.  This structure is used for variants
// still in their original form and at the position they
// appear in the vcf file.  It is possible that these
// variants can be reduced (e.g. CC -> CG could be rewritten
// as a C -> G SNP at the next position), but these are kept
// in a separate structure for use in intersections etc.
//struct originalVariants {

  // General information.
//  bool hasMultipleRecords;
//  unsigned int numberOfRecordsAtLocus;
//  int position;
//  int maxPosition;
//  double quality;
//  string info;
//  string filters;
//  string referenceSequence;
//  vector<int> reducedPosition;
 
  // Ref and alt allele information.
//  unsigned int numberAlts;
//  string ref;
//  string rsid;
//  string altString;
//  vector<bool> filtered;
//  vector<string> reducedRef;
//  vector<string> alts;
//  vector<string> reducedAlts;
//  vector<variantType> type;

  // Genotype information.
//  bool hasGenotypes;
//  string genotypeFormat;
//  string genotypes;
//};

// Define a structure that stores the information about the
// modified variants.  For example CG -> CA at position 10 will
// be stored as this in the originalVariants structure, but will
// be stored as a G -> A SNP at position 11 in this structure.
// This structure will point back to the originalVariants
// structure for genotypes and info etc.
//struct reducedVariants {

  // If there are multiple records at this locus, the recordNumber
  // indicates which record this variant comes from.
//  unsigned int recordNumber;
//  int originalPosition;
//  int altID; // If alts were A,G, this will indicate 1 for A and 2 for G etc.
//  string ref;
//  string alt;
//};

// At each locus, many different variants can exist.  The
// variantsAtLocus structure holds all of the variants
// present at this locus.
//struct variantsAtLocus {
//  string referenceSequence;
//  vector<reducedVariants> snps;
//  vector<reducedVariants> mnps;
//  vector<reducedVariants> insertions;
//  vector<reducedVariants> deletions;
//  vector<reducedVariants> complexVariants;
//};

// Define a structure that contains information about the different reference
// sequences encountered.  This includes information such as the number of
// records parsed for each reference sequence and if an intersection is
// performed, information about whether the reference sequences have been
// used in comparisons.
//struct refSeqInfo {
//  unsigned int numberRecords;
//  bool contiguous;
//  bool usedInComparison;
//};

} // namespace vcfCTools

#endif // STRUCTURES_H
