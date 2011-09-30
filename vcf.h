// ******************************************************
// vcfCTools (c) 2011 Alistair Ward
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ------------------------------------------------------
// Last modified: 18 February 2011
// ------------------------------------------------------
// Define the vcf class with all associated operations.
// ******************************************************

#ifndef VCF_H
#define VCF_H

#include "split.h"
#include "vcf_aux.h"

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <map>

using namespace std;

namespace vcfCTools {

// This structure contains a description of a particular
// variant.  A map containing a structure for each variant
// is created at each variant locus.  General information
// about the locus is contained in the following vInfo
// structure.
struct variantDescription {

  // Variant descriptions.
  string record;
  string referenceSequence;
  string rsid;
  string ref;
  string altString;
  double quality;
  string filters;
  string info;
  bool hasGenotypes;
  string genotypeFormatString;
  string genotypeString;
  
  // Boolean flags describing variant class.
  bool isBiallelicSnp;
  bool isTriallelicSnp;
  bool isQuadallelicSnp;
  bool isMnp;
  bool isInsertion;
  bool isDeletion;

  unsigned int variantClass; //DELETE
};

// Define the vcf class.
class vcf {
  public:
    vcf(void);
    ~vcf(void);

  public:

    // File opening and closing.
    bool openVcf(string);
    void closeVcf();

    // Variant reading and structures.
    bool getRecord();

  public:
    istream* input;
    ifstream file;
    string vcfFilename;
    
// Keep track of when a record is read successfully.
    bool success;
    bool update;

// variant information.
    string record;
    variantDescription variantRecord;
    string referenceSequence;
    vector<string> referenceSequenceVector;
    map<string, bool> referenceSequences;
    int position;

// Genotype information.
    bool hasGenotypes;
    bool processGenotypes;
    vector<string> samples;
    unsigned int numberSamples;
};

} // namespace vcfCTools

#endif // VCF_H
