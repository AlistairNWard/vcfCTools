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

struct headerInfoStruct {
  string number;
  string type;
  string description;
  bool success;
};

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




// Original structures and variables that may be disposable.
// When complete, this should all be removed.
struct information {
  string tag;
  unsigned int number;
  string type;
  vector<string> values;
  bool flag;

  information() {flag = false;}
};

struct variantGroup {
  string referenceSequence;
  unsigned int start;
  unsigned int end;
  string ref;
  unsigned int refLength;
  unsigned int noAlts;
  unsigned int noRecords;
  unsigned int noGroups;

  variantGroup() {
    noAlts = 0;
    noRecords = 0;
    noGroups = 0;
  }

  ~variantGroup()
  {}

  void clear() {
    noAlts = 0;
    noRecords = 0;
  }
};

// This structure contains general information about the
// variants described in the variantsDescription structure.
struct vInfo {
  string referenceSequence;

  // Logical flags.
  bool containsBiallelicSnp;
  bool containsMultipleSnp;
  bool containsTriallelicSnp;
  bool containsQuadallelicSnp;
  bool containsMnp;
  bool containsInsertion;
  bool containsDeletion;
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

    // Header parsing.
    void parseHeader(map<string, headerInfoStruct>&, map<string, headerInfoStruct>&, vector<string>&);
    bool headerInfo(string&, unsigned int);
    bool headerAdditionalInfo(string&);
    bool headerTitles(string&);
    bool noHeader();

    // Variant reading and structures.
    bool getRecord();

  public:
    istream* input;
    ifstream file;
    string vcfFilename;
    
// Keep track of when a record is read successfully.
    bool success;
    bool update;

// Header information and text.
    bool fromHeader;
    bool hasHeader;
    string headerLine;
    map<string, string> headerInfoLine;
    map<string, headerInfoStruct> headerInfoFields;
    map<string, string> headerFormatLine;
    map<string, headerInfoStruct> headerFormatFields;
    string headerText;
    string headerTitlesText;

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
