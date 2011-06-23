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
  int number;
  string type;
  string description;
  bool success;
};

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

class vcf {
  public:
    vcf(void);
    ~vcf(void);
  public:
    bool openVcf(string);
    void closeVcf();
    void parseHeader();
    bool headerInfo(string&, unsigned int);
    bool headerFiles(string&);
    bool headerAdditionalInfo(string&);
    bool headerTitles(string&);
    bool noHeader();
    bool getRecord(string&);
    void addVariantToStructure();
    unsigned int determineVariantClass(string&, string&);
    bool getVariantGroup(variantGroup&, string&);
    void processInfoFields(string&);
    information getInfo(string&);
    void processGenotypeFields(string&);
    information getGenotypeInfo(string&);
    bool parseVcf(string&, unsigned int, bool, ostream*, bool);
    bool parseVcfGroups(variantGroup&, string&, unsigned int, bool, ostream*, string&);
    string getDbsnpInfo();
    void writeRecord(ostream*);

  public:
    istream* input;
    ifstream file;
    string vcfFilename;

// General information.
    bool dbsnpVcf;

// Keep track of when a record is read successfully.
    bool success;
    bool update;
    bool removeGenotypes;

// Header information and text.
    bool fromHeader;
    bool hasHeader;
    bool processInfo;
    string headerLine;
    map<string, string> headerInfoLine;
    map<string, headerInfoStruct> headerInfoFields;
    map<string, string> headerFormatLine;
    map<string, headerInfoStruct> headerFormatFields;
    string headerText;
    string headerTitlesText;
    unsigned int numberDataSets;
    map<unsigned int, string> includedDataSets;

// variant information
    variantDescription variantRecord;
    int position;
    map<unsigned int, vector<variantDescription> > variants;
    map<unsigned int, vInfo> variantsInformation;
    map<unsigned int, vector<variantDescription> >::iterator variantsIter;
    map<string, string> infoTags;

    string record;
    string referenceSequence;
    vector<string> referenceSequenceVector;
    map<string, bool> referenceSequences;
    bool comparedReferenceSequence;
    string rsid;
    string ref;
    string altString;
    vector<string> alt;
    double quality;
    string sQuality;
    string filters;
    string info;
    bool hasMultipleAlternates;
    vector<bool> isSNP;
    vector<bool> isMNP;
    vector<bool> isDeletion;
    vector<bool> isInsertion;

// Genotype information.
    bool hasGenotypes;
    bool processGenotypes;
    string genotypeFormatString;
    vector<string> samples;
    unsigned int numberSamples;
    vector<string> genotypeFormat;
    vector<string> genotypes;
    bool phasedGenotype;
    map<string, string> genotypeTags;

    // Reference sequence information.
    string fasta;
};

} // namespace vcfCTools

#endif // VCF_H
