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
  unsigned int number;
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
    bool getRecord();
    void processInfoFields();
    information getInfo(string&);
    void processGenotypeFields(string&);
    information getGenotypeInfo(string&);
    bool parseVcf(string&, unsigned int, bool, ostream*);
    string getDbsnpInfo();
    string buildRecord(bool);

  public:
    istream* input;
    ifstream file;
    string vcfFilename;

// General information.
    bool dbsnpVcf;
    bool hapmapVcf;

// Header information and text.
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
    string record;
    string referenceSequence;
    vector<string> referenceSequenceVector;
    map<string, bool> referenceSequences;
    int position;
    string rsid;
    string ref;
    string alt;
    float quality;
    string sQuality;
    string filters;
    string info;
    map<string, string> infoTags;
    bool hasMultipleAlternates;
    bool isSNP;
    bool isMNP;
    bool isDeletion;
    bool isInsertion;

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
};

} // namespace vcfCTools

#endif // VCF_H
