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
    void getInfo(string&, int, string&, vector<string>&);
    void processGenotypeFields(string&);
    void getGenotypeInfo(string&, int, string&, vector<string>&);
    bool parseVcf(string&, unsigned int, bool, ostream*);

  public:
    istream* input;
    ifstream file;
    string vcfFilename;

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
    string filters;
    string info;
    map<string, string> infoTags;

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
