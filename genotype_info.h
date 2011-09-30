// ******************************************************
// vcfCTools (c) 2011 Alistair Ward
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ------------------------------------------------------
// Last modified: 18 February 2011
// ------------------------------------------------------
// Tools for manipulating the sample genotypes.
// ******************************************************

#ifndef GENOTYPE_INFO_H
#define GENOTYPE_INFO_H

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <map>
#include <vector>

#include "header.h"
#include "split.h"
#include "tools.h"
#include "vcf.h"

using namespace std;

namespace vcfCTools {

struct genoStruct {
  unsigned int ID;
  string number;
  string type;
  vector<string> values;
};

class genotypeInfo {
  public:
    genotypeInfo(string, string);
    ~genotypeInfo(void);
    void checkTypes(string&, int&, bool&);
    bool getAlleles();
    void modifyGenotypes(vcfHeader&, vector<int>&);
    void processFormats(vcfHeader&);
    void validateGenotypes(vcfHeader&, string&, int&, unsigned int&, vector<string>&, bool&);

  public:
    bool phased;
    bool unphased;
    unsigned int numberInGeno;
    string genotypeFormat;
    string genotypeString;
    vector<string> formats;
    vector<string> genotypeFormats;
    vector<string> genotypeStrings;
    vector<string> genotypes;
    vector<string> originalAlleles;
    vector<string> values;
    map<string, genoStruct> genotypeFields;
    map<string, genoStruct>::iterator genoIter;
};

} // namespace vcfCTools

#endif // GENOTYPE_INFO_H
