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

#include "vcf.h"

using namespace std;

namespace vcfCTools {

class genotypeInfo {
  public:
    genotypeInfo(string, string, map<string, headerInfoStruct>&);
    ~genotypeInfo(void);
    void modifyGenotypes(vector<int>&);

  public:
    string genotypeFormat;
    string genotypeString;
    vector<string> genotypeStrings;
    vector<string> formats;
    map<string, headerInfoStruct> headerInfo;
};

} // namespace vcfCTools

#endif // GENOTYPE_INFO_H
