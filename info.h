// ******************************************************
// vcfCTools (c) 2011 Alistair Ward
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ------------------------------------------------------
// Last modified: 18 February 2011
// ------------------------------------------------------
// Tools for manipulating the information string.
// ******************************************************

#ifndef INFO_H
#define INFO_H

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

struct infoStruct {
  string tag;
  string number;
  string type;
  string value;
  unsigned int numberValues;
};

class variantInfo {
  public:
    variantInfo(string&, map<string, headerInfoStruct>&);
    ~variantInfo(void);
    void modifyInfo(vector<int>&);
    void retrieveFields();

  public:
    string infoString;
    map<string, headerInfoStruct> headerInfo;
    vector<infoStruct> infoFields;
    vector<infoStruct>::iterator infoIter;
};

} // namespace vcfCTools

#endif // INFO_H
