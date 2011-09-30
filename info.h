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

#include "header.h"
#include "split.h"
#include "vcf.h"

using namespace std;

namespace vcfCTools {

struct infoStruct {
  string number;
  string type;
  vector<string> values;
  unsigned int numberValues;
};

class variantInfo {
  public:
    variantInfo(string&);
    ~variantInfo(void);
    void checkTypes(string&, int&, bool&);
    void modifyInfo(vector<int>&, vcfHeader&);
    void retrieveFields(vcfHeader&);
    void validateInfo(vcfHeader&, string&, int&, unsigned int&, bool&);

  public:
    string infoString;
    map<string, infoStruct> infoFields;
    map<string, infoStruct>::iterator infoIter;
};

} // namespace vcfCTools

#endif // INFO_H
