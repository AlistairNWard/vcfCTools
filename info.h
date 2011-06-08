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

class variantInfo {
  public:
    variantInfo(void);
    ~variantInfo(void);
    void processInfoFields(string&);
    void getInfo(string, string&, int);
    vector<string> buildAltInfo(string&, int, int);

  public:
    string tag;
    unsigned int number;
    string type;
    vector<string> values;
    map<string, string> infoTags;
    map<string, headerInfoStruct> header;
};

} // namespace vcfCTools

#endif // INFO_H
