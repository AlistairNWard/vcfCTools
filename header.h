// ******************************************************
// vcfCTools (c) 2011 Alistair Ward
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ------------------------------------------------------
// Last modified: 18 February 2011
// ------------------------------------------------------
// Tools for manipulating the vcf header.
// ******************************************************

#ifndef HEADER_H
#define HEADER_H

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <map>
#include <vector>

#include "split.h"

using namespace std;

namespace vcfCTools {

// Define header structures.
struct headerInfo {
  bool success;
  string description;
  string number;
  string type;
};

// Define the header class.
class vcfHeader {
  public:
    vcfHeader(void);
    ~vcfHeader(void);
    void parseAdditionalInfo();
    void parseHeader(istream*);
    void parseInfo(unsigned int);
    void parseTitles();
    void writeHeader(ostream*, bool, string&);

  public:
    unsigned int numberSamples;
    string line;
    string text;
    map<string, string> altLine;
    map<string, headerInfo> altFields;
    map<string, string> infoLine;
    map<string, headerInfo> infoFields;
    map<string, string> filterLine;
    map<string, headerInfo> filterFields;
    map<string, string> formatLine;
    map<string, headerInfo> formatFields;
    vector<string> samples;
};

} // namespace vcfCTools

#endif // HEADER_H
