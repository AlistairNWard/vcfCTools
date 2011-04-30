// ******************************************************
// vcfCTools (c) 2011 Alistair Ward
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ------------------------------------------------------
// Last modified: 18 February 2011
// ------------------------------------------------------
// Define the bed class with all associated operations.
// ******************************************************

#ifndef BED_H
#define BED_H

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

struct bedRecord {
  string referenceSequence;
  int start;
  int end;
  string info;
};

class bed {
  public:
    bed(void);
    ~bed(void);
  public:
    bool openBed(string& filename);
    void closeBed();
    void parseHeader();
    bool getRecord();

  public:
    istream* input;
    ifstream file;
    string bedFilename;

// variant information
    bool success;
    string record;
    string referenceSequence;
    vector<string> referenceSequenceVector;
    map<string, bool> referenceSequences;
    unsigned int numberTargets;
    unsigned int targetLength;
    unsigned int targetVariance;
    bedRecord bRecord;
};

} // namespace vcfCTools

#endif // BED_H
