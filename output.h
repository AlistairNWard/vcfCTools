// ******************************************************
// vcfCTools (c) 2011 Alistair Ward
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ------------------------------------------------------
// Last modified: 18 February 2011
// ------------------------------------------------------
// Define the output class with all associated operations.
// ******************************************************

#ifndef OUTPUT_H
#define OUTPUT_H

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <map>
#include <vector>

using namespace std;

namespace vcfCTools {

class output {
  public:
    output(void);
    ~output(void);
  public:
    ostream* openOutputFile(string&);
    void flushToBuffer(int, string&);
    void flushOutputBuffer();

  public:
    ostream* outputStream;
    string currentReferenceSequence;
    string outputRecord;
    map<int, vector<string> > outputBuffer;
    map<int, vector<string> >::iterator obIter;
    vector<string>::iterator recordIter;
};

} // namespace vcfCTools

#endif // OUTPUT_H
