// ******************************************************
// vcfCTools (c) 2011 Alistair Ward
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ------------------------------------------------------
// Last modified: 18 February 2011
// ------------------------------------------------------
// Create a class to handle lists of samples.
// ******************************************************

#ifndef SAMPLES_H
#define SAMPLES_H

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

class samples {
  public:
    samples(void);
    ~samples(void);
    bool openSamplesFile(string&);
    void getSamples(vcf&);

  public:
    string samplesFilename;
    istream* input;
    ifstream file;
    unsigned int noSamples;
    map<string, int> samplesMap;
};

} // namespace vcfCTools

#endif // SAMPLES_H
