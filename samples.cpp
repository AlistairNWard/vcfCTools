// ******************************************************
// vcfCTools (c) 2011 Alistair Ward
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ------------------------------------------------------
// Last modified: 18 February 2011
// ------------------------------------------------------
// Class for handling lists of samples.
// ******************************************************

#include "samples.h"

using namespace std;
using namespace vcfCTools;

// Constructor
samples::samples(void) {
  noSamples = 0;
};

// Desctructo.
samples::~samples(void) {};

// Open the samples list file.
bool samples::openSamplesFile(string& filename) {
  samplesFilename = filename;
  file.open(samplesFilename.c_str(), ifstream::in);
  input = &file;
  if (!file.is_open()) {
    cerr << "Failed to open file: " << samplesFilename << endl;
    exit(1);
  }
}

//
void samples::getSamples(vcf& v) {
  string sampleName;

  // Get the samples from the input file.
  while (getline(*input, sampleName)) {
    samplesMap[sampleName] = -1;
    noSamples++;
  }

  // Loop over all of the samples in the vcf file and if the sample
  // appears in the map, set the value equal to the position of this
  // sample in the list of genotype strings.
  unsigned int count = 0;
  for (vector<string>::iterator iter = v.samples.begin(); iter != v.samples.end(); iter++) {
    if (samplesMap.count(*iter) > 0) {samplesMap[*iter] = count;}
    count++;
  }

  // Now loop over all of the provided samples and check if there are
  // any that are not present in the vcf file.  If the map value is
  // -1 this indicates that it was not seen in the list of samples
  // from the vcf file.
  for (map<string, int>::iterator iter = samplesMap.begin(); iter != samplesMap.end(); iter++) {
    if (iter->second == -1) {
      cerr << "WARNING: Sample " << iter->first << " is not present in the vcf file." << endl;
      samplesMap.erase(iter);
      noSamples--;
    }
  }
}
