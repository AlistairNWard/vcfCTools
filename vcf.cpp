// ******************************************************
// vcfCTools (c) 2011 Alistair Ward
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ------------------------------------------------------
// Last modified: 18 February 2011
// ------------------------------------------------------
// vcfClass describes the vcf class and all operations.
// ******************************************************

#include "vcf.h"

#include <iostream>
#include <string>
#include <fstream>
#include <cstdlib>
#include <ctype.h>
#include <vector>

using namespace std;
using namespace vcfCTools;

// Constructor.
vcf::vcf(void) {
  hasGenotypes     = true;
  processGenotypes = false;
  success          = true;
}

// Destructor.
vcf::~vcf(void)
{}

// Open a vcf file.
bool vcf::openVcf(string filename) {
  vcfFilename = filename;
  if (vcfFilename != "-") {
    file.open(vcfFilename.c_str(), ifstream::in);
    input = &file;
    if (!file.is_open()) {
      cerr << "Failed to open file: " << vcfFilename << endl;
      exit(1);
    }
  }
  else {input = &cin;}
}

// Close the vcf file.
void vcf::closeVcf() {
  if (vcfFilename != "-") {
    file.close();
    if (file.is_open()) {
      cerr << "Failed to close file: " << vcfFilename << endl;
    }
  }
}

// Get the next record from the vcf file.
bool vcf::getRecord() {

// Read in the vcf record.
  //success = true;
  //if (fromHeader) {fromHeader = false;}
  //else {success = getline(*input, record);}
  success = getline(*input, record);

// Return false if no more records remain.
  if (!success) {return false;}

// Break the record up into its individual parts.  Leave the genotype fields
// as a string for now.  If the genotypes require parsing, this can be broken
// up when it is needed.
  vector<string> recordFields = split(record, '\t', 10);

// Resolve the information for this variant and add to a temporary structure.
// This will be added to the map of variants when all information has been
// collated.
  variantRecord.referenceSequence = recordFields[0];
  position                        = atoi(recordFields[1].c_str());
  variantRecord.rsid              = recordFields[2];
  variantRecord.ref               = recordFields[3];
  variantRecord.altString         = recordFields[4];
  variantRecord.quality           = atof(recordFields[5].c_str());
  variantRecord.filters           = recordFields[6];
  variantRecord.info              = recordFields[7];

  // Check that genotypes exist.
  if (recordFields.size() < 9) {
    hasGenotypes = false;
    variantRecord.hasGenotypes = false;
  } else {
    hasGenotypes = true;
    variantRecord.hasGenotypes = true;
    variantRecord.genotypeFormatString = recordFields[8];
    variantRecord.genotypeString = recordFields[9];
  }

  // If the position is not an integer, the conversion to an integer will have
  // failed and position = 0.  In this case, terminate with an error.
  if (position == 0 || variantRecord.quality == 0) {
    if (position == 0) {cerr << "ERROR: Unable to process variant position (not an integer)." << endl;}
    if (variantRecord.quality == 0 && recordFields[5] != "0" && recordFields[5] != ".") {
      cerr << "ERROR: Variant quality is not an integer or a floating point number." << endl;
    }
  }

// Add the reference sequence to the map.  If it didn't previously
// exist append the reference sequence to the end of the list as well. 
// This ensures that the order in which the reference sequences appeared
// in the header can be preserved.
  if (referenceSequences.count(variantRecord.referenceSequence) == 0) {
    referenceSequences[variantRecord.referenceSequence] = true;
    referenceSequenceVector.push_back(variantRecord.referenceSequence);
  }

  return success;
}
