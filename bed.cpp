// ******************************************************
// vcfCTools (c) 2011 Alistair Ward
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ------------------------------------------------------
// Last modified: 18 February 2011
// ------------------------------------------------------
// bedClass describes the bed class and all operations.
// ******************************************************

#include "bed.h"

using namespace std;
using namespace vcfCTools;

// Constructor.
bed::bed(void) {
  numberTargets = 0;
  targetLength = 0;
  targetVariance = 0;
}

// Destructor.
bed::~bed(void)
{}

// Open a bed file.
bool bed::openBed(string& bedFilename) {
  if (bedFilename != "-") {
    file.open(bedFilename.c_str(), ifstream::in);
    input = &file;
    if (!file.is_open()) {
      cerr << "Failed to open file: " << bedFilename << endl;
      exit(1);
    }
  }
  else {cerr << "bed file must be provided.  Cannot read from stdin." << endl;}
}

// Close the bed file.
void bed::closeBed() {
  if (bedFilename != "-") {
    file.close();
    if (file.is_open()) {
      cerr << "Failed to close file: " << bedFilename << endl;
    }
  }
}

// Parse the bed header.
void bed::parseHeader() {
  success = getRecord();
  while (success && record.substr(0, 1) == "#") {
    success = getline(*input, record);
  }
}

// Get the next record from the vcf file.
bool bed::getRecord() {
  success = getline(*input, record);

// Return false if no more records remain.
  if (!success) {return false;}

  vector<string> recordFields = split(record, '\t');

// Populate the variant values.
  bRecord.referenceSequence = recordFields[0];
  bRecord.start = atoi(recordFields[1].c_str()) + 1;
  bRecord.end   = atoi(recordFields[2].c_str());
  if (recordFields.size() > 3) {bRecord.info = recordFields[3];}

// Check that the start and end coordinates define a valid interval.
  if ( (bRecord.end - bRecord.start) < 0 ) {
    cerr << "Invalid target interval:" << endl;
    cerr << "\t" << record << endl;
    exit(1);
  }

// Update statistics on the targets.
  numberTargets++;
  targetLength += (bRecord.end - bRecord.start);
  targetVariance = 0;

// Add the reference sequence to the map.  If it didn't previously
// exist append the reference sequence to the end of the list as well. 
// This ensures that the order in which the reference sequences appeared
// in the header can be preserved.
  if (referenceSequences.count(bRecord.referenceSequence) == 0) {
    referenceSequences[bRecord.referenceSequence] = true;
    referenceSequenceVector.push_back(bRecord.referenceSequence);
  }

  return true;
}
