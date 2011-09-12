// ******************************************************
// vcfCTools (c) 2011 Alistair Ward
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ------------------------------------------------------
// Last modified: 18 February 2011
// ------------------------------------------------------
// output describes the output class and all operations.
// ******************************************************

#include "output.h"

using namespace std;
using namespace vcfCTools;

// Constructor.
output::output(void) 
{}

// Destructor.
output::~output(void)
{}

// Open the output file.
ostream* output::openOutputFile(string& outputFile) {
  ostream* outputStream;
  if (outputFile == "") {outputStream = &cout;}
  else {outputStream = new ofstream(outputFile.c_str());}

  return outputStream;
}

// Populate the output buffer with a record.
void output::flushToBuffer(int position, string& referenceSequence) {

  // If the reference sequence of the variant to add to the buffer
  // is not the same as the stored value and there are variants in the
  // buffer, flush the buffer to the output file.
  if (outputBuffer.size() != 0 && currentReferenceSequence != referenceSequence) {
    currentReferenceSequence = referenceSequence;
    for (obIter = outputBuffer.begin(); obIter != outputBuffer.end(); obIter++) {
      for (recordIter = obIter->second.begin(); recordIter != obIter->second.end(); recordIter++) {
        *outputStream << *recordIter << endl;
      }
      outputBuffer.erase(obIter);
    }
  }

  // If the output buffer contains more than 1000 entries, flush the
  // first entry to the output and erase it from the buffer.
  if (outputBuffer.size() > 1000) {
    obIter = outputBuffer.begin();
    for (recordIter = obIter->second.begin(); recordIter != obIter->second.end(); recordIter++) {
      *outputStream << *recordIter << endl;
    }
    outputBuffer.erase(obIter);
  }

  // Add the new built record to the buffer.
  outputBuffer[position].push_back(outputRecord);
}

// Clear all entries out of the output buffer.
void output::flushOutputBuffer() {
  for (obIter = outputBuffer.begin(); obIter != outputBuffer.end(); obIter++) {
    for (recordIter = obIter->second.begin(); recordIter != obIter->second.end(); recordIter++) {
      *outputStream << *recordIter << endl;
    }
    outputBuffer.erase(obIter);
  }
}
