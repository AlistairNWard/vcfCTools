// ******************************************************
// vcfCTools (c) 2011 Alistair Ward
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ------------------------------------------------------
// Last modified: 18 February 2011
// ------------------------------------------------------
// Tools for manipulating the information string.
// ******************************************************

#include "info.h"

using namespace std;
using namespace vcfCTools;

// Constructor
variantInfo::variantInfo(string& info, map<string, headerInfoStruct>& header) {
  infoString = info;
  headerInfo = header;
};

// Desctructor.
variantInfo::~variantInfo(void) {};

// If alleles have been removed, modify the info field so that all
// fields have the correct number of entries.
void variantInfo::modifyInfo(vector<int>& alleleIDs) {
  string modifiedInfo = "";

  retrieveFields();

  // Loop through each of the info fields and check for any fields that should
  // contain values corresponding to the different alternate alleles.
  infoIter = infoFields.begin();
  for (; infoIter != infoFields.end(); infoIter++) {
    if (infoIter->second.number == "A") {
      vector<string>::iterator vIter   = infoIter->second.values.begin();
      vector<int>::iterator alleleIter = alleleIDs.begin();

      // Skip the reference allele in the alleleIDs.
      alleleIter++;

      modifiedInfo += infoIter->first + "=";
      for (; alleleIter != alleleIDs.end(); alleleIter++) {
        if (*alleleIter != -1) {modifiedInfo += *vIter + ",";}
        vIter++;
      }

      // Strip the trailing ",".
      modifiedInfo = modifiedInfo.substr(0, modifiedInfo.size() - 1);
      modifiedInfo += ";";
    } else if (infoIter->second.number == "0") {
      modifiedInfo += infoIter->first + ";";
    } else {
      modifiedInfo += infoIter->first + "=" + infoIter->second.values[0] + ";";
    }
  }

  // Strip the trailing ";" and then replace the existing info string
  // with the modified one.
  modifiedInfo = modifiedInfo.substr(0, modifiedInfo.size() - 1);
  infoString   = modifiedInfo;
}

// Split the info string into its components and populate the arrays
// to store the information.
void variantInfo::retrieveFields() {
  string tag;
  vector<string> infoArray = split(infoString, ";");
  vector<string>::iterator infoIter = infoArray.begin();
  for (; infoIter != infoArray.end(); infoIter++) {

    // Split the entry on "=" to get the tag.
    vector<string> fields = split(*infoIter, "=");
    infoStruct info;
    tag = fields[0];

    // In order to process the info field, the header information is
    // required.  If this does not exist, terminate the program.
    if (headerInfo.count(tag) == 0) {
      cerr << "ERROR: No header description for info tag: " << tag << endl;
      exit(1);
    }

    info.number = headerInfo[tag].number;
    info.type   = headerInfo[tag].type;
    if (fields.size() > 1) {
      info.values        = split(fields[1], ",");
      info.numberValues  = info.values.size();
    } else {
      info.values.push_back("");
      info.numberValues = 0;
    }
    infoFields[tag] = info;
  }
}

// Parse through the entire info string and ensure that all entries have an
// explanation in the header and that the number and entry types are
// consistent with this header information.
void variantInfo::validateInfo(string& referenceSequence, int& position, unsigned int& noAlts, bool& error) {
  int integerValue;

  // Split the info field into its constituent parts.
  retrieveFields();
  for (infoIter = infoFields.begin(); infoIter != infoFields.end(); infoIter++) {

    // If the info number = A, there should be as many values as there are
    // alternate alleles.  Check that this is the case.
    if (infoIter->second.number == "A") {
      if (infoIter->second.numberValues != noAlts) {
        cerr << "ERROR: Incorrect number of entries at " << referenceSequence << ":";
        cerr << position << ", in info field for: " << infoIter->first << endl;
        error = true;
      }

      // Check the values conform to the expected types.
      checkTypes(referenceSequence, position, error);

    } else if (infoIter->second.number == "G") {

    // If the info number = '.', any number of entries is acceptable, so do not
    // check the number of values.
    } else if (infoIter->second.number == ".") {

      // Check the values conform to the expected types.
      checkTypes(referenceSequence, position, error);

    // If the info number = 0, then this info field is a flag.  Check that there is
    // only a single entry.
    } else if (infoIter->second.number == "0") {
      if (infoIter->second.numberValues != 0 || infoIter->second.type != "Flag") {
        cerr << "ERROR: Incorrect number of entries or not marked as a flag at " << referenceSequence << ":";
        cerr << position << ", in info field for: " << infoIter->first << "." << endl;
        error = true;
      }

      // Check the values conform to the expected types.
      checkTypes(referenceSequence, position, error);

    // Finally, if a specified number of entries is given, check that this number
    // is observed.  Also check that the entry is itself an integer as there are no
    // other allowed types.
    } else {

      // First check that the number is an integer.
      integerValue = atoi( infoIter->second.number.c_str() );
      if (integerValue == 0) {
        cerr << "ERROR: Invalid type in info field: " << infoIter->first << " at " << referenceSequence;
        cerr << ":" << position << "." << endl;
        error = true;
      }

      // Now check that the correct number is observed.
      if (infoIter->second.numberValues != integerValue) {
        cerr << "ERROR: Incorrect number of entries at " << referenceSequence << ":";
        cerr << position << ", in info field for: " << infoIter->first << endl;
        error = true;
      }

      // Check the values conform to the expected types.
      checkTypes(referenceSequence, position, error);
    }
  }
}

// Check that the observed values conform to the expected type.
void variantInfo::checkTypes(string& referenceSequence, int& position, bool& error) {
  int integerValue;
  float floatValue;
  vector<string>::iterator sIter;

  // Now check that the values are of the correct type.
  sIter = infoIter->second.values.begin();
  for (; sIter != infoIter->second.values.end(); sIter++) {

    // Check that integers appear if integers are expected.
    if (infoIter->second.type == "Integer") {
      integerValue = atoi( sIter->c_str() );
      if (integerValue == 0 && *sIter != "0") {
        cerr << "ERROR: Expected integer in info field: " << infoIter->first << " at " << referenceSequence;
        cerr << ":" << position << "." << endl;
        error = true;
      }

    // If floating point values are expected, check that this is seen.
    } else if (infoIter->second.type == "Float") {
      floatValue = atof( sIter->c_str() );
      if (floatValue == 0. && (*sIter != "0" && *sIter != "0." && *sIter != "0.0")) {
        cerr << "ERROR: Expected floats in info field: " << infoIter->first << " at " << referenceSequence;
        cerr << ":" << position << "." << endl;
        error = true;
      }
    } else if (infoIter->second.type == "Flag") {
    } else if (infoIter->second.type == "Character") {

      // Check that the values are of length 1.
      if (sIter->length() != 1) {
        cerr << "ERROR: Expected a single character for info field: " << infoIter->first << "at";
        cerr << referenceSequence << ":" << position << "." << endl;
        error = true;
      }
    } else if (infoIter->second.type == "String") {
    }
  }
}
