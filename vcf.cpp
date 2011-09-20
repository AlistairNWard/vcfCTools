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
  hasHeader = true;
  hasGenotypes = true;
  processGenotypes = false;
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

// Parse the vcf header.
void vcf::parseHeader(map<string, headerInfoStruct>& oInfo, map<string, headerInfoStruct>& oFormats, vector<string>& oSamples) {
  success = false;
  while(getline(*input, headerLine)) {
    if (headerLine.substr(0, 6) == "##INFO") {success = headerInfo(headerLine, 0);}
    else if (headerLine.substr(0, 8) == "##FORMAT") {success = headerInfo(headerLine, 1);}
    else if (headerLine.substr(0, 2) == "##") {success = headerAdditionalInfo(headerLine);}
    else if (headerLine.substr(0, 1) == "#") {
      success = headerTitles(headerLine);
      getline(*input, headerLine);
    }
    else {success = noHeader();}
    if (!success) {break;}
  }

// headerLine contains the first vcf record.  Process this record
// in preparation for the tools.  Check that this line contains
// data.  Empty files will have a blank line here.
  if (headerLine != "") {
    fromHeader = true;
    record = headerLine;
    string temp = "";
    success = getRecord();
  }

  oInfo    = headerInfoFields;
  oFormats = headerFormatFields;
  oSamples = samples;
}

// Parse information from the info and format descriptors.
bool vcf::headerInfo(string& headerLine, unsigned int headerLineType) {

// Break up the string and fond the info or format tag.
  size_t start = headerLine.find_first_of("<");
  size_t end = headerLine.find_first_of(">");
  string line = headerLine.substr(start, end - start + 1);

  // Find the tag, type and number of expected entries.
  size_t idPosition       = line.find("ID=");
  size_t numberPosition   = line.find("Number=");
  size_t typePosition     = line.find("Type=");
  size_t descPosition     = line.find("Description=");

  headerInfoStruct infoTag;
  string tag = line.substr(idPosition + 3, numberPosition - idPosition - 4);
  if (idPosition == string::npos || numberPosition == string::npos || 
    typePosition == string::npos || descPosition == string::npos) {
    infoTag.number  = ".";
    infoTag.type    = "unknown";
    infoTag.success = false;
  }
  else {
    string number = line.substr(numberPosition + 7, typePosition - numberPosition - 8);
    string type = line.substr(typePosition + 5, descPosition - typePosition - 6);
    string description = line.substr(descPosition + 12, line.size() - descPosition - 13);
    infoTag.success = true;

// The number field can take an integer value if the number of values for this is
// fixed.  If there is an entry per alternate, this should take the value 'A', one
// entry per genotype takes the value 'G' and if the number is unknown or variable
// then this should be '.'.
    if (number == "A") {
      infoTag.number = "A";
    } else if (number == "G") {
      infoTag.number = "G";
    } else if (number == ".") {
      infoTag.number = ".";
    } else if (isdigit(number[0])) {
      infoTag.number = number;
    } else {
      infoTag.success = false;
      infoTag.number  = ".";
      cerr << "WARNING: Malformed info/format in header for tag: ";
      cerr << tag << ".  Unknown value in Number field." << endl;
    }
    infoTag.type        = type;
    infoTag.description = description;
  }

  if (headerLineType == 0) {
    headerInfoFields[tag] = infoTag;
    headerInfoLine[tag] = headerLine;
  }
  else if (headerLineType == 1) {
    headerFormatFields[tag] = infoTag;
    headerFormatLine[tag] = headerLine;
  }

  return true;
}

// Parse additional information from the header.
bool vcf::headerAdditionalInfo(string& headerLine) {
  if (headerText == "") {headerText = headerLine;}
  else {headerText += "\n" + headerLine;}

  return true;
}

// Parse the header titles containing sample names if genotypes are present.
bool vcf::headerTitles(string& headerLine) {
  headerTitlesText = headerLine;

  vector<string> titles = split(headerLine,'\t');
  if (titles.size() > 8) {
    hasGenotypes = true;
    for (int i = 9; i < titles.size(); i++) {samples.push_back(titles[i]);}
    numberSamples = samples.size();
  }
  else if (titles.size() == 8) {
    hasGenotypes = false;
    processGenotypes = false;
    numberSamples = 0;
  }
  else {
    cerr << "Not all vcf standard fields are available." << endl;
    cerr << endl;
    cerr << "Got:     " << headerLine << endl;
    cerr << "Require: #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO" << endl;
    exit(1);
  }

  return false;
}

// No header is present.
bool vcf::noHeader() {

// Check if any of the other header strings have been populated.  It is possible
// that there is header information, but the title line (#CHROM POS ...) is missing
// or malformed.  If this is the case, terminate with a warning.

  if (!headerText.empty() || headerInfoFields.size() != 0 || headerFormatFields.size() != 0) {
    cerr << "No titles line in the vcf file (#CHROM POS etc.)" << endl;
    cerr << "vcfCTools requires a properly constructed header or no header at all to ensure correct operation." << endl;
    cerr << "Please add the required title line to the vcf file:" << endl;
    cerr << "	" << vcfFilename << endl;
    exit(1);
  }

  hasHeader = false;

  return false;
}

// Get the next record from the vcf file.
bool vcf::getRecord() {

// Read in the vcf record.
  success = true;
  if (fromHeader) {fromHeader = false;}
  else {success = getline(*input, record);}

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
    if (variantRecord.quality == 0 && recordFields[5] != "0") {cerr << "ERROR: Variant quality is not an integer or a floating point number." << endl;}
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
