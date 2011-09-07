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
  removeGenotypes = false;
  hasHeader = true;
  processInfo = false;
  hasGenotypes = true;
  processGenotypes = false;
  numberDataSets = 0;
  comparedReferenceSequence = false;
  fasta = "/d2/data/references/build_37/human_reference_v37.fa";
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
void vcf::parseHeader() {
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
    success = getRecord(temp);
  }
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
    headerInfoLine[tag]   = headerLine;
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

  if (!headerText.empty() || headerInfoFields.size() != 0 || headerFormatFields.size() != 0 || includedDataSets.size() != 0) {
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
bool vcf::getRecord(string& currentReferenceSequence) {

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
  if (position == 0) {
    cerr << "Error processing record." << endl;
    if (position == 0) {cerr << "Variant position is not an integer." << endl;}
    if (quality == 0) {cerr << "Variant quality is not an integer or a floating point number." << endl;}
    cerr << endl;
    cerr << "Record:" << endl << endl << record << endl;
    exit(1);
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

//
bool vcf::getVariantGroup(variantGroup& vg, string& refFa) {
  bool success = true, inGroup = true;
  string alRef, alAlt;
  size_t found;
  unsigned int end;

// Define the variant group structure.
  vg.clear();
  vg.referenceSequence = referenceSequence;
  vg.start = -1;

  //cout << "Variant: " << referenceSequence << ":" << position << endl;
  for (vector<string>::iterator iter = alt.begin(); iter != alt.end(); iter ++) {
    unsigned int start = alignAlternate(referenceSequence, position, ref, *iter, alRef, alAlt, refFa); // vcf_aux.cpp
    found = alRef.find("-");
    if (found == string::npos) {end = alRef.length() + start - 1;}
    else {end = found + start - 1;}
  }

  return success;
}

// Process the info entries.
void vcf::processInfoFields(string& infoString) {
  infoTags.clear();
  vector<string> infoEntries = split(infoString, ';');
  for (vector<string>::iterator iter = infoEntries.begin(); iter != infoEntries.end(); iter++) {
    size_t found = (*iter).find_first_of("=");
    string tag = (*iter).substr(0, found);
    if (found == string::npos) {infoTags[tag] = "0";}
    else {infoTags[tag] = (*iter).substr(found + 1, (*iter).size());}
  }
}

// Parse the genotype format string for this record and put the information into
// genotypeFormat.  Create a vector containing all the genotype information.
// There should be as many genotypes as samples.  Actual processing of the
// individual sample genotype strings is performed when required by
// vcf::getGenotypeInfo.
void vcf::processGenotypeFields(string& genotypeString) {
  genotypeTags.clear();
  vector<string> genotypeEntries = split(genotypeString, ':');

// Check that the genotype string has the correct number of entries.
  if (genotypeString == "./." || genotypeString == ".") {}
  else {
    if (genotypeEntries.size() != genotypeFormat.size()) {
      cerr << "Error processing genotypes." << endl;
      cerr << "The genotype string does not contain the correct number of entries:" << endl;
      cerr << "\tFormat  : " << genotypeFormatString << endl;
      cerr << "\tGenotype: " << genotypeString << endl;
      exit(1);
    }
    vector<string>::iterator formatIter = genotypeFormat.begin();
    for (vector<string>::iterator iter = genotypeEntries.begin(); iter != genotypeEntries.end(); iter++) {
      genotypeTags[*formatIter] = *iter;
      formatIter++;
    }
  }
}

// Get the information for a specific info tag.  Also check that it contains
// the correct number and type of entries.
//information vcf::getInfo(string& tag) {
//  information sInfo;
//  sInfo.tag = tag;
//
//// If this routine has been called and processInfoFields has no, or there are
//// no fields in the info string, terminate the program.  Information can only
//// be retrieved if the info fields have been processed and so entering this
//// routine without having procesed the info fields will result in a failue toi
//// extract information.
//  if (infoTags.size() == 0) {
//    cerr << "Routine vcf::getInfo called while there is no information in the" << endl;
//    cerr << "info string, or the routine processInfoFields has not previously been" << endl;
//    cerr << "called.  Please ensure that information exists in the vcf file" << endl;
//    cerr << "and that processInfoFields has been called." << endl;
//    cerr << endl;
//    exit(1);
//  }
//
//// Check if the tag exists in the header information.  If so,
//// determine the number and type of entries asscoiated with this
//// tag.
//  if (headerInfoFields.count(tag) > 0) {
//    sInfo.number = headerInfoFields[tag].number;
//    sInfo.type   = headerInfoFields[tag].type;
//
//// First check that the tag exists in the information string.  Then split
//// the entry on commas.  For flag entries, do not perform the split.
//    if (infoTags.count(tag) > 0) {
//      if (sInfo.number == 0 && sInfo.type == "Flag") {sInfo.values.push_back("true");}
//      else if (sInfo.number != 0 && sInfo.type == "Flag") {
//        cerr << "Error processing info string." << endl;
//        cerr << "Header inforamtion for entry: " << tag << " lists a flag with a non-zero number of entries." << endl;
//        exit(1);
//      } else {
//        sInfo.values = split(infoTags[tag],",");
//        if (sInfo.number != 0 && sInfo.values.size() != sInfo.number) {
//          cerr << "Error processing info string." << endl;
//          cerr << "Unexpected number of entries for info field " << tag << " at " << referenceSequence << ":" << position << endl;
//          exit(1);
//        }
//      }
//
//    // Requested info tag is not in the info string
//    } else {
//      //sInfo.number = 0;
//      cerr << "Tag: " << tag << " is not present in the info string." << endl;
//      cerr << variantsIter->first << endl;
//      cerr << "Terminating program." << endl;
//      exit(1);
//    }
//  }
//  else {
//    cerr << "Error processing info string." << endl;
//    cerr << "No information in the header for info entry: " << tag << endl;
//    exit(1);
//  }
//
//  return sInfo;
//}
//
//// Get the genotype information.
//information vcf::getGenotypeInfo(string& tag) {
//  information gInfo;
//  if (headerFormatFields.count(tag) > 0) {
//    gInfo.number = headerFormatFields[tag].number;
//    gInfo.type   = headerFormatFields[tag].type;
//    gInfo.values = split(genotypeTags[tag], ",");
//    if (gInfo.values.size() != gInfo.number) {
//      cerr << "Error processing info string." << endl;
//      cerr << "Unexpected number of entries for genotype entry " << tag << " at " << referenceSequence << ":" << position << endl;
//      exit(1);
//    }
//  }
//  else {
//    cerr << "Error processing info string." << endl;
//    cerr << "No information in the header for genotype format entry: " << tag << endl;
//    exit(1);
//  }
//
//  return gInfo;
//}

// Parse through the vcf file until the correct reference sequence is
// encountered and the position is greater than or equal to that requested.
bool vcf::parseVcf(string& compReferenceSequence, unsigned int compPosition, bool write, ostream* output, bool passFilters) {
  while (success && variantRecord.referenceSequence == compReferenceSequence) {
    variantsIter = variants.begin();
    if (variantsIter->first < compPosition) {
      if (write) {
        //OUTPUT
      }
      variants.erase(variantsIter);
      success = getRecord(compReferenceSequence);
    } else {
      break;
    }
  }

  return success;
}

// Parse through the vcf file until the correct reference sequence is
// encountered and then construct groups of variants occupying
// overlapping reference sequence and parse through these until the
// start position of the cluster is greater than or equal to the
// requested value.
bool vcf::parseVcfGroups(variantGroup& vc, string& compReferenceSequence, unsigned int compPosition, bool write, ostream* output, string& refFa) {
  bool success = true;
  if (vc.referenceSequence != compReferenceSequence) {
    while (referenceSequence != compReferenceSequence && success) {
      //if (write) {*output << record << endl;}
      success = getRecord(compReferenceSequence);
    }
  }
  if (success) {success = getVariantGroup(vc, refFa);}
  while ( (vc.referenceSequence == compReferenceSequence) && (vc.end < compPosition) && success) {
    //if (write) {*output << record << endl;}
    success = getVariantGroup(vc, refFa);
  }

  return success;
}

// Write a variant to the output stream.  Depending on the number of
// variants at this locus, different options for writing out are
// available.
void vcf::writeRecord(ostream* output) {
  if (variantsIter->second.size() == 1) {
    *output << variantsIter->second[0].record << endl;
  } else {
    for (vector<variantDescription>::iterator iter = variantsIter->second.begin(); iter != variantsIter->second.end(); iter++) {
      *output << iter->record << endl;
    }
  }
}
