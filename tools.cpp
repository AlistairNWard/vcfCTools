// ******************************************************
// vcfCTools (c) 2011 Alistair Ward
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ------------------------------------------------------
// Last modified: 25 February 2011
// ------------------------------------------------------
// Additional tools.
// ******************************************************

#include "tools.h"

using namespace std;

ostream* openOutputFile(string& outputFile) {
  ostream* output;
  if (outputFile == "") {output = &cout;}
  else {output = new ofstream(outputFile.c_str());}

  return output;
}

// Determine which file has priority for writing out records.
unsigned int setVcfPriority(string& priorityFile, vector<string>& vcfFiles) {
  unsigned int priority;
  if (priorityFile == "") {priority = 0;}
  else if (priorityFile == vcfFiles[0]) {priority = 1;}
  else if (priorityFile == vcfFiles[1]) {priority = 2;}
  else if (priorityFile == "merge") {priority = 3;}
  else {
    cerr << "vcf file give priority must be one of the two input vcf files or merge." << endl;
    exit(1);
  }

  return priority;
}

// Check that the two reference sequence lists are identical.
// If there are a different number or order, the results may
// not be as expected.
void checkReferenceSequences(vector<string>& vector1, vector<string>& vector2) {
  bool errorMessage = false;
  if ( vector1.size() != vector2.size() ) {
    cerr << "WARNING: Input files contain a different number of reference sequences." << endl;
    errorMessage = true;
  }
  else if (vector1 != vector2) {
    cerr <<  "WARNING: Input files contain different or differently ordered reference sequences." << endl;
    errorMessage = true;
  }
  if (errorMessage) {
    cerr << "Results may not be as expected." << endl;
    cerr << "Ensure that input files have the same reference sequences in the same order." << endl;
    cerr << "Reference sequence lists observed were:" << endl;
    cerr << "	vcf 1: ";
    for (vector<string>::iterator iter = vector1.begin(); iter != vector1.end(); iter++) {cerr << (*iter) << ", ";}
    cerr << endl;
    cerr << "	vcf 2: ";
    for (vector<string>::iterator iter = vector2.begin(); iter != vector2.end(); iter++) {cerr << (*iter) << ", ";}
    cerr << endl;
  }
}

// If the union or intersection of two vcf files is being performed
// and the output vcf file is to contain the information from both
// files, the headers need to be merged to ensure that all info and
// format entries have an explanation.
void mergeHeaders(vcf& v1, vcf& v2, vcf& v3) {

// If either file does not have a header, terminate the program.
// In order to merge the headers, the different fields must be
// checked to ensure the files are compatible.
  if (!v1.hasHeader || !v2.hasHeader) {
    cerr << "Both vcf files must have a header in order to merge data sets." << endl;
    exit(1);
  }

  v3.headerInfoFields   = v1.headerInfoFields;
  v3.headerFormatFields = v1.headerFormatFields;
  v3.numberDataSets     = v1.numberDataSets;
  v3.includedDataSets   = v1.includedDataSets;
  v3.headerText         = v1.headerText;
  v3.headerTitlesText   = v1.headerTitlesText;
  v3.headerInfoLine     = v1.headerInfoLine;
  v3.headerFormatLine   = v1.headerFormatLine;

// Merge the info field descriptions.
  for (map<string, headerInfoStruct>::iterator iter = v2.headerInfoFields.begin(); iter != v2.headerInfoFields.end(); iter++) {
    string tag = (*iter).first;
    if (v1.headerInfoFields.count(tag) != 0) {
      if (v1.headerInfoFields[tag].number != v2.headerInfoFields[tag].number || v1.headerInfoFields[tag].type != v2.headerInfoFields[tag].type) {
        cerr << v1.headerInfoFields[tag].number << endl;
        cerr << v1.headerInfoFields[tag].type << endl;
        cerr << v1.headerInfoFields[tag].description << endl;
        cerr << "Input vcf files have different definitions for " << tag << " field." << endl;
        exit(1);
      }
    }
    else {v3.headerInfoFields[tag] = v2.headerInfoFields[tag];}
  }

// Merge the format field descriptions.
  for (map<string, headerInfoStruct>::iterator iter = v2.headerFormatFields.begin(); iter != v2.headerFormatFields.end(); iter++) {
    string tag = (*iter).first;
    if (v1.headerFormatFields.count(tag) != 0) {
      if (v1.headerFormatFields[tag].number != v2.headerFormatFields[tag].number || \
          v1.headerFormatFields[tag].type != v2.headerFormatFields[tag].type) {
        cerr << "Input vcf files have different definitions for " << tag << " field." << endl;
        exit(1);
      }
    }
    else {v3.headerFormatFields[tag] = v2.headerFormatFields[tag];}
  }

// Now check to see if the vcf files contain information from multiple
// records themselves and create an ordered list in which the data
// will appear in the file.  For instance, of the first file has
// already got two sets of data and is being intersected with a file
// with one set of data, the order of data in the new vcf file will be
// the two sets from the first file followed by the second, e.g.
// AB=3/2/4, where the 3 and 2 are from the first file and the 4 is the
// value of AC from the second vcf.  The header will have a ##FILE for
// each of the three files, so the origin if the data can be recovered.
  if (v1.numberDataSets == 0) {
    v3.includedDataSets[v3.numberDataSets + 1] = v1.vcfFilename;
    v3.numberDataSets++;
  }
  if (v2.numberDataSets == 0) {
    v3.includedDataSets[v3.numberDataSets + 1] = v2.vcfFilename;
    v3.numberDataSets++;
  }
  else {
    for (unsigned int i = 0; i < v2.numberDataSets; i++) {
      v3.includedDataSets[v3.numberDataSets + 1] = v2.includedDataSets[i];
      v3.numberDataSets++;
    }
  }
}

// Write header to file.
void writeHeader(ostream* output, vcf& v, bool removeGenotypes, string& description) {
  if (!v.hasHeader) {
    v.headerText = "##fileformat=VCFv4.0\n##source=vcfCtools "; // + version
    v.headerTitlesText = "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO";
  }
  if (v.headerText != "") {*output << v.headerText << endl;}
  *output << description << endl;
  for (map<string, string>::iterator iter = v.headerInfoLine.begin(); iter != v.headerInfoLine.end(); iter++) {
    *output << (*iter).second << endl;
  }
  for (map<string, string>::iterator iter = v.headerFormatLine.begin(); iter != v.headerFormatLine.end(); iter++) {
    *output << (*iter).second << endl;
  }

// Write out a list of files indicating which data set belongs to which file.
  if (v.numberDataSets != 0) {
    for (unsigned int i = 1; i < v.numberDataSets + 1; i++) {
      *output << "##FILE=<ID=" << i << ",\"" << v.includedDataSets[i] << "\">" << endl;
    }
  }

  if (removeGenotypes) {
    vector<string> titleFields = split(v.headerTitlesText, '\t');
    string headerTitles = titleFields[0];;
    for (unsigned int i = 1; i < 8; i++) {headerTitles = headerTitles + "\t" + titleFields[i];}
    *output << headerTitles << endl;
  }
  else {*output << v.headerTitlesText << endl;}
}

// If either of the input files contain multiple data sets (e.g. multiple
// vcf files have undergone intersection or union calculations and all
// information has been retained) and the priority isn't set to 'merge',
// terminate the program.  This is to ensure that the origin of the data
// doesn't get confused.
void checkDataSets(vcf& v1, vcf& v2) {
  if (v1.numberDataSets + v2.numberDataSets != 0) {
    cerr << "\nERROR:" << endl;
    cerr << "input vcf file(s) contain data sets from multiple vcf files." << endl;
    cerr << "Further intersection or union operations must include --priority-file merge" << endl;
    cerr << "Other tools may be incompatible with this format." << endl;
    exit(1);
  }
}

// Build a variant record from its constituent parts.
void buildRecord(int position, variantDescription& v) {
  ostringstream sPosition, sQuality;
  sPosition << position;
  sQuality << v.quality;

  v.record = v.referenceSequence + "	" +
             sPosition.str() + "	" +
             v.rsid + "	" +
             v.ref + "	" +
             v.altString + "	" +
             sQuality.str() + "	" +
             v.filters + "	" + 
             v.info;

  // If genotypes exist, add them to the record.
  if (v.genotypeString != "") {
    v.record += "	" + 
                v.genotypeFormatString + "	" +
                v.genotypeString;
  }
}
