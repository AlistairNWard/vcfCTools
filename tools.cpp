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
void writeHeader(ostream* output, vcf& v, bool writeGenotypes, string& description) {
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

  if (!writeGenotypes) {
    vector<string> titleFields = split(v.headerTitlesText, '\t');
    string headerTitles = titleFields[0];
    for (unsigned int i = 0; i < 8; i++) {headerTitles = headerTitles + "\t" + titleFields[i];}
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

// Write out a vcf record to file.  The record written depends on the
// value of 'priority' and could therefore be the record from either
// of the vcf files, or a combination of them.
void writeVcfRecord(unsigned int priority, vcf& v1, vcf& v2, ostream* output) {
  if (priority == 0) {
    if (v1.quality >= v2.quality) {*output << v1.record << endl;}
    else {*output << v2.record << endl;}
  }
  else if (priority == 1) {*output << v1.record << endl;}
  else if (priority == 2) {*output << v2.record << endl;}
  else if (priority == 3) {

// Define the missing entry values (depends on the number of data sets
// in the file).
    string info = "";
    string missingEntry1 = ".";
    string missingEntry2 = ".";
    for (unsigned int i = 0; i < v1.numberDataSets; i++) {missingEntry1 += "/.";}
    for (unsigned int i = 0; i < v2.numberDataSets; i++) {missingEntry2 += "/.";}
    map<string, string> secondList = v2.infoTags;

// Build up the info field.
    for (map<string, string>::iterator iter = v1.infoTags.begin(); iter != v1.infoTags.end(); iter++) {
      string tag = (*iter).first;
      if (secondList.count(tag) != 0) {
        if (v1.headerInfoFields[tag].type != "Flag") {info += tag + "=" + v1.infoTags[tag] + "/" + v2.infoTags[tag] + ";";}
        secondList.erase(tag);
      }
      else if (v1.headerInfoFields[tag].type != "Flag") {info += tag + "=" + v1.infoTags[tag] + "/" + missingEntry2 + ";";}
    }

// Now include the info tags that are not populated in the first vcf file.
    for (map<string, string>::iterator iter = secondList.begin(); iter != secondList.end(); iter++) {
      string tag = (*iter).first;
      if (v2.headerInfoFields[tag].type != "Flag") {info += tag + "=" + missingEntry1 + "/" + v2.infoTags[tag] + ";";}
    }

// Build the complete record.
    info.erase(info.end() - 1, info.end());
    *output << v1.referenceSequence << "	" << v1.position << "	" << v1.rsid << "	" << v1.ref;
    *output << "	" << v1.alt << "	" << v2.alt << "	" << v1.quality << "/" << v2.quality;
    *output << "	.	" << info << endl;
  }
  else {
    cerr << "Unknown file priority." << endl;
    exit(1);
  }
}

// Populate the storedVariants structure with variant information.
storedVariants setStoredVariant(vcf& v) {
  storedVariants s;
  s.record = v.record;
  s.quality = v.quality;
  s.ref = v.ref;
  s.alt = v.alt;
  s.hasMultipleAlternates = v.hasMultipleAlternates;
  s.isSNP = v.isSNP;
  s.isMNP = v.isMNP;
  s.isDeletion = v.isDeletion;
  s.isInsertion = v.isInsertion;

  return s;
}

// Compare the variants at a single locus.
void compareVariants (vector<storedVariants>& var1, vector<storedVariants>& var2, bool findUnique, bool findUnion, string writeFrom, ostream* output){
  vector<storedVariants>::iterator iter1;
  vector<storedVariants>::iterator iter2;
  for (iter1 = var1.begin(); iter1 != var1.end(); iter1++) {
    bool unique = true;
    for (iter2 = var2.begin(); iter2 != var2.end(); iter2++) {
      if (iter1->ref == iter2->ref && iter1->alt == iter2->alt) {
        unique = false;
        if (!findUnique) {
          if (writeFrom == "a") {*output << iter1->record << endl;}
          else if (writeFrom == "b") {*output << iter2->record << endl;}
          else if (writeFrom == "q") {
            if (iter1->quality >= iter2->quality) {*output << iter1->record << endl;}
            else {*output << iter2->record << endl;}
          }
        }
        var2.erase(iter2);
        break;
      }
    }
    // If the variant from the first vcf file was not present in the second file, the Boolean
    // unique will be true.  If looking for the union of the files or the variants unique to
    // the first file, this variant should be written out.
    if (unique && ( (findUnique && writeFrom == "a") || findUnion) ) {*output << iter1->record << endl;}
  }
  // All elements of indelsAtLocus2 (the indels present in the second file), were deleted if the
  // same variant was present in the first file.  This means that all variants remaining in this
  // vector are unique to the second vcf file and should be written out if the union or variants
  // unique to the second file were requested.
  if (findUnion || (findUnique && writeFrom == "b")) {
    for (iter2 = var2.begin(); iter2 != var2.end(); iter2++) {*output << iter2->record << endl;}
  }
}
