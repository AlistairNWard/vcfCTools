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

// Calculate a factorial.
unsigned int fact(unsigned int& x) {
  unsigned int y;
  unsigned int z = 1;

  for (y = 0; y < x; y++) {z = z * (y + 1);}

  return z;
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

  if (removeGenotypes) {
    vector<string> titleFields = split(v.headerTitlesText, '\t');
    string headerTitles = titleFields[0];;
    for (unsigned int i = 1; i < 8; i++) {headerTitles = headerTitles + "\t" + titleFields[i];}
    *output << headerTitles << endl;
  }
  else {*output << v.headerTitlesText << endl;}
}
