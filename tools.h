// ******************************************************
// vcfCTools (c) 2011 Alistair Ward
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ------------------------------------------------------
// Last modified: 18 February 2011
// ------------------------------------------------------
// Additional tools
// ******************************************************

#ifndef TOOLS_H
#define TOOLS_H

#include <cstdio>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <getopt.h>
#include <stdlib.h>
#include <map>
#include <vector>
#include <algorithm>

#include "split.h"
#include "vcf.h"

using namespace std;
using namespace vcfCTools;

struct storedVariants {
  string record;
  string ref;
  string alt;
  float quality;
  bool hasMultipleAlternates;
  bool isSNP;
  bool isMNP;
  bool isDeletion;
  bool isInsertion;
};

ostream* openOutputFile(string&);
unsigned int setVcfPriority(string&, vector<string>&);
void checkReferenceSequences(vector<string>&, vector<string>&);
void mergeHeaders(vcf&, vcf&, vcf&);
void writeHeader(ostream*, vcf&, bool, string&);
void checkDataSets(vcf&, vcf&);
void writeVcfRecord(unsigned int, vcf&, vcf&, ostream*);
storedVariants setStoredVariant(vcf&);
void compareVariants (vector<storedVariants>&, vector<storedVariants>&, bool, bool, string, ostream*);

#endif
