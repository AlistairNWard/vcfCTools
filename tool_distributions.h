// ******************************************************
// vcfCTools (c) 2011 Alistair Ward
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ------------------------------------------------------
// Last modified: 18 February 2011
// ------------------------------------------------------
// Collate information.
// ******************************************************

#ifndef TOOL_DISTRIBUTIONS_H
#define TOOL_DISTRIBUTIONS_H

#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include <getopt.h>
#include <stdlib.h>

#include "header.h"
#include "info.h"
#include "tools.h"
#include "variant.h"
#include "vcf.h"
#include "vcfCTools_tool.h"

using namespace std;

namespace vcfCTools {

struct distributionsStruct {
  unsigned int number;
  map<string, double> secondary;
};

class distributionsTool : public AbstractTool {
  public:
    distributionsTool( void );
    ~distributionsTool( void );
    int Help( void );
    int Run( int argc, char* argv[] );
    int parseCommandLine( int argc, char* argv[] );
    void distributions(vcf&, variant&);
    void performCollate(vcf&, int, variantDescription&);
    void writePrimaryInfo();
    void writeDistributions();

  private:
    string commandLine;
    string vcfFile;
    string outputFile;
    ostream* output;
    string currentReferenceSequence;
    string primaryInfo;
    string secondaryInfoString;
    vector<string> secondaryInfo;
    vector<string> distFields;
    bool secondaryQuality;
    bool usePrimary;
    bool useDistributions;
    bool useDistQ;
    string distString;
    map<string, map<string, unsigned int> > dist;
    map<double, unsigned int> distQ;
    map<string, distributionsStruct> distributionsdInfo;

    // Boolean flags.
    bool processComplex;
    bool processSnps;
    bool processMnps;
    bool processIndels;
};

} // namespace vcfCTools

#endif
