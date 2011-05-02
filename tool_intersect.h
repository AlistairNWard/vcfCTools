// ******************************************************
// vcfCTools (c) 2011 Alistair Ward
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ------------------------------------------------------
// Last modified: 18 February 2011
// ------------------------------------------------------
// Calculate the intersection of two vcf files, or a vcf
// file and a bed file.
// ******************************************************

#ifndef TOOL_INTERSECT_H
#define TOOL_INTERSECT_H

#include <cstdio>
#include <iostream>
#include <string>
#include <getopt.h>
#include <stdlib.h>

#include "bed.h"
#include "bedStructure.h"
#include "intersect.h"
#include "tools.h"
#include "variant.h"
#include "vcf.h"
#include "vcfCTools_tool.h"

using namespace std;

namespace vcfCTools {

class intersectTool : public AbstractTool {
  public:
    intersectTool(void);
    ~intersectTool(void);
    int Help(void);
    int Run(int argc, char* argv[]);
    int parseCommandLine(int argc, char* argv[]);

  private:
    string commandLine;
    string bedFile;
    vector<string> vcfFiles;
    string outputFile;
    ostream* output;

    unsigned int recordsInMemory;
    bool passFilters;
    bool findCommon;
    bool findUnion;
    bool findUnique;
    string writeFrom;
    string currentReferenceSequence;

    // Boolean flags.
    bool processSnps;
    bool processMnps;
    bool processIndels;
};

} // namespace vcfCTools

#endif
