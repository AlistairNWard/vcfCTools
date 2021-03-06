// ******************************************************
// vcfCTools (c) 2011 Alistair Ward
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ------------------------------------------------------
// Last modified: 18 February 2011
// ------------------------------------------------------
// Template for tool creation
// ******************************************************

#ifndef TOOL_MERGE_H
#define TOOL_MERGE_H

#include "header.h"
#include "output.h"
#include "tools.h"
#include "variant.h"
#include "vcfCTools_tool.h"
#include "vcf.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <getopt.h>

using namespace std;

namespace vcfCTools {

class mergeTool : public AbstractTool {
  public:
    mergeTool(void);
    ~mergeTool(void);
    int Help(void);
    int Run(int argc, char* argv[]);
    int parseCommandLine(int argc, char* argv[]);

  private:
    string commandLine;
    vector<string> vcfFiles;
    string outputFile;
    string currentReferenceSequence;

    // Boolean flags.
    bool processComplex;
    bool processIndels;
    bool processMnps;
    bool processRearrangements;
    bool processSnps;
    bool processSvs;
};

} // namespace vcfCTools

#endif
