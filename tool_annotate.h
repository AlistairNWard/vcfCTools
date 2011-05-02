// ******************************************************
// vcfCTools (c) 2011 Alistair Ward
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ------------------------------------------------------
// Last modified: 18 February 2011
// ------------------------------------------------------
// Annotate a vcf file with dbsnp or hapmap membership.
// The input dbsnp or hapmap files need to be in vcf
// format.
// ******************************************************

#ifndef TOOL_ANNOTATE_H
#define TOOL_ANNOTATE_H

#include "tools.h"
#include "vcf.h"
#include "bed.h"
#include "vcfCTools_tool.h"

#include <iostream>
#include <fstream>
#include <string>
#include <getopt.h>

using namespace std;

namespace vcfCTools {

class annotateTool : public AbstractTool {
  public:
    annotateTool(void);
    ~annotateTool(void);
    int Help(void);
    int Run(int argc, char* argv[]);
    int parseCommandLine(int argc, char* argv[]);

  private:
    string commandLine;
    string vcfFile;
    string outputFile;
    string annVcfFile;
    string bedFile;
    string currentReferenceSequence;
    string newRecord;
    bool annotateDbsnp;
    bool annotateVcf;
    bool annotateBed;
    ostream* output;
    unsigned int recordsInMemory;

    // Boolean flags.
    bool processSnps;
    bool processMnps;
    bool processIndels;
};

} // namespace vcfCTools

#endif
