// ******************************************************
// vcfCTools (c) 2011 Alistair Ward
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ------------------------------------------------------
// Last modified: 18 February 2011
// ------------------------------------------------------
// Create the union of two vcf files.
// ******************************************************

#ifndef TOOL_UNION_H
#define TOOL_UNION_H

#include "tools.h"
#include "vcf.h"
#include "vcfCTools_tool.h"

#include <iostream>
#include <fstream>
#include <string>
#include <getopt.h>
#include <vector>

using namespace std;

namespace vcfCTools {

class unionTool : public AbstractTool {
  public:
    unionTool(void);
    ~unionTool(void);
    int Help(void);
    int Run(int argc, char* argv[]);
    int parseCommandLine(int argc, char* argv[]);

  private:
    void unionVcf(vcf&, vcf&, unsigned int, ostream*);

  private:
    string commandLine;
    vector<string> vcfFiles;
    string outputFile;
    ostream* output;
    string priorityFile;
    unsigned int priority;
};

} // namespace vcfCTools

#endif
