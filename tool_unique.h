// ******************************************************
// vcfCTools (c) 2011 Alistair Ward
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ------------------------------------------------------
// Last modified: 18 February 2011
// ------------------------------------------------------
// Calculate the unique fraction of a vcf file.
// ******************************************************

#ifndef TOOL_UNIQUE_H
#define TOOL_UNIQUE_H

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

class uniqueTool : public AbstractTool {
  public:
    uniqueTool(void);
    ~uniqueTool(void);
    int Help(void);
    int Run(int argc, char* argv[]);
    int parseCommandLine(int argc, char* argv[]);

  private:
    void uniqueVcf(vcf&, vcf&, ostream*);

  private:
    string commandLine;
    vector<string> vcfFiles;
    string outputFile;
    ostream* output;
};

} // namespace vcfCTools

#endif
