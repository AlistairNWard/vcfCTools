// ******************************************************
// vcfCTools (c) 2011 Alistair Ward
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ------------------------------------------------------
// Last modified: 18 February 2011
// ------------------------------------------------------
// Template for tool creation
// ******************************************************

#ifndef TOOL_TEMPLATE_H
#define TOOL_TEMPLATE_H

#include "tools.h"
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
    ostream* output;
};

} // namespace vcfCTools

#endif
