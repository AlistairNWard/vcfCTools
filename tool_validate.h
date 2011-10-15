// ******************************************************
// vcfCTools (c) 2011 Alistair Ward
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ------------------------------------------------------
// Last modified: 18 February 2011
// ------------------------------------------------------
// vcf file validation.
//
// Check for missing data, incomplete data, or data that
// is inconsistent with the information in the header.
// ******************************************************

#ifndef TOOL_VALIDATE_H
#define TOOL_VALIDATE_H

#include <cstdio>
#include <iostream>
#include <string>
#include <getopt.h>
#include <stdlib.h>

#include "genotype_info.h"
#include "header.h"
#include "info.h"
#include "symbolic_alternates.h"
#include "variant.h"
#include "vcf.h"
#include "vcfCTools_tool.h"

using namespace std;

namespace vcfCTools {

class validateTool : public AbstractTool {
  public:
    validateTool(void);
    ~validateTool(void);
    int Help(void);
    int Run(int argc, char* argv[]);
    int parseCommandLine(int argc, char* argv[]);
    //void validateAlternateAlleles(vcfHeader&, variant&);

  private:
    bool error;
    string commandLine;
    string currentReferenceSequence;
    string vcfFile;
};

} // namespace vcfCTools

#endif
