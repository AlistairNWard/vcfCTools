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

#include "vcfCTools_tool.h"

using namespace std;

namespace vcfCTools {

class validateTool : public AbstractTool {
  public:
    validateTool( void );
    ~validateTool( void );
  public:
    int Help( void );
    int Run( int argc, char* argv[] );
    int parseCommandLine( int argc, char* argv[] );

  private:
    string commandLine;
    string vcfFile;
    string currentReferenceSequence;

    // Boolean flags.
    bool processSnps;
    bool processMnps;
    bool processIndels;
};

} // namespace vcfCTools

#endif
