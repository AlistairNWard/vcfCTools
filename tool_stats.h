// ******************************************************
// vcfCTools (c) 2011 Alistair Ward
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ------------------------------------------------------
// Last modified: 18 February 2011
// ------------------------------------------------------
// Generate statistics about the input vcf file.
// ******************************************************

#ifndef TOOL_STATS_H
#define TOOL_STATS_H

#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include <getopt.h>
#include <stdlib.h>

#include "header.h"
#include "output.h"
#include "stats.h"
#include "tools.h"
#include "variant.h"
#include "vcf.h"
#include "vcfCTools_tool.h"

using namespace std;

namespace vcfCTools {

class statsTool : public AbstractTool {
  public:
    statsTool( void );
    ~statsTool( void );
    int Help( void );
    int Run( int argc, char* argv[] );
    int parseCommandLine( int argc, char* argv[] );

  private:
    string commandLine;
    string vcfFile;
    string outputFile;
    string currentReferenceSequence;
    string annotationFlagsString;
    vector<string> annotationFlags;
    string detailedGenotypeQualityString;
    string genotypeQualityString;

    // Boolean flags.
    bool generateAfs;
    bool generateDetailed;
    bool processComplex;
    bool processIndels;
    bool processMnps;
    bool processSnps;
    bool generateSampleStats;
    bool splitMnps;
    bool useAnnotations;
};

} // namespace vcfCTools

#endif
