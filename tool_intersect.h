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
#include "tools.h"
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
    //void intersectVcf(vcf&, vcf&, ostream*, bool, bool, bool, string);
    void intersectVcf(vcf&, vcf&);
    //void intersectVariantGroups(vcf&, vcf&, ostream*, bool, bool, bool, string, string&);
    void intersectVariantGroups(vcf&, vcf&, string&);
    //void intersectVcfBed(vcf&, bed&, ostream*);
    void intersectVcfBed(vcf&, bed&);
    //void intersectVariantGroupsBed(vcf&, bed&, ostream*);
    void intersectVariantGroupsBed(vcf&, bed&);

  private:
    string commandLine;
    string bedFile;
    vector<string> vcfFiles;
    string outputFile;
    ostream* output;
    storedVariants s;
    vector<struct storedVariants> snpsAtLocus1;
    vector<struct storedVariants> mnpsAtLocus1;
    vector<struct storedVariants> indelsAtLocus1;
    vector<struct storedVariants> snpsAtLocus2;
    vector<struct storedVariants> mnpsAtLocus2;
    vector<struct storedVariants> indelsAtLocus2;

    bool passFilters;
    bool findCommon;
    bool findUnion;
    bool findUnique;
    bool groupVariants;
    string referenceFasta;
    string writeFrom;
    string currentReferenceSequence;
    unsigned int currentPosition;
};

} // namespace vcfCTools

#endif
