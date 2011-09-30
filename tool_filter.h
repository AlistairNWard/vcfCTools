// ******************************************************
// vcfCTools (c) 2011 Alistair Ward
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ------------------------------------------------------
// Last modified: 1 March 2011
// ------------------------------------------------------
// Filter the vcf file on user specific criteria.  vcf
// records that fail the filters have the filter field
// populated with a semi-colon seperated list of failed
// filters.
// ******************************************************

#ifndef TOOL_FILTER_H
#define TOOL_FILTER_H

#include "header.h"
#include "info.h"
#include "output.h"
#include "samples.h"
#include "tools.h"
#include "variant.h"
#include "vcf.h"
#include "vcfCTools_tool.h"

#include <iostream>
#include <fstream>
#include <string>
#include <getopt.h>

using namespace std;

namespace vcfCTools {

class filterTool : public AbstractTool {
  public:
    filterTool(void);
    ~filterTool(void);
    int Help(void);
    int Run(int argc, char* argv[]);
    int parseCommandLine(int argc, char* argv[]);
    vector<string> checkInfoFields(vcfHeader&, vcf&, string&);
    void filter(variant&);
    void performFilter(vcf&, int, variantDescription&);

  private:
    string commandLine;
    string vcfFile;
    string outputFile;
    double filterQualityValue;
    string removeInfoString;
    vector<string> removeInfoList;
    string keepInfoFields;
    vector<string> keepInfoList;
    string stripInfo;
    vector<string> stripInfoList;
    bool writeRecord;
    bool conditionalFilter;
    string currentReferenceSequence;
    string filterString;
    unsigned int genotypePosition;
    string samplesListFile;

    // Boolean flags.
    bool appliedFilters;
    bool cleardbSnp;
    bool filterFail;
    bool filterQuality;
    bool findHets;
    bool keepRecords;
    bool markPass;
    bool processComplex;
    bool processIndels;
    bool processMnps;
    bool processRearrangements;
    bool processSnps;
    bool processSvs;
    bool removeInfo;
    bool removeGenotypes;
    bool stripRecords;
    bool splitMnps;
    bool useSampleList;
};

} // namespace vcfCTools

#endif
