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

#include "info.h"
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
    vector<string> checkInfoFields(vcf&, string&);
    void filter(vcf&, variant&);
    void performFilter(vcf&, int, variantDescription&);

  private:
    string commandLine;
    string vcfFile;
    string outputFile;
    ostream* output;
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
    bool splitMnps;
    bool processSnps;
    bool processMnps;
    bool processIndels;
    bool cleardbSnp;
    bool filterQuality;
    bool markPass;
    bool filterFail;
    bool keepRecords;
    bool removeGenotypes;
    bool removeInfo;
    bool stripRecords;
    bool findHets;
    bool useSampleList;
};

} // namespace vcfCTools

#endif
