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

#include "tools.h"
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
    void filter(vcf&);

  private:
    string commandLine;
    string vcfFile;
    string outputFile;
    ostream* output;
    bool filterQuality;
    double filterQualityValue;
    bool markPass;
    bool filterFail;
    bool removeGenotypes;
    bool removeInfo;
    string removeInfoString;
    bool stripRecords;
    string stripInfo;
    vector<string> stripInfoList;
    bool writeRecord;
    unsigned int recordsInMemory;
    bool conditionalFilter;
    string filterString;
    bool findHets;
    unsigned int genotypePosition;
};

} // namespace vcfCTools

#endif
