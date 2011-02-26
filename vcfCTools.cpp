// ******************************************************
// vcfCTools (c) 2011 Alistair Ward
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ------------------------------------------------------
// Last modified: 18 February 2011
// ------------------------------------------------------
// Toolkit for vcf file manipulation
// ******************************************************

// includes
#include "tool_intersect.h"
#include "tool_stats.h"
#include "tool_validate.h"
#include "vcfCTools_version.h"

#include <cstdio>
#include <iostream>
#include <string>

using namespace vcfCTools;
using namespace std;


// vcfCTools tool list
static const string INTERSECT = "intersect";
static const string STATS     = "stats";
static const string VALIDATE  = "validate";

// help and version
static const string HELP          = "help";
static const string LONG_HELP     = "--help";
static const string SHORT_HELP    = "-h";
static const string VERSION       = "version";
static const string LONG_VERSION  = "--version";
static const string SHORT_VERSION = "-v";

// determine if string is a help constant
static bool IsHelp(char* str) {
  return (str == HELP || str == LONG_HELP || str == SHORT_HELP);
}

// determine if string is a version constant
static bool IsVersion(char* str) {
  return (str == VERSION || str == LONG_VERSION || str == SHORT_VERSION);
}

// Determine the tool.
AbstractTool* CreateTool(const string& arg) {
  if (arg == INTERSECT) return new intersectTool;
  if (arg == STATS    ) return new statsTool;
  if (arg == VALIDATE ) return new validateTool;

  return 0;
}

// Print help information.
int Help(int argc, char* argv[]) {
  // Check if the requested help is for a specific tool.
  if (argc > 2) {
    AbstractTool* tool = CreateTool(argv[2]);
    if (tool) return tool->Help();
  }

  // General help information.
  cout << endl;
  cout << "Usage: vcfCTools [tool] [options]" << endl << endl;
  cout << "Available tools:" << endl;
  cout << "  intersect:\n\tCalculate the intersection of two vcf files (or a vcf and a bed file)." << endl;
  cout << "  stats:\n\tGenerate statistics on a vcf file." << endl;
  cout << "  validate:\n\tValidate a vcf file." << endl;
  cout << endl;
  cout << "vcfCTools help tool for help on a specific tool." << endl << endl;
  return 0;
}

// Print out version information.
int Version(void) {
  cout << endl;
  cout << "vcfCTools   - vcf file manipulation tools." << endl;
  cout << "Version " << VCFCTOOLS_VERSION_MAJOR << "." << VCFCTOOLS_VERSION_MINOR << " - " << VCFCTOOLS_DATE << endl;
  cout << "Alistair Ward, Marth Lab, Department of Biology, Boston College." << endl;
  cout << endl;
  return 0;
}

// vcfCTools tool select.
int main(int argc, char* argv[]) {
  // If no tool is selected, show the help.
  if (argc == 1) return Help(argc, argv);

 //vcfCTools help
  if (IsHelp(argv[1])) return Help(argc, argv);

  //vcfCTools version.
  if (IsVersion(argv[1])) return Version();

  // If a tool is specified, determine and run the tool.  If the tool
  // does not exist, show the help.
  AbstractTool* tool = CreateTool(argv[1]);
  if (tool ) return tool->Run(argc, argv);
  else return Help(argc, argv);
}
