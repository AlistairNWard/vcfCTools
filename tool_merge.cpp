// ******************************************************
// vcfCTools (c) 2011 Alistair Ward
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ------------------------------------------------------
// Last modified: 18 February 2011
// ------------------------------------------------------
// Template for tool creation
// ******************************************************

#include "tool_merge.h"

using namespace std;
using namespace vcfCTools;

// intersectTool imlementation.
mergeTool::mergeTool(void)
  : AbstractTool()
{}

// Destructor.
mergeTool::~mergeTool(void) {}

// Help
int mergeTool::Help(void) {
  cout << "Merge help" << endl;
  cout << "Usage: ./vcfCTools merge [options]." << endl;
  cout << endl;
  cout << "Options:" << endl;
  cout << "  -h, --help" << endl;
  cout << "	display intersect help." << endl;
  cout << "  -i, --in" << endl;
  cout << "	input vcf files to merge (minimum two files)." << endl;
  cout << "  -o, --out" << endl;
  cout << "	output vcf file." << endl;
  cout << endl;

  return 0;
}

// Parse the command line and get all required and optional arguments.
int mergeTool::parseCommandLine(int argc, char* argv[]) {
  commandLine = argv[0];
  for (int i = 2; i < argc; i++) {
    commandLine += " ";
    commandLine += argv[i];
  }

  int argument; // Counter for getopt.

  // Define the long options.
  while (true) {
  static struct option long_options[] = {
    {"help", no_argument, 0, 'h'},
    {"in", required_argument, 0, 'i'},
    {"out", required_argument, 0, 'o'},

    {0, 0, 0, 0}
  };

    int option_index = 0;
    argument = getopt_long(argc, argv, "hi:o:", long_options, &option_index);

    if (argument == -1) {break;}
    switch (argument) {

      // Input vcf file - required input.
      case 'i':
        vcfFiles.push_back(optarg);
        break;

      // Output vcf file - required input.
      case 'o':
        outputFile = optarg;
        break;

      // Help.
      case 'h':
        return Help();

      //
      case '?':
        cout << "Unknown option: " << argv[optind - 1];
        exit(1);
 
      // default
      default:
        abort ();

    }
  }

// Remaining arguments are unknown, so terminate with an error.
  if (optind < argc - 1) {
    cerr << "Unknown options." << endl;
    exit(1);
  }

// Check that at least two  vcf files were specified.
  if (vcfFiles.size() < 2) {
    cerr << "At least two vcf files must be supplied in order to perform a merge." << endl;
    exit(1);
  }

  return 0;
}

// Run the tool.
int mergeTool::Run(int argc, char* argv[]) {
  int getOptions = mergeTool::parseCommandLine(argc, argv);

  output = openOutputFile(outputFile);
  string taskDescription = "##vcfCtools=merge ";
  vector<string> samples;
  for (vector<string>::iterator iter = vcfFiles.begin(); iter != vcfFiles.end(); iter++) {
    taskDescription += (*iter) + ", ";
  }
  taskDescription.erase(taskDescription.end() - 2, taskDescription.end());

  unsigned int index = 0;
  for (vector<string>::iterator iter = vcfFiles.begin(); iter != vcfFiles.end(); iter++) {
    vcf v; // Create a vcf object.
    v.openVcf(vcfFiles[index]);
    v.parseHeader();

// Store the samples list from the first vcf file.  The samplesList from 
// all other vcf files being merged will be checked against this.
// Also, print out the header.
    if (index == 0) {
      samples = v.samples;
      writeHeader(output, v, true, taskDescription); // tools.py
    }
    else {
      if (v.samples != samples) {cerr << "WARNING: Different samples in file: " << v.vcfFilename << endl;}
    }

// Print out the records.
    while (v.getRecord()) {*output << v.record << endl;}

// Close the vcf file.
    v.closeVcf(); // Close the vcf file
  }

  return 0;
}
