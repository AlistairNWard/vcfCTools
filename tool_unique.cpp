// ******************************************************
// vcfCTools (c) 2011 Alistair Ward
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ------------------------------------------------------
// Last modified: 18 February 2011
// ------------------------------------------------------
// Calculate the unique fraction of a vcf file.
// ******************************************************

#include "tool_unique.h"

using namespace std;
using namespace vcfCTools;

// intersectTool imlementation.
uniqueTool::uniqueTool(void)
  : AbstractTool()
{}

// Destructor.
uniqueTool::~uniqueTool(void) {}

// Help
int uniqueTool::Help(void) {
  cout << "Unique  help" << endl;
  cout << "Usage: ./vcfCTools unique [options]." << endl;
  cout << endl;
  cout << "Options:" << endl;
  cout << "  -h, --help" << endl;
  cout << "	display intersect help." << endl;
  cout << "  -i, --in" << endl;
  cout << "	input vcf files.." << endl;
  cout << endl;
  cout << "  -o, --out" << endl;
  cout << "	output vcf file." << endl;
  cout << endl;

  return 0;
}

// Parse the command line and get all required and optional arguments.
int uniqueTool::parseCommandLine(int argc, char* argv[]) {
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

      // Output vcf file.
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
  if (vcfFiles.size() != 2) {
    cerr << "Two vcf files are required for performing the union." << endl;
    exit(1);
  }

  return 0;
}

//Calculate the unique fraction.
void uniqueTool::uniqueVcf(vcf& v1, vcf& v2, ostream* output) {
  bool success1 = v1.getRecord();
  bool success2 = v2.getRecord();
  string currentReferenceSequence = v1.referenceSequence;

// If the end of the first vcf file is reached, it can
// have no more unique records, so terminate.
  while (success1) {

// If the end of the second file is reached, output all
// of the remaining records from the first vcf file as
// they must all be unique.
    if (!success2) {
      *output << v1.record << endl;
      success1 = v1.getRecord();
    }

    if (v1.referenceSequence == v2.referenceSequence && v1.referenceSequence == currentReferenceSequence) {
      if (v1.position == v2.position) {
        success1 = v1.getRecord();
        success2 = v2.getRecord();
      }
      else if (v2.position > v1.position) {success1 = v1.parseVcf(v2.referenceSequence, v2.position, true, output);}
      else if (v1.position > v2.position) {success2 = v2.parseVcf(v1.referenceSequence, v1.position, true, NULL);}
    }
    else {
      if (v1.referenceSequence == currentReferenceSequence) {success1 = v1.parseVcf(v2.referenceSequence, v2.position, true, output);}
      else if (v2.referenceSequence == currentReferenceSequence) {success2 = v2.parseVcf(v1.referenceSequence, v1.position, true, NULL);}

// If the last record for a reference sequence is the same for both vcf
// files, they will both have referenceSequences different from the
// current reference sequence.  Change the reference sequence to reflect
// this and proceed.
      else {
        if (v1.referenceSequence != v2.referenceSequence) {
          cerr << "ERROR: Reference sequences for both files are unexpectedly different." << endl;
          cerr << "Check that both files contain records for the following reference sequences:" << endl;
          cerr << "	" << v1.referenceSequence << " and " << v2.referenceSequence << endl;
          exit(1);
        }
      }
      currentReferenceSequence = v1.referenceSequence;
    }
  }
}

// Run the tool.
int uniqueTool::Run(int argc, char* argv[]) {
  int getOptions = uniqueTool::parseCommandLine(argc, argv);

  output = openOutputFile(outputFile);

  vcf v1; // Define vcf object.
  vcf v2; // Define vcf object.

// Open the vcf files.
  v1.openVcf(vcfFiles[0]);
  v2.openVcf(vcfFiles[1]);

// Read in the header information.
  v1.parseHeader();
  v2.parseHeader();

// Make it clear to the user which unique fraction is being
// calculated.  It is always the first vcf file inputted.
  cout << "Generating records unique to: " << vcfFiles[0] << endl;

// Check that the header for the two files contain the same samples.
  if (v1.samples != v2.samples) {
    cerr << "vcf files contain different samples (or sample order)." << endl;
    exit(1);
  }
  else {
    string taskDescription = "##vcfCTools=generate variants unique to " + vcfFiles[0] + " when compared to " + vcfFiles[1];
    if (v1.hasHeader) {writeHeader(output, v1, true, taskDescription);} // tools.cpp
    else {writeHeader(output, v2, true, taskDescription);} // tools.cpp
  }

// Calculate the unique fraction.
  uniqueTool::uniqueVcf(v1, v2, output);

// Check that the input files had the same list of reference sequences.
// If not, it is possible that there were some problems.
  checkReferenceSequences(v1.referenceSequenceVector, v2.referenceSequenceVector); // tools.py

// Close the vcf objects.
  v1.closeVcf();
  v2.closeVcf();

  return 0;
}
