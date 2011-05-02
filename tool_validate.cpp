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

#include "tool_validate.h"
#include "vcf.h"

#include <iostream>
#include <string>

using namespace std;
using namespace vcfCTools;

// validateTool imlementation.
validateTool::validateTool(void)
  : AbstractTool()
{
  currentReferenceSequence = true;
}

// Destructor.
validateTool::~validateTool(void) {}

// Help
int validateTool::Help(void) {
  cout << "Validation help" << endl;
  cout << "Usage: ./vcfCTools validate [options]." << endl;
  cout << endl;
  cout << "Options:" << endl;
  cout << "  -h, --help" << endl;
  cout << "     display intersect help." << endl;
  cout << "  -i, --in" << endl;
  cout << "     input vcf file." << endl;
  cout << "  -o, --output" << endl;
  cout << "     output file." << endl;
  cout << "  -1, --snps" << endl;
  cout << "	analyse SNPs." << endl;
  cout << "  -2, --mnps" << endl;
  cout << "	analyse MNPs." << endl;
  cout << "  -3, --indels" << endl;
  cout << "	analyse indels." << endl;
  return 0;
}

// Parse the command line and get all required and optional arguments.
int validateTool::parseCommandLine(int argc, char* argv[]) {
  commandLine = argv[0];
  for (int i = 2; i < argc; i++) {
    commandLine += " ";
    commandLine += argv[i];
  }

  int argument; // Counter for getopt.
  // Define the long options.
  static struct option long_options[] = {
    {"help", no_argument, 0, 'h'},
    {"in", required_argument, 0, 'i'},
    {"snps", no_argument, 0, '1'},
    {"mnps", no_argument, 0, '2'},
    {"indels", no_argument, 0, '3'},

    {0, 0, 0, 0}
  };

  while (true) {
    int option_index = 0;
    argument = getopt_long(argc, argv, "hi:", long_options, &option_index);

    if (argument == -1)
      break;

    switch (argument) {
      // Input vcf file - required input.
      case 'i':
        vcfFile = optarg;
        break;

      // Analyse SNPs.
      case '1':
        processSnps = true;
        break;

      // Analyse MNPs.
      case '2':
        processMnps = true;
        break;

      // Analyse indels.
      case '3':
        processIndels = true;
        break;
      
      // Help.
      case 'h':
        return Help();

      //
      case '?':
        cout << "Unknown option: " << argv[optind - 1] << endl;
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

// Check that a vcf file was specified.
  if (vcfFile == "") {
    cerr << "A vcf file must be specified (--in, -i)." << endl;
    exit(1);
  }

  return 0;
}

// Run the tool.
int validateTool::Run(int argc, char* argv[]) {
  int getOptions = validateTool::parseCommandLine(argc, argv);

// The info/format fields and the genotypes need to be processed to
// ensure that there are no errors.
  vcf v; // Create a vcf object.
  v.processInfo = true;
  v.processGenotypes = true;

  v.openVcf(vcfFile);
  v.parseHeader();

// Check that all of the info and format fields in the header were
// succesfully parsed.
  map<string, headerInfoStruct>::iterator iter;
  for (iter = v.headerInfoFields.begin(); iter != v.headerInfoFields.end(); iter++) {
    if ( !(iter->second).success) {
      cerr << "Error parsing the info header lines." << endl;
      cerr << "Failed to read: " << (iter->first) << endl;
      cerr << endl;
      cerr << "Info header lines must conform to the spec to allow file validation." << endl;
      exit(1);
    }
  }

  for (iter = v.headerFormatFields.begin(); iter != v.headerFormatFields.end(); iter++) {
    if ( !(iter->second).success) {
      cerr << "Error parsing the format header lines." << endl;
      cerr << "Failed to read: " << (iter->first) << endl;
      cerr << endl;
      cerr << "Info header lines must conform to the spec to allow file validation." << endl;
      exit(1);
    }
  }

  int previousPosition = 0;
  string previousReferenceSequence = "";
  bool positionSorted = true;
  bool referenceSequenceSorted = true;
  map<string, bool> parsedReferenceSequences;

// Read through all the entries in the file.
  while(v.getRecord(currentReferenceSequence)) {

// Check that the current record is not before the previous one (i.e.
// check that the vcf file is sorted).
    if (previousReferenceSequence != v.referenceSequence) {
      if (parsedReferenceSequences[v.referenceSequence]) {referenceSequenceSorted = false;}
      previousReferenceSequence = v.referenceSequence;
    }
    else if (v.position < previousPosition) {positionSorted = false;}

    previousPosition = v.position;
    parsedReferenceSequences[v.referenceSequence] = true;

// For each field in the info string, check that there exists a
// line in the header explaining the field and that the field
// contains the correct number and type of entries.
    for (map<string, string>::iterator iter = v.infoTags.begin(); iter != v.infoTags.end(); iter++) {
      string tag = (*iter).first;
      information sInfo = v.getInfo(tag);
    }

// Parse all of the genotype information and check that it is complete and
// consistent with the header information.
    for (vector<string>::iterator iter = v.genotypes.begin(); iter != v.genotypes.end(); iter++) {
      v.processGenotypeFields(*iter);
      for (vector<string>::iterator formatIter = v.genotypeFormat.begin(); formatIter != v.genotypeFormat.end(); formatIter++) {
        information sInfo = v.getGenotypeInfo(*formatIter);
      }
    }
  }

// Close the vcf file and return.
  v.closeVcf();

// Write to screen that the file is unsorted or that no errors were found.
  if (!referenceSequenceSorted) {cerr << "vcf file has unsorted referenceSequences" << endl;}
  if (!positionSorted) {cerr << "vcf file has unsorted positions" << endl;}
  if (positionSorted && referenceSequenceSorted) {cout << "No errors found with the vcf file." << endl;}

  return 0;
}
