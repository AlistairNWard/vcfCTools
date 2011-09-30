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

using namespace std;
using namespace vcfCTools;

// validateTool imlementation.
validateTool::validateTool(void)
  : AbstractTool()
{
  currentReferenceSequence = "";
  error                    = false;
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

      // Help.
      case 'h':
        return Help();

      //
      case '?':
        cerr << "Unknown option: " << argv[optind - 1] << endl;
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

  // Create a vcf object.
  vcf v; // Create a vcf object.
  v.openVcf(vcfFile);

  // Define a variant object.
  variant var; // Define variant object.
  var.determineVariantsToProcess(true, true, true, true, true, true, false, true, false);

  // Define a header object and parse the header information.
  vcfHeader header;
  header.parseHeader(v.input);

  // Check that all of the info descriptions in the header are in the correct form.
  map<string, headerInfo>::iterator iter;
  for (iter = header.infoFields.begin(); iter != header.infoFields.end(); iter++) {
    if ( !(iter->second.success) ) {
      cerr << "ERROR: Malformed info string in the header: " << iter->first << endl;
      exit(1);
    }
  }

  for (iter = header.formatFields.begin(); iter != header.formatFields.end(); iter++) {
    if ( !(iter->second.success) ) {
      cerr << "ERROR: Malformed format string in the header: " << iter->first << endl;
      exit(1);
    }
  }

  // Read through all the entries in the file.
  v.success = v.getRecord();
  while (v.success) {

    // Build the variant structure for this reference sequence.
    if (var.originalVariantsMap.size() == 0) {
      currentReferenceSequence = v.variantRecord.referenceSequence;
      v.success = var.buildVariantStructure(v);
    }

    // Loop over the variant structure until it is empty.  While v.update is true,
    // i.e. when the reference sequence is still the current reference sequence,
    // keep adding variants to the structure.
    while (var.originalVariantsMap.size() != 0) {
      if (v.variantRecord.referenceSequence == currentReferenceSequence && v.success) {
        var.addVariantToStructure(v.position, v.variantRecord);
        v.success = v.getRecord();
      }
      var.ovmIter = var.originalVariantsMap.begin();

      // Loop over all records at this locus.
      var.ovIter = var.ovmIter->second.begin();
      for (; var.ovIter != var.ovmIter->second.end(); var.ovIter++) {

        // Check the info string for inconsistencies.
        variantInfo info(var.ovIter->info);
        info.validateInfo(header, var.ovIter->referenceSequence, var.ovIter->position, var.ovIter->numberAlts, error);

        // Check the genotypes for inconsistencies.
        if (var.ovIter->hasGenotypes) {
          genotypeInfo gen(var.ovIter->genotypeFormat, var.ovIter->genotypes);
          gen.validateGenotypes(header, var.ovIter->referenceSequence, var.ovIter->position, var.ovIter->numberAlts, var.samples, error);
        }
      }
      var.originalVariantsMap.erase(var.ovmIter);
    }
  }

  // Close the vcf files.
  v.closeVcf();

  // If no errors were found, indicate that this was the case.
  if (!error) {cerr << "No errors found with vcf file." << endl;}

  return 0;
}
