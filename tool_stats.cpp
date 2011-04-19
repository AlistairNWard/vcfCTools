// ******************************************************
// vcfCTools (c) 2011 Alistair Ward
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ------------------------------------------------------
// Last modified: 18 February 2011
// ------------------------------------------------------
// Generate statistcs on an input vcf file.
// ******************************************************

#include "tool_stats.h"

using namespace std;
using namespace vcfCTools;

// statsTool imlementation.
statsTool::statsTool(void)
  : AbstractTool()
{
  recordsInMemory = 100;
  generateAfs = false;
  groupVariants = false;
}

// Destructor.
statsTool::~statsTool(void) {}

// Help
int statsTool::Help(void) {
  cout << "Stats help" << endl;
  return 0;
}

// Parse the command line and get all required and optional arguments.
int statsTool::parseCommandLine(int argc, char* argv[]) {
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
    {"out", required_argument, 0, 'o'},
    {"allele-frequency-spectrum", no_argument, 0, 'a'},
    {"group-variants", required_argument, 0, 'g'},

    {0, 0, 0, 0}
  };

  while (true) {
    int option_index = 0;
    argument = getopt_long(argc, argv, "hi:o:ag:", long_options, &option_index);

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

      // Output file.
      case 'o':
        outputFile = optarg;
        break;
 
      // Generate the allele frequency spectrum.
      case 'a':
        generateAfs = true;
        break;
      // Group variants before generating stats.
      case 'g':
        groupVariants = true;
        referenceFasta = optarg;
        break;

      //
      case '?':
        cout << "Unknown option: " << argv[optind - 2] << endl;
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
int statsTool::Run(int argc, char* argv[]) {
  int getOptions = statsTool::parseCommandLine(argc, argv);

  vcf v; // Create a vcf object.
  statistics stats; // Create a statistics object.

  v.openVcf(vcfFile);
  output = openOutputFile(outputFile);
  v.parseHeader();

// If the vcf records are to be concatanated into groups, a variantGroup
// structure is required.
  //variantGroup vc;

// Read through all the entries in the file.  First construct the
// structure to contain the variants in memory and populate.
  v.success = v.getRecord();
  string currentReferenceSequence = v.variantRecord.referenceSequence;
  v.success = v.buildVariantStructure(recordsInMemory, currentReferenceSequence, false, output);
  while (v.success) {
    if (v.variantRecord.referenceSequence == currentReferenceSequence) {

      // Add the new variant to the structure and process the first element,
      // finally erasing this element.  The structure should consequently
      // maintain the same size.
      v.addVariantToStructure();
      v.variantsIter = v.variants.begin();
      stats.generateStatistics(v, v.variantsIter->first, v.variantsIter->second, v.variantsInformation[v.variantsIter->first], generateAfs);
      v.variants.erase(v.variantsIter);
      v.success = v.getRecord();
    } else {

      // Process all of the elements remaining in the variants structure (i.e.
      // those for the current reference sequence), deleting them as they are
      // processed, then set the currentReferenceSequence to the new value and
      // rebuild the variant structure.
      for (v.variantsIter = v.variants.begin(); v.variantsIter != v.variants.end(); v.variantsIter++) {
        stats.generateStatistics(v, v.variantsIter->first, v.variantsIter->second, v.variantsInformation[v.variantsIter->first], generateAfs);
        v.variants.erase(v.variantsIter);
      }
      currentReferenceSequence = v.variantRecord.referenceSequence;
      v.success = v.buildVariantStructure(recordsInMemory, currentReferenceSequence, false, output);
    }
  }

// When the end of the vcf file is reached, process the remaining records stored
// in the variants structure.
  for (v.variantsIter = v.variants.begin(); v.variantsIter != v.variants.end(); v.variantsIter++) {
    stats.generateStatistics(v, v.variantsIter->first, v.variantsIter->second, v.variantsInformation[v.variantsIter->first], generateAfs);
    v.variants.erase(v.variantsIter);
  }
  //unsigned int noGroups = 0;
  //while(v.getRecord()) {
  //  if (groupVariants) {
  //    success = v.getVariantGroup(vc, referenceFasta);
  //    //cout << vc.start << ", Number of records=" << vc.noRecords << ", Number of alternates=" << vc.noAlts << endl;
  //  } else {
  //    stats.generateStatistics(v, generateAfs);
  //  }
  //}
  //if (groupVariants) {cout << "No groups: " << vc.noGroups << endl;}

// Count the total number of variants in each class and then rint out the
// statistics.
  stats.countByFilter();
  if (stats.hasSnp) {
    stats.printSnpStatistics(output);
    if (stats.hasAnnotations) {stats.printSnpAnnotations(output);}
    if (generateAfs) {
      stats.printAcs(output);
      stats.printAfs(output);
    }
  }
  if (stats.hasMnp) {stats.printMnpStatistics(output);}
  if (stats.hasIndel) {stats.printIndelStatistics(output);}

// Close the vcf file and return.
  v.closeVcf();

  return 0;
}
