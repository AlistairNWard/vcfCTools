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
  generateAfs = false;
  currentReferenceSequence = "";
}

// Destructor.
statsTool::~statsTool(void) {}

// Help
int statsTool::Help(void) {
  cout << "Stats help" << endl;
  cout << "Usage: ./vcfCTools stats [options]." << endl;
  cout << endl;
  cout << "Options:" << endl;
  cout << "  -h, --help" << endl;
  cout << "     display intersect help." << endl;
  cout << "  -i, --in" << endl;
  cout << "     input vcf file." << endl;
  cout << "  -o, --output" << endl;
  cout << "     output vcf file." << endl;
  cout << "  -a, --allele-frequency-spectrum" << endl;
  cout << "     generate statistics as a function of the AFS." << endl;
  cout << "  -1, --snps" << endl;
  cout << "     analyse SNPs." << endl;
  cout << "  -2, --mnps" << endl;
  cout << "     analyse MNPs." << endl;
  cout << "  -3, --indels" << endl;
  cout << "     analyse indels." << endl;
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
    {"snps", no_argument, 0, '1'},
    {"mnps", no_argument, 0, '2'},
    {"indels", no_argument, 0, '3'},

    {0, 0, 0, 0}
  };

  while (true) {
    int option_index = 0;
    argument = getopt_long(argc, argv, "hi:o:a123", long_options, &option_index);

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
  variant var; // Create a variant structure to hold the variants.
  var.determineVariantsToProcess(processSnps, processMnps, processIndels);
  statistics stats; // Create a statistics object.

  v.openVcf(vcfFile);
  output = openOutputFile(outputFile);
  v.parseHeader();

// Read through all the entries in the file.  First construct the
// structure to contain the variants in memory and populate.
  v.update = true;
  while (v.success) {
    // Build the variant structure for this reference sequence.
    if (var.variantMap.size() == 0) {
      currentReferenceSequence = v.variantRecord.referenceSequence;
      v.success = var.buildVariantStructure(v, currentReferenceSequence, false, output);
    }

    // Loop over the variant structure until it is empty.  While v.update is true,
    // i.e. when the reference sequence is still the current reference sequence,
    // keep adding variants to the structre.
    //while (v.variants.size() != 0) {
    while (var.variantMap.size() != 0) {
      if (v.update && v.success) {
        var.addVariantToStructure(v.position, v.variantRecord);
        v.success = v.getRecord(currentReferenceSequence);
      }
      var.vmIter = var.variantMap.begin();
      stats.generateStatistics(var, v, var.vmIter->first, generateAfs);
      var.variantMap.erase(var.vmIter);
    } 
  }

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
