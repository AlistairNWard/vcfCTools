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
  generateAfs              = false;
  generateDetailed         = false;
  useAnnotations           = false;
  annotationFlagsString    = "";
  currentReferenceSequence = "";
  sampleSnps               = false;
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
  cout << "  -d, --detailed" << endl;
  cout << "     generate detailed statistics for each SNP considering samples with genotype quality greater than value specified.." << endl;
  cout << "  -n, --annotation" << endl;
  cout << "     include statistics on listed annotations (comma separated list or 'all')." << endl;
  cout << "  -s, --sample-snps" << endl;
  cout << "     include SNP statistics on all individual samples (requires genotypes to be present and s cut-off genotype quality to be specified)." << endl;
  cout << "  -1, --snps" << endl;
  cout << "	analyse SNPs." << endl;
  cout << "  -2, --mnps" << endl;
  cout << "	analyse MNPs." << endl;
  cout << "  -3, --indels" << endl;
  cout << "	analyse indels." << endl;
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
    {"annotations", required_argument, 0, 'n'},
    {"detailed", required_argument, 0, 'd'},
    {"sample-snps", required_argument, 0, 's'},
    {"snps", no_argument, 0, '1'},
    {"mnps", no_argument, 0, '2'},
    {"indels", no_argument, 0, '3'},

    {0, 0, 0, 0}
  };

  while (true) {
    int option_index = 0;
    argument = getopt_long(argc, argv, "hi:o:ad:n:s:124", long_options, &option_index);

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

      // Generate detailed SNP statistics.
      case 'd':
        generateDetailed = true;
        detailedGenotypeQualityString = optarg;
        break;

      // Determine whether to output stats on annotations and
      // if so, which flags.
      case 'n':
        useAnnotations = true;
        annotationFlagsString = optarg;
        break;

      // Generate SNP statistics for all samples.
      case 's':
        sampleSnps = true;
        genotypeQualityString = optarg;
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
        cerr << "Unknown option: " << argv[optind - 2] << endl;
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
  var.determineVariantsToProcess(processSnps, processMnps, processIndels, false);
  statistics stats; // Create a statistics object.

  v.openVcf(vcfFile);
  output = openOutputFile(outputFile);
  v.parseHeader();

// If statistics are being generated on a per-sample basis (or detailed
// statistics are being generate, check that genotypes exist.
  if (sampleSnps || generateDetailed) {
    if (sampleSnps) {stats.minGenotypeQuality = atof(genotypeQualityString.c_str());}
    if (generateDetailed) {stats.minDetailedGenotypeQuality = atof(detailedGenotypeQualityString.c_str());}
    if (stats.minGenotypeQuality == 0 && ( genotypeQualityString != "0" && genotypeQualityString != "0." && genotypeQualityString != "0.0") ) {
      cerr << "ERROR: genotype quality for --sample-snps (-s) must be a double." << endl;
      exit(1);
    }
    if (stats.minDetailedGenotypeQuality == 0 && ( detailedGenotypeQualityString != "0" && detailedGenotypeQualityString != "0." && 
        detailedGenotypeQualityString != "0.0") ) {
      cerr << "ERROR: genotype quality for --detailed (-d) must be a double." << endl;
      exit(1);
    }
    if (!v.hasGenotypes) {
      cerr << "ERROR: Genotype information must be present to perform sample level statistics." << endl;
      exit(1);
    } else {
      if (sampleSnps) {stats.processSampleSnps = true;}
      if (generateDetailed) {stats.generateDetailed  = true;}
    }
  }

// If statistics on annotations are required, generate a list of flags to get
// statistics on.  Provide a warning if the flags do not appear in the header.
  if (useAnnotations) {
    size_t found = annotationFlagsString.find(",");
    annotationFlags.clear();
    if (found == string::npos) {annotationFlags.push_back(annotationFlagsString);}
    else {annotationFlags = split(annotationFlagsString, ",");}
    for (vector<string>::iterator iter = annotationFlags.begin(); iter != annotationFlags.end(); iter++) {
      if (*iter != "all" && v.headerInfoFields.count(*iter) == 0) {
        cerr << "WARNING: Info ID " << *iter << " is used for annotation stats, but does not appear in the header." << endl;
      }
    }
  }

// Print the header for detailed statistics if necessary.
  if (generateDetailed) {stats.printDetailedHeader(output);}

// Read through all the entries in the file.  First construct the
// structure to contain the variants in memory and populate.
  while (v.success) {
    // Build the variant structure for this reference sequence.
    if (var.variantMap.size() == 0) {
      currentReferenceSequence = v.variantRecord.referenceSequence;
      v.success = var.buildVariantStructure(v);
    }

    // Loop over the variant structure until it is empty.  While v.update is true,
    // i.e. when the reference sequence is still the current reference sequence,
    // keep adding variants to the structre.
    //while (v.variants.size() != 0) {
    while (var.variantMap.size() != 0) {
      if (v.variantRecord.referenceSequence == currentReferenceSequence && v.success) {
        var.addVariantToStructure(v.position, v.variantRecord, v.dbsnpVcf);
        v.success = v.getRecord(currentReferenceSequence);
      }
      var.vmIter = var.variantMap.begin();
      stats.generateStatistics(var, v, var.vmIter->first, useAnnotations, annotationFlags, generateAfs, output);
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
  if (sampleSnps) {stats.printSampleSnps(v, output);}

// Close the vcf file and return.
  v.closeVcf();

  return 0;
}
