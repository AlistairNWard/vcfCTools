// ******************************************************
// vcfCTools (c) 2011 Alistair Ward
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ------------------------------------------------------
// Last modified: 18 February 2011
// ------------------------------------------------------
// Calculate the intersection of two vcf files or a vcf
// file and a bed file.
// ******************************************************

#include "tool_intersect.h"

using namespace std;
using namespace vcfCTools;

// intersectTool imlementation.
intersectTool::intersectTool(void)
  : AbstractTool()
{
  recordsInMemory = 100;
  passFilters     = false;
  findCommon      = false;
  findUnion       = false;
  findUnique      = false;
  sitesOnly       = false;
  currentReferenceSequence = "";
}

// Destructor.
intersectTool::~intersectTool(void) {}

// Help
int intersectTool::Help(void) {
  cout << "Intersect help" << endl;
  cout << "Usage: ./vcfCTools intersect [options]." << endl;
  cout << endl;
  cout << "Options:" << endl;
  cout << "  -h, --help" << endl;
  cout << "	display intersect help." << endl;
  cout << "  -i, --in" << endl;
  cout << "	input vcf files (two, or one if intersecting with bed file)." << endl;
  cout << "  -b, --bed" << endl;
  cout << "	input bed file." << endl;
  cout << "  -o, --output" << endl;
  cout << "	output vcf file." << endl;
  cout << "  -c, --common" << endl;
  cout << "	output variants present in both files." << endl;
  cout << "  -u, --union" << endl;
  cout << "	output variants present in either file." << endl;
  cout << "  -q, --unique" << endl;
  cout << "	output variants unique to one of the files." << endl;
  cout << "  -s, --sites-only" << endl;
  cout << "	only compare files based on sites.  Do not evaluate the alleles." << endl;
  cout << "  -p, --pass-filters" << endl;
  cout << "	Only variants that pass filters are considered." << endl;
  cout << "  -1, --snps" << endl;
  cout << "	analyse SNPs." << endl;
  cout << "  -2, --mnps" << endl;
  cout << "	analyse MNPs." << endl;
  cout << "  -3, --indels" << endl;
  cout << "	analyse indels." << endl;
  cout << endl;
  cout << "Additional information:" << endl;
  cout << "  The -c, -u and -q options require either 'a', 'b' or 'q' (not valid for -q) as an argument." << endl;
  cout << endl;
  cout << "  a: Write out records from the first file." << endl;
  cout << "  b: Write out records from the second file." << endl;
  cout << "  c: Only write out variants if they share ref and alt alleles (info from first file)." << endl;
  cout << "  d: Only write out variants if they share ref and alt alleles (info from second file)." << endl;
  cout << "  u: Write out all records from both files." << endl;
  cout << "  q: Write out records with the highest variant quality." << endl;
  cout << endl;
  exit(0);

  return 0;
}

// Parse the command line and get all required and optional arguments.
int intersectTool::parseCommandLine(int argc, char* argv[]) {
  commandLine = argv[0];
  for (int i = 2; i < argc; i++) {
    commandLine += " ";
    commandLine += argv[i];
  }

  int argument; // Counter for getopt.

  // Define the long options.
  static struct option long_options[] = {
    {"help", no_argument, 0, 'h'},
    {"bed", required_argument, 0, 'b'},
    {"in", required_argument, 0, 'i'},
    {"out", required_argument, 0, 'o'},
    {"pass-filters", no_argument, 0, 'p'},
    {"common", required_argument, 0, 'c'},
    {"union", required_argument, 0, 'u'},
    {"unique", required_argument, 0, 'q'},
    {"sites-only", no_argument, 0, 's'},
    {"snps", no_argument, 0, '1'},
    {"mnps", no_argument, 0, '2'},
    {"indels", no_argument, 0, '3'},

    {0, 0, 0, 0}
  };

  while (true) {
    int option_index = 0;
    argument = getopt_long(argc, argv, "hb:i:o:pc:u:q:s123", long_options, &option_index);

    if (argument == -1) {break;}
    switch (argument) {
      // Input bed file.
      case 'b':
        bedFile = optarg;
        break;

      // Input vcf file - required input.
      case 'i':
        vcfFiles.push_back(optarg);
        break;

      // Help.
      case 'h':
        return Help();

      // Output file.
      case 'o':
        outputFile = optarg;
        break;

      // Common variants.
      case 'c':
        findCommon = true;
        writeFrom = optarg;
        break;

      // Find the union.
      case 'u':
        findUnion = true;
        writeFrom = optarg;
        break;

      // Unique variants.
      case 'q':
        findUnique = true;
        writeFrom = optarg;
        break;

      // Only compare variants based on the position.  Do not
      // interrogate the alleles.
      case 's':
        sitesOnly = true;
        break;

      // Only consider variants if they pass filters.
      case 'p':
        passFilters = true;
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

// Check that either two vcf files or one vcf and one bed file is specified.
  if (vcfFiles.size() == 0 || (vcfFiles.size() == 1 and bedFile == "") || (bedFile != "" && vcfFiles.size() != 1) || vcfFiles.size() > 2) {
    cerr << "Two vcf files or a vcf and a bed file must be specified (--in, -i, --bed, -b)." << endl;
    exit(1);
  }

// Check whether finding common, union or unique variants.  If more than one of
// these options are selected, terminate the program.
  if ( (findCommon + findUnion + findUnique) > 1) {
    cerr << "Only one operation (-c [--common], -u [--union] or -q [--unique]) can be performed at once." << endl;
    exit(1);
  } else if (!findCommon && !findUnion && !findUnique) {
    cerr << "One operation (-c [--common], -u [--union] or -q [--unique]) must be selected." << endl;
    exit(1);
  } else {
    if (writeFrom == "a") {cerr << "Writing out records from file: " << vcfFiles[0] << endl;}
    else if (writeFrom == "b") {cerr << "Writing out records from file: " << vcfFiles[1] << endl;}
    else if (writeFrom == "c" && findUnion || findCommon) {cerr << "Writing out records from file: " << vcfFiles[0] << endl;}
    else if (writeFrom == "d" && findUnion || findCommon) {cerr << "Writing out records from file: " << vcfFiles[1] << endl;}
    else if (writeFrom == "u" && findUnion || findCommon) {cerr << "Writing out all records." << endl;}
    else if ( (findCommon || findUnion) && writeFrom == "q") {cerr << "Writing out records with the highest quality value." << endl;}
    else {
      cerr << "The file from which the records are to be written needs to be selected." << endl;
      cerr << "    a - write records from the first inputted file." << endl;
      cerr << "    b - write records from the second inputted file." << endl;
      cout << "    c - write out variants if they share ref and alt alleles (info from first file)." << endl;
      cerr << "    d - write out variants if they share ref and alt alleles (info from second file)." << endl;
      cerr << "    u - write out all records from both files." << endl;
      if (findCommon) {cerr << "    q - write records with the highest variant quality." << endl;}
      exit(1);
    }
  }
    
  return 0;
}

// Run the tool.
int intersectTool::Run(int argc, char* argv[]) {
  int getOptions = intersectTool::parseCommandLine(argc, argv);
  output = openOutputFile(outputFile);

  intersect ints; // Define an intersection object.
  ints.setBooleanFlags(findCommon, findUnion, findUnique, sitesOnly, false);  // Set the flags required for performing intersections.
  ints.writeFrom = writeFrom;

// If intersection is between a vcf file and a bed file, create a vcf and a bed object
// and intersect.
  if (bedFile != "") {
    cerr << "DISABLED UNTIL clearReferenceSequence and variant structure comparisons" << endl;
    cerr << "have been updated to handle bed files." << endl;
    exit(0);
    vcf v; // Create a vcf object.
    variant var; // Create a variant object.
    var.determineVariantsToProcess(processSnps, processMnps, processIndels);

    bed b; // Create a bed object.
    bedStructure bs; // Create a bed structure.

    v.openVcf(vcfFiles[0]);
    b.openBed(bedFile);

// Parse the headers.
    v.parseHeader();
    b.parseHeader();

// Write the header to the output file.
    string taskDescription = "##vcfCtools=intersect " + vcfFiles[0] + ", " + bedFile;
    writeHeader(output, v, false, taskDescription);

// Intersect the files.
    ints.intersectVcfBed(v, var, b, bs, findUnique, false, output);

// Check that the input files had the same list of reference sequences.
// If not, it is possible that there were some problems.
    checkReferenceSequences(v.referenceSequenceVector, b.referenceSequenceVector); // tools.cpp

// Close the vcf file and return.
    v.closeVcf();

// Output some brief statistics on the targets.
    cout << "Number: " << b.numberTargets << endl;
    cout << "Total target length: " << b.targetLength << endl;
    cout << "Average target length: " << b.targetLength / b.numberTargets << endl;
    
// Close the bed object.
    b.closeBed();

// Intersection of two vcf files.
  } else {
    vcf v1; // Create a vcf object.
    variant var1; // Create a variant object.
    var1.determineVariantsToProcess(processSnps, processMnps, processIndels);

    vcf v2; // Create a vcf object.
    variant var2;
    var2.determineVariantsToProcess(processSnps, processMnps, processIndels);
    
// Open the vcf files.
    v1.openVcf(vcfFiles[0]);
    v2.openVcf(vcfFiles[1]);

// Read in the header information.
    v1.parseHeader();
    v2.parseHeader();
    checkDataSets(v1, v2); // tools.cpp

// Check that the header for the two files contain the same samples.
    if (v1.samples != v2.samples) {
      cerr << "vcf files contain different samples (or sample order)." << endl;
      exit(1);
    } else {
      string taskDescription = "##vcfCTools=intersect " + vcfFiles[0] + ", " + vcfFiles[1];
      writeHeader(output, v1, false, taskDescription); // tools.cpp
    }

// Intersect the two vcf files.
    ints.intersectVcf(v1, var1, v2, var2, output);

// Check that the input files had the same list of reference sequences.
// If not, it is possible that there were some problems.
    checkReferenceSequences(v1.referenceSequenceVector, v2.referenceSequenceVector); // tools.cpp

// Close the vcf files.
    v1.closeVcf();
    v2.closeVcf();
  }

  return 0;
}
