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
  distanceDistribution     = false;
  allowMismatch            = false;
  passFilters              = false;
  findCommon               = false;
  findUnion                = false;
  findUnique               = false;
  sitesOnly                = false;
  processComplex           = false;
  processIndels            = false;
  processMnps              = false;
  processSnps              = false;
  whollyWithin             = false;
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
  cout << "  -d, --distance" << endl;
  cout << "	calculate the distribution of distances of each variant from the nearest bed interval." << endl;
  cout << "  -c, --common" << endl;
  cout << "	output variants present in both files." << endl;
  cout << "  -m, --mismatch" << endl;
  cout << "     variants of the same type at the same locus are considered a match." << endl;
  cout << "  -q, --unique" << endl;
  cout << "	output variants unique to one of the files." << endl;
  cout << "  -u, --union" << endl;
  cout << "	output variants present in either file." << endl;
  cout << "  -s, --sites-only" << endl;
  cout << "	only compare files based on sites.  Do not evaluate the alleles." << endl;
  cout << "  -p, --pass-filters" << endl;
  cout << "	Only variants that pass filters are considered." << endl;
  cout << "  -w, --wholly-within-interval" << endl;
  cout << "	For bed-intersections, start and end of ref/variant allele must fall within interval." << endl;
  cout << "  -1, --snps" << endl;
  cout << "	analyse SNPs." << endl;
  cout << "  -2, --mnps" << endl;
  cout << "	analyse MNPs." << endl;
  cout << "  -3, --indels" << endl;
  cout << "	analyse indels." << endl;
  cout << "  -4, --complex" << endl;
  cout << "	analyse complex events." << endl;
  cout << endl;
  cout << "Additional information:" << endl;
  cout << "  The -c, -u and -q options require either 'a', 'b' or 'q' (not valid for -q) as an argument." << endl;
  cout << endl;
  cout << "  a: Write out records from the first file." << endl;
  cout << "  b: Write out records from the second file." << endl;
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
    {"common", required_argument, 0, 'c'},
    {"in", required_argument, 0, 'i'},
    {"out", required_argument, 0, 'o'},
    {"distance", no_argument, 0, 'd'},
    {"mismatch", no_argument, 0, 'm'},
    {"pass-filters", no_argument, 0, 'p'},
    {"unique", required_argument, 0, 'q'},
    {"union", required_argument, 0, 'u'},
    {"sites-only", no_argument, 0, 's'},
    {"wholly-within-interval", no_argument, 0, 'w'},
    {"snps", no_argument, 0, '1'},
    {"mnps", no_argument, 0, '2'},
    {"indels", no_argument, 0, '3'},
    {"complex", no_argument, 0, '4'},

    {0, 0, 0, 0}
  };

  while (true) {
    int option_index = 0;
    argument = getopt_long(argc, argv, "hb:i:o:dmpc:u:q:sw1234", long_options, &option_index);

    if (argument == -1) {break;}
    switch (argument) {
      // Input bed file.
      case 'b':
        bedFile = optarg;
        break;

      // Determine if the distribution of the variant distance from the
      // nearest bed interval is required.
      case 'd':
        distanceDistribution = true;
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

      // Determine if comparing for exact matches or not.
      case 'm':
        allowMismatch = true;
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

      // An allele must fall wholly within bed interval.
      case 'w':
        whollyWithin = true;
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

      // Analyse complex events.
      case '4':
        processComplex = true;
        break;
      
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
    else {
      cerr << "The file from which the records are to be written needs to be selected." << endl;
      cerr << "    a - write records from the first inputted file." << endl;
      cerr << "    b - write records from the second inputted file." << endl;
      exit(1);
    }
  }
    
  return 0;
}

// Run the tool.
int intersectTool::Run(int argc, char* argv[]) {
  int getOptions = intersectTool::parseCommandLine(argc, argv);

  // Define an output object and open the output file.
  output ofile;
  ofile.outputStream = ofile.openOutputFile(outputFile);

  intersect ints; // Define an intersection object.
  ints.setBooleanFlags(findCommon, findUnion, findUnique, sitesOnly, false, whollyWithin);  // Set the flags required for performing intersections.
  if (writeFrom == "a") {ints.flags.writeFromFirst = true;}
  else if (writeFrom == "b") {ints.flags.writeFromFirst = false;}

  // If intersection is between a vcf file and a bed file, create a vcf and a bed object
  // and intersect.
  if (bedFile != "") {
    vcf v; // Create a vcf object.
    variant var; // Create a variant object.
    var.determineVariantsToProcess(processSnps, processMnps, processIndels, processComplex, false, true, false);

    bed b; // Create a bed object.
    bedStructure bs; // Create a bed structure.

    v.openVcf(vcfFiles[0]);
    b.openBed(bedFile);

    // Parse the headers.
    v.parseHeader(var.headerInfoFields, var.headerFormatFields, var.samples);

    b.parseHeader();

    // Write the header to the output file.
    string taskDescription = "##vcfCtools=intersect " + vcfFiles[0] + ", " + bedFile;
    writeHeader(ofile.outputStream, v, false, taskDescription);

    // Intersect the files.
    ints.intersectVcfBed(v, var, b, bs, ofile);

    // Check that the input files had the same list of reference sequences.
    // If not, it is possible that there were some problems.
    //checkReferenceSequences(v.referenceSequenceVector, b.referenceSequenceVector); // tools.cpp

    // Close the vcf file and return.
    v.closeVcf();

    // Output some brief statistics on the targets.
    cerr << "Number: " << b.numberTargets << endl;
    cerr << "Total target length: " << b.targetLength << endl;
    cerr << "Average target length: " << b.targetLength / b.numberTargets << endl;
    cerr << endl;

    // Print out the distribution of distances to targets.
    cerr << "Histogram describing the distribution of variant distances from targets:" << endl;
    cerr << endl;

    map<string, map<int, unsigned int> >::iterator rIter = ints.distanceDist.begin();
    map<int, unsigned int>::iterator dIter;
    for (; rIter != ints.distanceDist.end(); rIter++) {
      cerr << "Reference sequence: " << rIter->first << endl;
      for (dIter = rIter->second.begin(); dIter != rIter->second.end(); dIter++) {
        cerr << dIter->first << "	" << dIter->second << endl;
      }
    }

    // Close the bed object.
    b.closeBed();

  // Intersection of two vcf files.
  } else {
    vcf v1; // Create a vcf object.
    variant var1; // Create a variant object.
    var1.determineVariantsToProcess(processSnps, processMnps, processIndels, processComplex, false, true, true);

    vcf v2; // Create a vcf object.
    variant var2;
    var2.determineVariantsToProcess(processSnps, processMnps, processIndels, processComplex, false, true, true);
    
    // Open the vcf files.
    v1.openVcf(vcfFiles[0]);
    v2.openVcf(vcfFiles[1]);

    // Read in the header information.
    v1.parseHeader(var1.headerInfoFields, var1.headerFormatFields, var1.samples);
    v2.parseHeader(var2.headerInfoFields, var2.headerFormatFields, var2.samples);

    // Check that the header for the two files contain the same samples.
    if (v1.samples != v2.samples) {
      cerr << "vcf files contain different samples (or sample order)." << endl;
      exit(1);
    } else {
      string taskDescription = "##vcfCTools=intersect " + vcfFiles[0] + ", " + vcfFiles[1];
      if (ints.flags.writeFromFirst) {writeHeader(ofile.outputStream, v1, false, taskDescription);} // tools.cpp
      else {writeHeader(ofile.outputStream, v2, false, taskDescription);} // tools.cpp
    }

    // Intersect the two vcf files.
    ints.intersectVcf(v1, var1, v2, var2, ofile);

    // Check that the input files had the same list of reference sequences.
    // If not, it is possible that there were some problems.
    ints.checkReferenceSequences(var1, var2);

    // Close the vcf files.
    v1.closeVcf();
    v2.closeVcf();
  }

  return 0;
}
