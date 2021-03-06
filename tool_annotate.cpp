// ******************************************************
// vcfCTools (c) 2011 Alistair Ward
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ------------------------------------------------------
// Last modified: 18 February 2011
// ------------------------------------------------------
// Annotate a vcf file with dbsnp or hapmap membership.
// The input dbsnp or hapmap files need to be in vcf
// format.
// ******************************************************

#include "tool_annotate.h"

using namespace std;
using namespace vcfCTools;

// intersectTool imlementation.
annotateTool::annotateTool(void)
  : AbstractTool()
{
  currentReferenceSequence = "";
  annotateDbsnp            = false;
  annotateVcf              = false;
  annotateBed              = false;
  processComplex           = false;
  processIndels            = false;
  processMnps              = false;
  processRearrangements    = false;
  processSnps              = false;
  processSvs               = false;
  sitesOnly                = false;
  whollyWithin             = false;
}

// Destructor.
annotateTool::~annotateTool(void) {}

// Help
int annotateTool::Help(void) {
  cout << "Annotate help" << endl;
  cout << "Usage: ./vcfCTools annotate [options]." << endl;
  cout << endl;
  cout << "Options:" << endl;
  cout << "  -h, --help" << endl;
  cout << "	display intersect help." << endl;
  cout << "  -i, --in" << endl;
  cout << "	input vcf file." << endl;
  cout << "  -o, --output" << endl;
  cout << "	output vcf file." << endl;
  cout << "  -a, --annotation-vcf" << endl;
  cout << "	input annotation vcf file." << endl;
  cout << "  -d, --dbsnp" << endl;
  cout << "	input dbsnp vcf file." << endl;
  cout << "  -b, --bed" << endl;
  cout << "	input bed file." << endl;
  cout << "  -s, --sites-only" << endl;
  cout << "	only compare files based on sites.  Do not evaluate the alleles." << endl;
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
  cout << "  -5, --structural-variants" << endl;
  cout << "	analyse structural variantion events." << endl;
  cout << "  -6, --rearrangements" << endl;
  cout << "	analyse complex rearrangement events." << endl;
  cout << endl;

  return 0;
}

// Parse the command line and get all required and optional arguments.
int annotateTool::parseCommandLine(int argc, char* argv[]) {
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
    {"dbsnp", required_argument, 0, 'd'},
    {"annotation-vcf", required_argument, 0, 'a'},
    {"bed", required_argument, 0, 'b'},
    {"sites-only", no_argument, 0, 'b'},
    {"wholly-within", no_argument, 0, 'b'},
    {"snps", no_argument, 0, '1'},
    {"mnps", no_argument, 0, '2'},
    {"indels", no_argument, 0, '3'},
    {"complex", no_argument, 0, '4'},
    {"structural-variants", no_argument, 0, '5'},
    {"rearrangements", no_argument, 0, '6'},

    {0, 0, 0, 0}
  };

    int option_index = 0;
    argument = getopt_long(argc, argv, "hi:o:a:b:d:w:s123456", long_options, &option_index);

    if (argument == -1) {break;}
    switch (argument) {

      // Input vcf file - required input.
      case 'i':
        vcfFile = optarg;
        break;

      // Ouput vcf file.
      case 'o':
        outputFile = optarg;
        break;

      // Input dbsnp vcf file.
      case 'd':
        annVcfFile = optarg;
        annotateDbsnp = true;
        break;

      // Input hapmap vcf file.
      case 'a':
        annVcfFile = optarg;
        annotateVcf = true;
        break;

      // Input bed file.
      case 'b':
        bedFile = optarg;
        annotateBed = true;
        break;

      // Only compare variants based on the position.  Do not
      // interrogate the alleles.
        case 's':
        sitesOnly = true;
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

      // Analyse indels.
      case '4':
        processComplex = true;
        break;

      // Analyse structural variants.
      case '5':
        processSvs = true;
        break;

      // Analyse complex rearrangements.
      case '6':
        processRearrangements = true;
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

// Check that a either a dbsnp, hapmap vcf and/or bed file was specified.
 if ( !annotateDbsnp && !annotateVcf && !annotateBed ) {
    cerr << "A dbsnp, annotation vcf file(s) or a bed file must be specified (-d, --dbsnp, -a, --annotation-vcf, -b, --bed)." << endl;
    exit(1);
  }

// Check that a dbsnp or hapmap file have been

  return 0;
}

// Run the tool.
int annotateTool::Run(int argc, char* argv[]) {
  int getOptions = annotateTool::parseCommandLine(argc, argv);

  // Define an output object and open the output file.
  output ofile;
  ofile.outputStream = ofile.openOutputFile(outputFile);

  // Define the vcf object.
  vcf v;
  v.openVcf(vcfFile);

  // Define the variant object.
  variant var;
  var.determineVariantsToProcess(processSnps, processMnps, processIndels, processComplex, processSvs, processRearrangements, false, true, true);

  intersect ints; // Define an intersection object.
  ints.setBooleanFlags(true, false, false, sitesOnly, true, whollyWithin);  // Set the flags required for performing intersections.
  ints.flags.writeFromFirst = true;

  // Define the header object and read in the header information.
  vcfHeader header;
  header.parseHeader(v.input);

  // Add an extra line to the vcf header to indicate the file used for
  // performing dbsnp annotation.
  string taskDescription = "##vcfCTools=annotated vcf file with ";

  if (annotateDbsnp || annotateVcf) {
    if (annotateDbsnp) {
      header.infoLine["dbSNP"] = "##INFO=<ID=dbSNP,Number=0,Type=Flag,Description=\"Membership in dbSNP file " + annVcfFile;
      header.infoLine["dbSNP"] += " with common alleles.\">";
      //v.headerInfoLine["dbSNPX"] = "##INFO=<ID=dbSNPX,Number=0,Type=Flag,Description=\"Membership in dbSNP file " + annVcfFile;
      //v.headerInfoLine["dbSNPX"] += " with different alleles.\">";
      //v.headerInfoLine["dbSNPM"] = "##INFO=<ID=dbSNPM,Number=0,Type=Flag,Description=\"Membership in dbSNP file " + annVcfFile;
      //v.headerInfoLine["dbSNPM"] += ". Either the vcf or dbSNP entry show a variant with multiple alternate alleles.\">";
      taskDescription += "vcf file " + annVcfFile;
    }
    taskDescription += "vcf file " + annVcfFile;
    header.writeHeader(ofile.outputStream, false, taskDescription);

    // Either a vcf file, a dbsnp vcf file or a bed file can be provided for
    // annotation.  To annotate from multiple files, piping should be used.
    vcf annVcf; // Define a vcf object.
    annVcf.openVcf(annVcfFile);

    // Define a variant object
    variant annVar; // Define a variant object
    annVar.determineVariantsToProcess(processSnps, processMnps, processIndels, processComplex, processSvs, processRearrangements, false, true, true);

    // Define a header object and parse the header of the annotation vcf file.
    vcfHeader annHeader;
    annHeader.parseHeader(annVcf.input);
    if (annotateDbsnp) {annVar.isDbsnp = true;}

    // Perform the annotation by intersecting the two vcf files.
    ints.intersectVcf(header, annHeader, v, var, annVcf, annVar, ofile);

    // Check that the input files had the same list of reference sequences.
    // If not, it is possible that there were some problems.
    ints.checkReferenceSequences(var, annVar);

    // Close the vcf files.
    v.closeVcf();
    annVcf.closeVcf();
  }

  if (annotateBed) {
    taskDescription += "bed file " + bedFile;
    header.writeHeader(ofile.outputStream, false, taskDescription);

    bed b; // Define a bed object.
    bedStructure bs; // Define a bed structure object.

    // Open the bed file and parse the header.
    b.openBed(bedFile);
    b.parseHeader(bedFile);

    // Perform the annotation by intersecting the vcf file with the bed file.
    ints.intersectVcfBed(header, v, var, b, bs, ofile);

    //checkReferenceSequences(v.referenceSequenceVector, b.referenceSequenceVector);} // tools.cpp
    // Close the vcf and bed files.
    b.closeBed();
    v.closeVcf();
  }

  return 0;
}
