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
{}

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
  cout << "	input vcf file," << endl;
  cout << "  -o, --output" << endl;
  cout << "	output vcf file," << endl;
  cout << "  -d, --dbsnp" << endl;
  cout << "	input dbsnp vcf file," << endl;
  cout << "  -h, --hapmap" << endl;
  cout << "	input hapmap vcf file," << endl;
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
    {"hapmap", required_argument, 0, 'm'},

    {0, 0, 0, 0}
  };

    int option_index = 0;
    argument = getopt_long(argc, argv, "hi:o:d:m:", long_options, &option_index);

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
        annotationFile = optarg;
        annotateDbsnp = true;
        break;

      // Input hapmap vcf file.
      case 'm':
        annotationFile = optarg;
        annotateHapmap = true;
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

// Check that a vcf file was specified.
 if (vcfFile == "") {
    cerr << "A vcf file must be specified (--in, -i)." << endl;
    exit(1);
  }

// Check that a either a dbsnp or a hapmap vcf file was specified.
 if ( (annotateDbsnp && annotateHapmap) || (!annotateDbsnp && !annotateHapmap) ) {
    cerr << "Either a dbsnp or a hapmap vcf file must be specified (--dbsnp, -d, --hapmap, -h)." << endl;
    exit(1);
  }

// Check that a dbsnp or hapmap file have been

  return 0;
}

// Intersect two vcf files.  It is assumed that the two files are
// sorted by genomic coordinates and the reference sequences are
// in the same order.
void annotateTool::annotateVcf(vcf& v, vcf& d, ostream* output) {
  bool success1 = v.getRecord();
  bool success2 = d.getRecord();
  string currentReferenceSequence = v.referenceSequence;

// Finish when the end of the first file has been reached.
  while (success1) {

// If the end of the dbsnp vcf file is reached, write out the
// remaining records from the vcf file.
    if (!success2) {
      *output << v.record << endl;
      success1 = v.getRecord();
    }

    if (v.referenceSequence == d.referenceSequence && v.referenceSequence == currentReferenceSequence) {
      if (v.position == d.position) {
        if (d.dbsnpVcf) {
          int number;
          string type;
          vector<string> values;
          string tag = "VC";
          d.getInfo(tag, number, type, values);
          if (values[0] == "SNP") {v.rsid = d.rsid;}
        }
        else if (d.hapmapVcf) {v.info += ";HM3";}
        string record = v.buildRecord(true);
        *output << record << endl;

        success1 = v.getRecord();
        success2 = d.getRecord();
      }
      else if (d.position > v.position) {success1 = v.parseVcf(d.referenceSequence, d.position, true, output);}
      else if (v.position > d.position) {success2 = d.parseVcf(v.referenceSequence, v.position, false, NULL);}
    }
    else {
      if (v.referenceSequence == currentReferenceSequence) {success1 = v.parseVcf(d.referenceSequence, d.position, true, output);}
      else if (d.referenceSequence == currentReferenceSequence) {success2 = d.parseVcf(v.referenceSequence, v.position, false, NULL);}

// If the last record for a reference sequence is the same for both vcf
// files, they will both have referenceSequences different from the
// current reference sequence.  Change the reference sequence to reflect
// this and proceed.
      else {
        if (v.referenceSequence != d.referenceSequence) {
          cerr << "ERROR: Reference sequences for both files are unexpectedly different." << endl;
          cerr << "Check that both files contain records for the following reference sequences:" << endl;
          cerr << "\t" << v.referenceSequence << " and " << d.referenceSequence << endl;
          exit(1);
        }
      }
      currentReferenceSequence = v.referenceSequence;
    }
  }
}

// Run the tool.
int annotateTool::Run(int argc, char* argv[]) {
  annotateDbsnp = false;
  annotateHapmap = false;
  int getOptions = annotateTool::parseCommandLine(argc, argv);
  output = openOutputFile(outputFile);

  vcf v; // Define vcf object.
  vcf d; // Define dbsnp/hapmap vcf object.

  if (annotateDbsnp) {
    d.dbsnpVcf = true;
    d.processInfo = true;
  }
  else if (annotateHapmap) {d.hapmapVcf = true;}

// Open the vcf files.
  v.openVcf(vcfFile);
  d.openVcf(annotationFile);

// Read in the header information.
  v.parseHeader();
  d.parseHeader();

// Add an extra line to the vcf header to indicate the file used for
// performing dbsnp annotation.
  string taskDescription = "##vcfCTools=annotated vcf file with ";
  if (d.dbsnpVcf) {taskDescription += "dbSNP file " + annotationFile;}
  else if (d.hapmapVcf) {
    taskDescription += "hapmap file " + annotationFile;
    v.headerInfoLine["HM3"] = "##INFO=<ID=HM3,Number=0,Type=Flag,Description=\"Hapmap3.2 membership determined from file " + \
                                annotationFile + "\">";
    v.headerInfoLine["HM3A"] = "##INFO=<ID=HM3A,Number=0,Type=Flag,Description=\"Hapmap3.2 membership (with different alleles) \
                                , determined from file " + annotationFile + "\">";
  }
  writeHeader(output, v, true, taskDescription); // tools.cpp

// Annotate the vcf file.
  annotateVcf(v, d, output);

// Check that the input files had the same list of reference sequences.
// If not, it is possible that there were some problems.
  checkReferenceSequences(v.referenceSequenceVector, d.referenceSequenceVector); // tools.cpp

// Close the vcf files.
  v.closeVcf();
  d.closeVcf();

  return 0;
}
