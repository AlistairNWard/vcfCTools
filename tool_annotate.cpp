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
  cout << "	input vcf file." << endl;
  cout << "  -o, --output" << endl;
  cout << "	output vcf file." << endl;
  cout << "  -a, --annotation-vcf" << endl;
  cout << "	input annotation vcf file." << endl;
  cout << "  -d, --dbsnp" << endl;
  cout << "	input dbsnp vcf file." << endl;
  cout << "  -b, --bed" << endl;
  cout << "	input bed file." << endl;
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

    {0, 0, 0, 0}
  };

    int option_index = 0;
    argument = getopt_long(argc, argv, "hi:o:d:a:b:", long_options, &option_index);

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
        dbsnpFile = optarg;
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

// Check that a either a dbsnp, hapmap vcf and/or bed file was specified.
 if ( !annotateDbsnp && !annotateVcf && !annotateBed ) {
    cerr << "A dbsnp, annotation vcf file(s) or a bed file must be specified (-d, --dbsnp, -a, --annotation-vcf, -b, --bed)." << endl;
    exit(1);
  }

// Check that a dbsnp or hapmap file have been

  return 0;
}

// Intersect two vcf files.  It is assumed that the two files are
// sorted by genomic coordinates and the reference sequences are
// in the same order.
void annotateTool::annotate(vcf& v, vcf& dbsnp, vcf& annVcf, bed& b, bool haveDbsnp, bool haveVcf, bool haveBed, ostream* output) {
  bool successVcf;
  successVcf = v.getRecord();
  currentReferenceSequence = v.referenceSequence;

// Get the next record from the requested annotation files.
  bool successDbsnp = false;
  bool successAnnVcf = false;
  bool successBed = false;
  bool build = false;
  if (haveDbsnp) {successDbsnp = dbsnp.getRecord();}
  if (haveVcf) {successAnnVcf = annVcf.getRecord();}
  if (haveBed) {successBed = b.getRecord();}

// Finish when the end of the first file has been reached.
  while (successVcf) {

// If the end of the dbsnp vcf file is reached, write out the
// remaining records from the vcf file.
    if (!successDbsnp && !successAnnVcf && !successBed) {
      *output << v.record << endl;
      successVcf = v.getRecord();
      if (!successVcf) {break;}
    }

// If a dbSNP file was provided, parse until the position in the vcf file
// is found or passed.  Similarly for hapmap and bed files.
    if (haveDbsnp) {
      if (dbsnp.referenceSequence == currentReferenceSequence && v.position > dbsnp.position) {
        successDbsnp = dbsnp.parseVcf(v.referenceSequence, v.position, false, output);
      }
      if (dbsnp.referenceSequence == v.referenceSequence && v.position == dbsnp.position) {
        string tag = "VC";
        information sInfo = dbsnp.getInfo(tag);
        if (sInfo.values[0] == "SNP") {
          v.rsid = dbsnp.rsid;
          v.info += ";DBSNP";
        }
        build = true;
        successDbsnp = dbsnp.getRecord();
      }
    }
    if (haveVcf) {
      if (annVcf.referenceSequence == currentReferenceSequence && v.position > annVcf.position) {
        successAnnVcf = annVcf.parseVcf(v.referenceSequence, v.position, false, output);
      }
      if (v.referenceSequence == annVcf.referenceSequence && v.position == annVcf.position) {
        v.info += ";ANN=" + annVcf.info;;
        build = true;
        successAnnVcf = annVcf.getRecord();
      }
    }
    if (haveBed) {
      if (b.referenceSequence == currentReferenceSequence && v.position > b.end) {
        successBed = b.parseBed(v.referenceSequence, v.position);
      }
      if (v.referenceSequence == b.referenceSequence && v.position >= b.start && v.position <= b.end) {
        v.info += ";ANN=" + b.info;
        build = true;
      }
    }

    if (build) {
      newRecord = v.buildRecord(false);
      build = false;
    } else {
      newRecord = v.record;
    }
    *output << newRecord << endl;
    successVcf = v.getRecord();

// If the reference sequence in the vcf file changes, ensure that the annotation files are on the same
// reference sequence.
    if (v.referenceSequence != currentReferenceSequence) {
      currentReferenceSequence = v.referenceSequence;
      // dbSNP
      if (haveDbsnp) {
        if (dbsnp.referenceSequence != v.referenceSequence) {successDbsnp = dbsnp.parseVcf(v.referenceSequence, v.position, false, output);}
      }
      // Annotation vcf file
      if (haveVcf) {
        if (annVcf.referenceSequence != v.referenceSequence) {successAnnVcf = annVcf.parseVcf(v.referenceSequence, v.position, false, output);}
      }
      // bed file
      if (haveBed) {
        if (b.referenceSequence != v.referenceSequence) {successBed = b.parseBed(v.referenceSequence, v.position);}
      }
    }
  }
}

// Run the tool.
int annotateTool::Run(int argc, char* argv[]) {
  annotateDbsnp = false;
  annotateVcf = false;
  annotateBed = false;
  int getOptions = annotateTool::parseCommandLine(argc, argv);
  output = openOutputFile(outputFile);

  vcf v; // Define vcf object.
  v.openVcf(vcfFile);
  v.parseHeader();

  vcf dbsnp; // Define dbsnp vcf object.
  vcf annVcf; // Define dbsnp vcf object.
  bed b; // Define dbsnp vcf object.

  if (annotateDbsnp) {
    dbsnp.dbsnpVcf = true;
    dbsnp.processInfo = true;
    dbsnp.openVcf(dbsnpFile);
    dbsnp.parseHeader();
  }
  if (annotateVcf) {
    annVcf.openVcf(annVcfFile);
    annVcf.parseHeader();
  }
  if (annotateBed) {b.openBed(bedFile);}

// Add an extra line to the vcf header to indicate the file used for
// performing dbsnp annotation.
  string taskDescription = "##vcfCTools=annotated vcf file with ";
  if (annotateDbsnp) {taskDescription += "dbSNP file " + dbsnpFile;}
  if (annotateVcf || annotateBed) {

// Add the info line for the ANN tag in the header.  If an intersection with
// the input bed file is found and the fourth column of the bed file is CDS, 
// the string ANN=CDS will be added to the info string.  This needs to be
// included in the header to ensure the tool works correctly.
    if (annotateVcf) {
      v.headerInfoLine["ANN"] = "##INFO=<ID=ANN,Number=1,Type=String,Description=\"Annotation from vcf file " + annVcfFile + "\">";
    }
    if (annotateBed) {
      v.headerInfoLine["ANN"] = "##INFO=<ID=ANN,Number=1,Type=String,Description=\"Annotation from bed file " + bedFile + "\">";
    }

    if (annotateDbsnp) {taskDescription += ", ";}
    if (annotateVcf) {taskDescription += "vcf file " + annVcfFile;}
    if (annotateBed && annotateVcf) {taskDescription += ", bed file " + bedFile;}
    else if (annotateBed) {taskDescription += "bed file " + bedFile;}
  }
  writeHeader(output, v, false, taskDescription); // tools.cpp

// Annotate the vcf file.
  annotate(v, dbsnp, annVcf, b, annotateDbsnp, annotateVcf, annotateBed, output);

// Check that the input files had the same list of reference sequences.
// If not, it is possible that there were some problems.
  if (annotateDbsnp) {checkReferenceSequences(v.referenceSequenceVector, dbsnp.referenceSequenceVector);} // tools.cpp
  if (annotateVcf) {checkReferenceSequences(v.referenceSequenceVector, annVcf.referenceSequenceVector);} // tools.cpp
  if (annotateBed) {checkReferenceSequences(v.referenceSequenceVector, b.referenceSequenceVector);} // tools.cpp

// Close the vcf files.
  v.closeVcf();
  dbsnp.closeVcf();
  annVcf.closeVcf();
  b.closeBed();

  return 0;
}
