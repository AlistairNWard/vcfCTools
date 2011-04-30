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
  recordsInMemory = 100;
  currentReferenceSequence = "";
  annotateDbsnp = false;
  annotateVcf = false;
  annotateBed = false;
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
  cout << "  -1, --snps" << endl;
  cout << "     analyse SNPs." << endl;
  cout << "  -2, --mnps" << endl;
  cout << "     analyse MNPs." << endl;
  cout << "  -3, --indels" << endl;
  cout << "     analyse indels." << endl;
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
    {"snps", no_argument, 0, '1'},
    {"mnps", no_argument, 0, '2'},
    {"indels", no_argument, 0, '3'},

    {0, 0, 0, 0}
  };

    int option_index = 0;
    argument = getopt_long(argc, argv, "hi:o:d:a:b:123", long_options, &option_index);

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
//  bool successVcf;
//  currentReferenceSequence = "";
//  successVcf = v.getRecord(currentReferenceSequence);
//  if (successVcf) {v.update = true;}
//  currentReferenceSequence = v.referenceSequence;
//  string tag;
//  information annInfo;
//
//// Get the next record from the requested annotation files.
//  bool successDbsnp = false;
//  bool successAnnVcf = false;
//  bool successBed = false;
//  bool build = false;
//  if (haveDbsnp) {successDbsnp = dbsnp.getRecord(currentReferenceSequence);}
//  if (haveVcf) {successAnnVcf = annVcf.getRecord(currentReferenceSequence);}
//  if (haveBed) {successBed = b.getRecord();}
//
//// Finish when the end of the first file has been reached.
//  while (successVcf) {
//    string alleles;
//    string dbsnpAlleles;
//
//// If the end of the annotation files is reached, write out the
//// remaining records from the vcf file.
//    if (!successDbsnp && !successAnnVcf && !successBed) {
//      *output << v.record << endl;
//      successVcf = v.getRecord(currentReferenceSequence);
//      if (!successVcf) {break;}
//    }
//
//// If a dbSNP file was provided, parse until the position in the vcf file
//// is found or passed.  Similarly for hapmap and bed files.
//    if (haveDbsnp) {
//      if (v.hasMultipleAlternates) {
//        cerr << "Not yet able to handle checking multiple alternate alleles." << endl;
//        cerr << "Reference sequence: " << v.referenceSequence << ", position: " << v.position << endl;
//        exit(1);
//      } else {
//
//        //Annotate biallelic SNPs.
//        if (v.isSNP[0]) {
//          if (dbsnp.referenceSequence == currentReferenceSequence && v.position > dbsnp.position) {
//            successDbsnp = dbsnp.parseVcf(v.referenceSequence, v.position, false, output, false);
//          }
//          if (dbsnp.referenceSequence == v.referenceSequence && v.position == dbsnp.position) {
//            tag = "VC";
//            annInfo = dbsnp.getInfo(tag);
//            if (annInfo.values[0] == "SNP") {
//              v.rsid = dbsnp.rsid;
//
//              // Check that dbSNP and the vcf have the same alleles.  Tag the info field
//              // with dbSNP if they match, dnSNPX otherwise.  If the site has multiple
//              // alternates, do not include in the statistics.
//              if (v.hasMultipleAlternates || dbsnp.hasMultipleAlternates) {
//                v.info += ";dbSNPM";
//              } else {
//                alleles = v.ref + v.alt[0];
//                dbsnpAlleles = dbsnp.ref + dbsnp.alt[0];
//                for (int i = 0; i < 2; i++) {
//                  alleles[i] = tolower(alleles[i]);
//                  dbsnpAlleles[i] = tolower(dbsnpAlleles[i]);
//                }
//                sort(alleles.begin(), alleles.end());
//                sort(dbsnpAlleles.begin(), dbsnpAlleles.end());
//                v.info = (dbsnpAlleles == alleles) ? v.info + ";dbSNP" : v.info + ";dbSNPX";
//              }
//            } else {
//              v.rsid = ".";
//            }
//            build = true;
//            successDbsnp = dbsnp.getRecord(currentReferenceSequence);
//          }
//        }
//      }
//    }
//
//    // If another vcf file is provided for annotations, parse this file.
//    if (haveVcf) {
//      if (annVcf.referenceSequence == currentReferenceSequence && v.position > annVcf.position) {
//        successAnnVcf = annVcf.parseVcf(v.referenceSequence, v.position, false, output, false);
//      }
//      if (v.referenceSequence == annVcf.referenceSequence && v.position == annVcf.position) {
//        tag = "ANN";
//        if (v.infoTags.count(tag) != 0) {
//
//          // If an annotation exists, add the new annotation to the end of the comma separated list.
//          size_t found = v.info.find("ANN=");
//          size_t end = v.info.find_first_of(";", found + 1);
//
//          // Find current values in ANN.
//          annInfo = v.getInfo(tag);
//          string newField = "ANN=";
//          for (vector<string>::iterator iter = annInfo.values.begin(); iter != annInfo.values.end(); iter++) {
//            newField += (*iter) + ",";
//          }
//          newField += annVcf.filters;
//          v.info.replace(found, end - found - 1, newField);
//        } else {
//          v.info += ";ANN=" + annVcf.filters;
//        }
//        build = true;
//        successAnnVcf = annVcf.getRecord(currentReferenceSequence);
//      }
//    }
//
//    //Finally, if a bed file is provided, parse this file and compare with the current 
//    //vcf record.
//    if (haveBed) {
//      if (b.referenceSequence == currentReferenceSequence && v.position > b.end) {
//        successBed = b.parseBed(v.referenceSequence, v.position);
//      }
//      if (v.referenceSequence == b.referenceSequence && v.position >= b.start && v.position <= b.end) {
//        tag = "ANN";
//        if (v.infoTags.count(tag) != 0) {
//
//          // If an annotation exists, add the new annotation to the end of the comma separated list.
//          size_t found = v.info.find("ANN=");
//          size_t end = v.info.find_first_of(";", found + 1);
//
//          // Find current values in ANN.
//          annInfo = v.getInfo(tag);
//          string newField = "ANN=";
//          for (vector<string>::iterator iter = annInfo.values.begin(); iter != annInfo.values.end(); iter++) {
//            newField += (*iter) + ",";
//          }
//          newField += b.info;
//          v.info.replace(found, end - found - 1, newField);
//        } else {
//          v.info += ";ANN=" + b.info;
//        }
//        build = true;
//      }
//    }
//
//    if (build) {
//      //newRecord = v.buildRecord(false);
//      build = false;
//    } else {
//      newRecord = v.record;
//    }
//    *output << newRecord << endl;
//    successVcf = v.getRecord(currentReferenceSequence);
//
//// If the reference sequence in the vcf file changes, ensure that the annotation files are on the same
//// reference sequence.
//    if (v.referenceSequence != currentReferenceSequence) {
//      currentReferenceSequence = v.referenceSequence;
//      // dbSNP
//      if (haveDbsnp) {
//        if (dbsnp.referenceSequence != v.referenceSequence) {successDbsnp = dbsnp.parseVcf(v.referenceSequence, v.position, false, output, false);}
//      }
//      // Annotation vcf file
//      if (haveVcf) {
//        if (annVcf.referenceSequence != v.referenceSequence) {successAnnVcf = annVcf.parseVcf(v.referenceSequence, v.position, false, output, false);}
//      }
//      // bed file
//      if (haveBed) {
//        if (b.referenceSequence != v.referenceSequence) {successBed = b.parseBed(v.referenceSequence, v.position);}
//      }
//    }
//  }
}

// Run the tool.
int annotateTool::Run(int argc, char* argv[]) {
  int getOptions = annotateTool::parseCommandLine(argc, argv);
  output = openOutputFile(outputFile);

  vcf v; // Define vcf object.
  v.openVcf(vcfFile);
  v.parseHeader();

  vcf dbsnp; // Define dbsnp vcf object.
  vcf annVcf; // Define dbsnp vcf object.
  bed b; // Define dbsnp vcf object.
  v.processInfo = true;

  if (annotateDbsnp) {
    dbsnp.processInfo = true;
    dbsnp.dbsnpVcf = true;
    dbsnp.openVcf(dbsnpFile);
    dbsnp.parseHeader();
  }
  if (annotateVcf) {
    annVcf.openVcf(annVcfFile);
    annVcf.parseHeader();
  }
  if (annotateBed) {
    b.openBed(bedFile);
  }

// Add an extra line to the vcf header to indicate the file used for
// performing dbsnp annotation.
  string taskDescription = "##vcfCTools=annotated vcf file with ";
  if (annotateDbsnp) {
    v.headerInfoLine["dbSNP"] = "##INFO=<ID=dbSNP,Number=0,Type=Flag,Description=\"Membership in dbSNP file " + dbsnpFile;
    v.headerInfoLine["dbSNP"] += " with common alleles.\">";
    v.headerInfoLine["dbSNPX"] = "##INFO=<ID=dbSNPX,Number=0,Type=Flag,Description=\"Membership in dbSNP file " + dbsnpFile;
    v.headerInfoLine["dbSNPX"] += " with different alleles.\">";
    v.headerInfoLine["dbSNPM"] = "##INFO=<ID=dbSNPM,Number=0,Type=Flag,Description=\"Membership in dbSNP file " + dbsnpFile;
    v.headerInfoLine["dbSNPM"] += ". Either the vcf or dbSNP entry show a variant with multiple alternate alleles.\">";
  }
  if (annotateVcf || annotateBed) {

// Add the info line for the ANN tag in the header.  If an intersection with
// the input bed file is found and the fourth column of the bed file is CDS, 
// the string ANN=CDS will be added to the info string.  This needs to be
// included in the header to ensure the tool works correctly.
    if (annotateVcf) {
      v.headerInfoLine["ANN"] = "##INFO=<ID=ANN,Number=.,Type=String,Description=\"Annotation from vcf file " + annVcfFile + "\">";
    }
    if (annotateBed) {
      v.headerInfoLine["ANN"] = "##INFO=<ID=ANN,Number=.,Type=String,Description=\"Annotation from bed file " + bedFile + "\">";
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
