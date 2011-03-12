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
{}

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
  while (true) {
    static struct option long_options[] = {
      {"help", no_argument, 0, 'h'},
      {"bed", required_argument, 0, 'b'},
      {"in", required_argument, 0, 'i'},
      {"out", required_argument, 0, 'o'},
      {"common", required_argument, 0, 'c'},
      {"union", required_argument, 0, 'u'},
      {"unique", required_argument, 0, 'q'},

      {0, 0, 0, 0}
    };

    int option_index = 0;
    argument = getopt_long(argc, argv, "hb:i:o:c:u:q:", long_options, &option_index);

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
    else if ( (findCommon || findUnion) && writeFrom == "q") {cerr << "Writing out records with the highest quality value." << endl;}
    else {
      cerr << "The file from which the records are to be written needs to be selected." << endl;
      cerr << "    a - write records from the first inputted file." << endl;
      cerr << "    b - write records from the second inputted file." << endl;
      if (findCommon) {cerr << "    q - write records with the highest variant quality." << endl;}
      exit(1);
    }
  }
    
  return 0;
}

void intersectTool::intersectVcf(vcf& v1, vcf& v2, ostream* output, bool findCommon, bool findUnion, bool findUnique, string writeFrom) {
  bool success1 = v1.getRecord();
  bool success2 = v2.getRecord();
  string currentReferenceSequence = v1.referenceSequence;

// Define writeWhileParse.  If parsing through one of the vcf files to find the next
// common coordinate, this Boolean will instruct the parseVcf routine on whether or
// not to write the records to the output file.
  bool writeWhileParse;

// As soon as the end of either file is reached, there can be no
// more intersecting SNPs, so terminate.
  while (success1 && success2) {
    if (v1.referenceSequence == v2.referenceSequence && v1.referenceSequence == currentReferenceSequence) {
      if (v1.position == v2.position) {

// It is possible that there are multiple records for this locus in each file.
// Store the record and relevant info (e.g. variant type, length, alleles, quality)
// and move to the next record and check if it corresponds to the same locus.
//
// For intersections, only output records that occur in both files and have the same
// ref and alt alleles.
//
// For unions, include all variants, but ensure that each variant is only output once.
//
// For unique fractions, only output those variants unique to one of the vcf files.
// A deletion in one file that has a deletion in the other file will be included if the
// ref and alt alleles are different.
        currentPosition = v1.position;

// Check the first vcf file.
        while (v1.referenceSequence == currentReferenceSequence && v1.position == currentPosition && success1) {
          s = setStoredVariant(v1);
          if (v1.isSNP) {snpsAtLocus1.push_back(s);}
          if (v1.isMNP) {mnpsAtLocus1.push_back(s);}
          if (v1.isDeletion || v1.isInsertion) {indelsAtLocus1.push_back(s);}
          success1 = v1.getRecord();
        }

// and the second vcf file.
        while (v2.referenceSequence == currentReferenceSequence && v2.position == currentPosition && success2) {
          s = setStoredVariant(v2);
          if (v2.isSNP) {snpsAtLocus2.push_back(s);}
          if (v2.isMNP) {mnpsAtLocus2.push_back(s);}
          if (v2.isDeletion || v2.isInsertion) {indelsAtLocus2.push_back(s);}
          success2 = v2.getRecord();
        }

// Now compare the contents of the two sets of variants at this locus.  Only compare 
// variants of the same class.

        vector<storedVariants>::iterator sVariant1;
        vector<storedVariants>::iterator sVariant2;
        //SNPs
        if (snpsAtLocus1.size() != 0 && snpsAtLocus2.size() != 0) {
          if (snpsAtLocus1.size() > 1 || snpsAtLocus2.size() > 1) {
            cerr << "Cannot handle multiple SNP alleles in separate records yet.  Coordinate: " << currentReferenceSequence;
            cerr << ":" << currentPosition << endl;
            exit(1);
          }
          else {
            // Both files have a SNP at this locus.  If the reference alleles
            // are different, there is an error.
            if (snpsAtLocus1[0].ref != snpsAtLocus2[0].ref) {
              cerr << "SNPs at " << currentReferenceSequence << ":" << currentPosition;
              cerr << " have different reference alleles (" << snpsAtLocus1[0].ref << "/";
              cerr << snpsAtLocus2[0].ref << endl;
              cerr << "Please check input vcf files." << endl;
              exit(1);
            }
            // Different alternates.
            if (snpsAtLocus1[0].alt != snpsAtLocus2[0].alt) {
              cerr << "Different SNP alt alleles not yet handled " << currentReferenceSequence << ":" << currentPosition << endl;
              exit(1);
            }
            // Both SNPs have the same reference and alternate alleles.
            else {
              if (!findUnique) {
                if (writeFrom == "a") {*output << snpsAtLocus1[0].record << endl;}
                else if (writeFrom == "b") {*output << snpsAtLocus2[0].record << endl;}
                else {
                  if (snpsAtLocus1[0].quality >= snpsAtLocus2[0].quality) {*output << snpsAtLocus1[0].record << endl;}
                  else {*output << snpsAtLocus2[0].record << endl;}
                }
              }
            }
          }
          snpsAtLocus1.clear();
          snpsAtLocus2.clear();
        }
        //MNPs
        if (mnpsAtLocus1.size() != 0 && mnpsAtLocus2.size() != 0) {
          cerr << "Not yet handling MNPs" << endl;
          mnpsAtLocus1.clear();
          mnpsAtLocus2.clear();
        }
        //indels
        if (indelsAtLocus1.size() != 0 && indelsAtLocus2.size() != 0) {
          for (sVariant1 = indelsAtLocus1.begin(); sVariant1 != indelsAtLocus1.end(); sVariant1++) {
            bool unique = true;
            for (sVariant2 = indelsAtLocus2.begin(); sVariant2 != indelsAtLocus2.end(); sVariant2++) {
              if ((*sVariant1).ref == (*sVariant2).ref && (*sVariant1).alt == (*sVariant2).alt) {
                unique = false;
                if (!findUnique) {
                  if (writeFrom == "a") {*output << sVariant1->record << endl;}
                  else if (writeFrom == "b") {*output << sVariant2->record << endl;}
                  else if (writeFrom == "q") {
                    if (sVariant1->quality >= sVariant2->quality) {*output << sVariant1->record << endl;}
                    else {*output << sVariant2->record << endl;}
                  }
                }
                indelsAtLocus2.erase(sVariant2);
                break;
              }
            }
            // If the variant from the first vcf file was not present in the second file, the Boolean
            // unique will be true.  If looking for the union of the files or the variants unique to
            // the first file, this variant should be written out.
            if (unique && ( (findUnique && writeFrom == "a") || findUnion) ) {*output << sVariant1->record << endl;}
          }

          // All elements of indelsAtLocus2 (the indels present in the second file), were deleted if the
          // same variant was present in the first file.  This means that all variants remaining in this
          // vector are unique to the second vcf file and should be written out if the union or variants
          // unique to the second file were requested.
          if (findUnion || (findUnique && writeFrom == "b")) {
            for (sVariant2 = indelsAtLocus2.begin(); sVariant2 != indelsAtLocus2.end(); sVariant2++) {
              *output << sVariant2->record << endl;
            }
          }
          indelsAtLocus1.clear();
          indelsAtLocus2.clear();
        }
      } else if (v2.position > v1.position) {
        writeWhileParse = false;
        if (findUnion || (findUnique && (writeFrom == "a" || writeFrom == "q") ) ) {writeWhileParse = true;}
        success1 = v1.parseVcf(v2.referenceSequence, v2.position, writeWhileParse, output);
      } else if (v1.position > v2.position) {
        writeWhileParse = false;
        if (findUnion || (findUnique && (writeFrom == "b" || writeFrom == "q") ) ) {writeWhileParse = true;}
        success2 = v2.parseVcf(v1.referenceSequence, v1.position, writeWhileParse, output);
      }
    }
    else {
      if (v1.referenceSequence == currentReferenceSequence) {
        writeWhileParse = false;
        if (findUnion || (findUnique && (writeFrom == "a" || writeFrom == "q") ) ) {writeWhileParse = true;}
        success1 = v1.parseVcf(v2.referenceSequence, v2.position, writeWhileParse, output);
      } else if (v2.referenceSequence == currentReferenceSequence) {
        writeWhileParse = false;
        if (findUnion || (findUnique && (writeFrom == "b" || writeFrom == "q") ) ) {writeWhileParse = true;}
        success2 = v2.parseVcf(v1.referenceSequence, v1.position, writeWhileParse, output);
      }

// If the last record for a reference sequence is the same for both vcf
// files, they will both have referenceSequences different from the
// current reference sequence.  Change the reference sequence to reflect
// this and proceed.
      else {
        if (v1.referenceSequence != v2.referenceSequence) {
          cerr << "ERROR: Reference sequences for both files are unexpectedly different." << endl;
          cerr << "Check that both files contain records for the following reference sequences:" << endl;
          cerr << "\t" << v1.referenceSequence << " and " << v2.referenceSequence << endl;
          exit(1);
        }
      }
      currentReferenceSequence = v1.referenceSequence;
    }
  }
}

// Intersect a vcf file and a bed file.  It is assumed that the 
// two files are sorted by genomic coordinates and the reference
// sequences are in the same order.
void intersectTool::intersectVcfBed(vcf& v, bed& b, ostream* output) {
  bool successb = b.getRecord();
  bool successv = v.getRecord();
  string currentReferenceSequence = v.referenceSequence;

// As soon as the end of the first file is reached, there are no
// more intersections and the program can terminate.
  while (successv && successb) {
    if (v.referenceSequence == b.referenceSequence) {
      if (v.position < b.start) {successv = v.parseVcf(b.referenceSequence, b.start, false, output);}
      else if (v.position > b.end) {successb = b.parseBed(v.referenceSequence, v.position);}
      else {
        *output << v.record << endl;
        successv = v.getRecord();
      }
    }
    else {
      if (v.referenceSequence == currentReferenceSequence) {successv = v.parseVcf(b.referenceSequence, b.start, false, output);}
      if (b.referenceSequence == currentReferenceSequence) {successb = b.parseBed(v.referenceSequence, v.position);}
      currentReferenceSequence = v.referenceSequence;
    }
  }
}

// Run the tool.
int intersectTool::Run(int argc, char* argv[]) {
  findCommon = false;
  findUnion  = false;
  findUnique = false;
  int getOptions = intersectTool::parseCommandLine(argc, argv);
  output = openOutputFile(outputFile);

// If intersection is between a vcf file and a bed file, create a vcf and a bed object
// and intersect.
  if (bedFile != "") {
    vcf v; // Create a vcf object.
    bed b; // Create a bed object.

    v.openVcf(vcfFiles[0]);
    b.openBed(bedFile);
    v.parseHeader();

// Write the header to the output file.
    string taskDescription = "##vcfCtools=intersect " + vcfFiles[0] + ", " + bedFile;
    writeHeader(output, v, true, taskDescription);

// Intersect the files.
    intersectVcfBed(v, b, output);

// Check that the input files had the same list of reference sequences.
// If not, it is possible that there were some problems.
    checkReferenceSequences(v.referenceSequenceVector, b.referenceSequenceVector); // tools.cpp

// Close the vcf file and return.
    v.closeVcf();
    b.closeBed();
  }
  else {
    vcf v1; // Create a vcf object.
    vcf v2; // Create a vcf object.
    
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
    }
    else {
      string taskDescription = "##vcfCTools=intersect " + vcfFiles[0] + ", " + vcfFiles[1];
      writeHeader(output, v1, true, taskDescription); // tools.cpp
    }

// Intersect the two vcf files.
    intersectVcf(v1, v2, output, findCommon, findUnion, findUnique, writeFrom);

// Check that the input files had the same list of reference sequences.
// If not, it is possible that there were some problems.
    checkReferenceSequences(v1.referenceSequenceVector, v2.referenceSequenceVector); // tools.cpp

// Close the vcf files.
    v1.closeVcf();
    v2.closeVcf();
  }

  return 0;
}
