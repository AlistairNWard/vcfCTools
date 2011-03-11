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
  cout << "  -f, --priority-file" << endl;
  cout << "	output records from the vcf file (default: record with highest quality)." << endl;
  cout << endl;

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
    {"priority-file", required_argument, 0, 'f'},
    {"in", required_argument, 0, 'i'},
    {"out", required_argument, 0, 'o'},

    {0, 0, 0, 0}
  };

    int option_index = 0;
    argument = getopt_long(argc, argv, "hb:f:i:o:", long_options, &option_index);

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

      // Determine which file takes priority.
      case 'f':
        priorityFile = optarg;
        break;

      // Output file.
      case 'o':
        outputFile = optarg;
        break;

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

// Check that either two vcf files or one vcf and one bed file is specified.
  if (vcfFiles.size() == 0 || (vcfFiles.size() == 1 and bedFile == "") || (bedFile != "" && vcfFiles.size() != 1) || vcfFiles.size() > 2) {
    cerr << "Two vcf files or a vcf and a bed file must be specified (--in, -i, --bed, -b)." << endl;
    exit(1);
  }

  return 0;
}

void intersectTool::intersectVcf(vcf& v1, vcf& v2, unsigned int priority, ostream* output) {
  bool success1 = v1.getRecord();
  bool success2 = v2.getRecord();
  string currentReferenceSequence = v1.referenceSequence;

// s soon as the end of either file is reached, there can be no
// more intersecting SNPs, so terminate.
  while (success1 && success2) {
    cout << v1.referenceSequence << ":" << v1.position << " " << v2.referenceSequence << ":" << v2.position << endl;
    if (v1.referenceSequence == v2.referenceSequence && v1.referenceSequence == currentReferenceSequence) {
      cout << "EQUAL" << endl;
      if (v1.position == v2.position) {

// It is possible that there are multiple records for this locus in each file.
// Store the record and relevant info (e.g. variant type, length, alleles, quality)
// and move to the next record and check if it corresponds to the same locus.
// The intersection is only true for variants of the same class (i.e. if two records
// have a variant at ref:pos, but one is an indel and the other a SNP, the record
// will not be included in the intersection).  After gathering all variants for this
// locus, check the intersections.
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
            exit(0);
          }
          else {
            // Both files have a SNP at this locus.  If the reference alleles
            // are different, there is an error.  If the alternate alleles are
            // the same, use priority to determine which record to output.
            if (snpsAtLocus1[0].ref != snpsAtLocus2[0].ref) {
              cerr << "SNPs at " << currentReferenceSequence << ":" << currentPosition;
              cerr << " have different reference alleles (" << snpsAtLocus1[0].ref << "/";
              cerr << snpsAtLocus2[0].ref << endl;
              cerr << "Please check input vcf files." << endl;
              exit(0);
            }
            // Different alternates.
            if (snpsAtLocus1[0].alt != snpsAtLocus2[0].alt) {
              cerr << "Different SNP alt alleles not yet handled " << currentReferenceSequence << ":" << currentPosition << endl;
              exit(0);
            }
            // Both SNPs have the same reference and alternate alleles.
            else {
              *output << snpsAtLocus1[0].record << endl;
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
            for (sVariant2 = indelsAtLocus2.begin(); sVariant2 != indelsAtLocus2.end(); sVariant2++) {
              if ((*sVariant1).ref == (*sVariant2).ref && (*sVariant1).alt == (*sVariant2).alt) {
                *output << (*sVariant1).record << endl;
                indelsAtLocus2.erase(sVariant2);
                break;
              }
            }
          }
          indelsAtLocus1.clear();
          indelsAtLocus2.clear();
        }

        //writeVcfRecord(priority, v1, v2, output); // tools.cpp
      }
      else if (v2.position > v1.position) {success1 = v1.parseVcf(v2.referenceSequence, v2.position, false, output);}
      else if (v1.position > v2.position) {success2 = v2.parseVcf(v1.referenceSequence, v1.position, false, output);}
    }
    else {
      if (v1.referenceSequence == currentReferenceSequence) {success1 = v1.parseVcf(v2.referenceSequence, v2.position, false, output);}
      else if (v2.referenceSequence == currentReferenceSequence) {success2 = v2.parseVcf(v1.referenceSequence, v1.position, false, output);}

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
    priority = setVcfPriority(priorityFile, vcfFiles); // tools.cpp
    vcf v1; // Create a vcf object.
    vcf v2; // Create a vcf object.
    vcf v3; // Generate a new vcf object that will contain the header information of the new file.
    
// Open the vcf files.
    v1.openVcf(vcfFiles[0]);
    v2.openVcf(vcfFiles[1]);

// Read in the header information.
    v1.parseHeader();
    v2.parseHeader();
    if (priority == 3) {
      mergeHeaders(v1, v2, v3); // tools.cpp
      v1.processInfo = true;
      v2.processInfo = true;
    }
    else {checkDataSets(v1, v2);}

// Check that the header for the two files contain the same samples.
    if (v1.samples != v2.samples) {
      cerr << "vcf files contain different samples (or sample order)." << endl;
      exit(1);
    }
    else {
      string taskDescription = "##vcfCTools=intersect " + vcfFiles[0] + ", " + vcfFiles[1];
      if (priority == 3) {writeHeader(output, v3, true, taskDescription);}
      else if ( (priority == 2 && v2.hasHeader) || !v1.hasHeader) {writeHeader(output, v2, true, taskDescription);} // tools.cpp
      else {writeHeader(output, v1, true, taskDescription);} // tools.cpp
    }

// Intersect the two vcf files.
    intersectVcf(v1, v2, priority, output);

// Check that the input files had the same list of reference sequences.
// If not, it is possible that there were some problems.
    checkReferenceSequences(v1.referenceSequenceVector, v2.referenceSequenceVector); // tools.cpp

// Close the vcf files.
    v1.closeVcf();
    v2.closeVcf();
  }

  return 0;
}
