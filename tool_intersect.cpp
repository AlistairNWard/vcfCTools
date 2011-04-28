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
  passFilters = false;
  findCommon = false;
  findUnion  = false;
  findUnique = false;
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
  cout << "  -p, --pass-filters" << endl;
  cout << "	Only variants that pass filters are considered." << endl;
  cout << "  -1, --snps" << endl;
  cout << "     analyse SNPs." << endl;
  cout << "  -2, --mnps" << endl;
  cout << "     analyse MNPs." << endl;
  cout << "  -3, --indels" << endl;
  cout << "     analyse indels." << endl;
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
    {"in", required_argument, 0, 'i'},
    {"out", required_argument, 0, 'o'},
    {"pass-filters", no_argument, 0, 'p'},
    {"common", required_argument, 0, 'c'},
    {"union", required_argument, 0, 'u'},
    {"unique", required_argument, 0, 'q'},
    {"snps", no_argument, 0, '1'},
    {"mnps", no_argument, 0, '2'},
    {"indels", no_argument, 0, '3'},

    {0, 0, 0, 0}
  };

  while (true) {
    int option_index = 0;
    argument = getopt_long(argc, argv, "hb:i:o:pc:u:q:123", long_options, &option_index);

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

// Intersect two vcf files.  Intersect by variant position only.
void intersectTool::intersectVcf(vcf& v1, variant& var1, vcf& v2, variant& var2) {
  v1.update = true;
  v2.update = true;

  // The following Boolean flags are used when adding variants to the map.  They
  // determine whether or not to write out records that are unique to one of the
  // files.
  bool write1 = ( (findUnique && writeFrom == "a") || findUnion) ? true : false;
  bool write2 = ( (findUnique && writeFrom == "b") || findUnion) ? true : false;

  // Build the variant structures for each 
  var1.buildVariantStructure(v1, v1.variantRecord.referenceSequence, write1, output);
  var2.buildVariantStructure(v2, v2.variantRecord.referenceSequence, write2, output);

  // Set the pointers to the start of each variant map.
  var1.vmIter = var1.variantMap.begin();
  var2.vmIter = var2.variantMap.begin();

  // Parse and compare the two variant structures until the end of one of the files
  // is reached and the variant structure for that file is empty.
  while ( !(var1.variantMap.size() == 0 && !v1.success) && !(var2.variantMap.size() == 0 && !v2.success) ) {

    // If the two variant structures are built with the same reference sequence, compare
    // the contents and parse through all varians for this reference sequence.
    if (var1.vmIter->second.referenceSequence == var2.vmIter->second.referenceSequence) {
      currentReferenceSequence = var1.vmIter->second.referenceSequence;

      while (var1.variantMap.size() != 0 && var2.variantMap.size() != 0) {

        // Variants at the same locus.
        if (var1.vmIter->first == var2.vmIter->first) {
          if (!findUnique) {var1.writeVariants(output);}

          // Clear the compared variants from the structure and add the next one from 
          // the file into the structure if it is from the same reference sequence.
          var1.variantMap.erase(var1.vmIter);
          if (v1.variantRecord.referenceSequence == currentReferenceSequence && v1.success) {
            var1.addVariantToStructure(v1.position, v1.variantRecord);
            v1.success = v1.getRecord(currentReferenceSequence);
          }
          if (var1.variantMap.size() != 0) {var1.vmIter = var1.variantMap.begin();}

          var2.variantMap.erase(var2.vmIter);
          if (v2.variantRecord.referenceSequence == currentReferenceSequence && v2.success) {
            var2.addVariantToStructure(v2.position, v2.variantRecord);
            v2.success = v2.getRecord(currentReferenceSequence);
          }
          if (var2.variantMap.size() != 0) {var2.vmIter = var2.variantMap.begin();}

        // Variant from the first vcf file is at a larger coordinate than that in the
        // second vcf file.  Parse through the second file until the position is greater
        // than or equal to that in the second file.
        } else if (var1.vmIter->first > var2.vmIter->first) {
          if (write2) {var2.writeVariants(output);}
          var2.variantMap.erase(var2.vmIter);
          if (v2.variantRecord.referenceSequence == currentReferenceSequence && v2.success) {
            var2.addVariantToStructure(v2.position, v2.variantRecord);
            v2.success = v2.getRecord(currentReferenceSequence);
          }
          if (var2.variantMap.size() != 0) {var2.vmIter = var2.variantMap.begin();}
          
        // Variant from the first vcf file is at a smaller coordinate than that in the
        // second vcf file.
        } else if (var1.vmIter->first < var2.vmIter->first) {
          if (write1) {var1.writeVariants(output);}
          var1.variantMap.erase(var1.vmIter);
          if (v1.variantRecord.referenceSequence == currentReferenceSequence && v1.success) {
            var1.addVariantToStructure(v1.position, v1.variantRecord);
            v1.success = v1.getRecord(currentReferenceSequence);
          }
          if (var1.variantMap.size() != 0) {var1.vmIter = var1.variantMap.begin();}

        }
      }

      // Having finished comparing, there may still be variants left from one of the two files.
      // Check that the two variant structures are empty and if not, finish processing the
      // remaining variants for this reference sequence.
      if (var1.variantMap.size() != 0) {var1.clearReferenceSequence(v1, currentReferenceSequence, write1, output);}
      if (var2.variantMap.size() != 0) {var2.clearReferenceSequence(v2, currentReferenceSequence, write2, output);}

      // Now both variant maps are exhausted, so rebuild the maps with the variants from the
      // next reference sequence in the file.
      var1.buildVariantStructure(v1, v1.variantRecord.referenceSequence, write1, output);
      var2.buildVariantStructure(v2, v2.variantRecord.referenceSequence, write2, output);
 
      if (var1.variantMap.size() != 0) {var1.vmIter = var1.variantMap.begin();}
      if (var2.variantMap.size() != 0) {var2.vmIter = var2.variantMap.begin();}

    // If the variant structures are from different reference sequences, parse through the
    // second vcf file until the next reference sequence is found.  If finding the union
    // of the variants unique to the second vcf file, write them out.
    } else {
      if (var2.variantMap.size() != 0) {var2.clearReferenceSequence(v2, var2.vmIter->second.referenceSequence, write2, output);}
      var2.buildVariantStructure(v2, v2.variantRecord.referenceSequence, write2, output);
      if (var2.variantMap.size() != 0) {var2.vmIter = var2.variantMap.begin();}
    }
  }
}

// Intersect a vcf file and a bed file.  It is assumed that the 
// two files are sorted by genomic coordinates and the reference
// sequences are in the same order.  Do not group together variants
// in common reference sequence.
void intersectTool::intersectVcfBed(vcf& v, bed& b) {
  bool successBed = b.getRecord();
  bool successVcf = v.getRecord(currentReferenceSequence);
  string currentReferenceSequence = v.referenceSequence;
  v.update = true;

// As soon as the end of the first file is reached, there are no
// more intersections and the program can terminate.
  while (successVcf && successBed) {
    if (v.referenceSequence == b.referenceSequence) {
      if (v.position < b.start) {successVcf = v.parseVcf(b.referenceSequence, b.start, false, output, passFilters);}
      else if (v.position > b.end) {successBed = b.parseBed(v.referenceSequence, v.position);}
      else {
        *output << v.record << endl;
        successVcf = v.getRecord(currentReferenceSequence);
      }
    } else {
      if (v.referenceSequence == currentReferenceSequence) {successVcf = v.parseVcf(b.referenceSequence, b.start, false, output, passFilters);}
      if (b.referenceSequence == currentReferenceSequence) {successBed = b.parseBed(v.referenceSequence, v.position);}
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
    cerr << "Not updated bed intersection tool." << endl;
    exit(1);
    //vcf v; // Create a vcf object.
    //bed b; // Create a bed object.

    //v.openVcf(vcfFiles[0]);
    //b.openBed(bedFile);
    //v.parseHeader();

// Write the header to the output file.
    //string taskDescription = "##vcfCtools=intersect " + vcfFiles[0] + ", " + bedFile;
    //writeHeader(output, v, false, taskDescription);

// Intersect the files.
    //if (groupVariants) {
      //intersectVariantGroupsBed(v, b, output);
      //intersectVariantGroupsBed(v, b);
    //} else {
      //intersectVcfBed(v, b, output);
      //intersectVcfBed(v, b);
    //}

// Check that the input files had the same list of reference sequences.
// If not, it is possible that there were some problems.
    //checkReferenceSequences(v.referenceSequenceVector, b.referenceSequenceVector); // tools.cpp

// Close the vcf file and return.
    //v.closeVcf();
    //b.closeBed();
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
    intersectVcf(v1, var1, v2, var2);

// Check that the input files had the same list of reference sequences.
// If not, it is possible that there were some problems.
    checkReferenceSequences(v1.referenceSequenceVector, v2.referenceSequenceVector); // tools.cpp

// Close the vcf files.
    v1.closeVcf();
    v2.closeVcf();
  }

  return 0;
}
