// ******************************************************
// vcfCTools (c) 2011 Alistair Ward
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ------------------------------------------------------
// Last modified: 1 March 2011
// ------------------------------------------------------
// Filter the vcf file on user specific criteria.  vcf
// records that fail the filters have the filter field
// populated with a semi-colon seperated list of failed
// filters.
// ******************************************************

#include "tool_filter.h"

using namespace std;
using namespace vcfCTools;

// intersectTool imlementation.
filterTool::filterTool(void)
  : AbstractTool()
{
  recordsInMemory = 100;
  markPass = false;
  filterFail = false;
  filterQuality = false;
  findHets = false;
  removeGenotypes = false;
  removeInfo = false;
  stripRecords = false;
  conditionalFilter = false;
  filterString = "";
}

// Destructor.
filterTool::~filterTool(void) {}

// Help
int filterTool::Help(void) {
  cout << "Template help" << endl;
  cout << "Usage: ./vcfCTools filter [options]." << endl;
  cout << endl;
  cout << "Options:" << endl;
  cout << "  -h, --help" << endl;
  cout << "	display intersect help." << endl;
  cout << "  -i, --in" << endl;
  cout << "	input vcf files (two, or one if intersecting with bed file)." << endl;
  cout << "  -o, --out" << endl;
  cout << "	output vcf file." << endl;
  cout << "  -d, --delete-info" << endl;
  cout << "	remove the entry from the info field." << endl;
  cout << "  -e, --find-hets" << endl;
  cout << "	output the samples that are hets to the info string." << endl;
  cout << "  -f, --filter-string" << endl;
  cout << "	a conditional statement on which to filter (enclosed in quotation marks)." << endl;
  cout << "  -l, --fail-filter" << endl;
  cout << "	filter out all variants that are not marked as 'PASS'." << endl;
  cout << "  -m, --mark-as-pass" << endl;
  cout << "	mark all records as 'PASS'." << endl;
  cout << "  -q, --quality" << endl;
  cout << "	filter on variant quality." << endl;
  cout << "  -r, --remove-genotypes" << endl;
  cout << "	do not include genotypes in the output vcf file." << endl;
  cout << "  -s, --strip-records" << endl;
  cout << "	strip out records containing the specified info field (comma separated list)." << endl;
  cout << endl;
  exit(0);

  return 0;
}

// Parse the command line and get all required and optional arguments.
int filterTool::parseCommandLine(int argc, char* argv[]) {
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
      {"delete-info", required_argument, 0, 'd'},
      {"filter-string", required_argument, 0, 'f'},
      {"fail-filter",no_argument, 0, 'l'},
      {"find-hets",no_argument, 0, 'e'},
      {"mark-as-pass", no_argument, 0, 'm'},
      {"quality", required_argument, 0, 'q'},
      {"remove-genotypes", required_argument, 0, 'r'},
      {"strip-records", required_argument, 0, 's'},

      {0, 0, 0, 0}
    };

    int option_index = 0;
    argument = getopt_long(argc, argv, "hi:o:d:f:elq:mrs:", long_options, &option_index);

    if (argument == -1) {break;}
    switch (argument) {

      // Input vcf file - required input.
      case 'i':
        vcfFile = optarg;
        break;

      // Output vcf file.
      case 'o':
        outputFile = optarg;
        break;

      // Remove an entry from the info field.
      case 'd':
        removeInfo = true;
        removeInfoString = optarg;
        break;

      // Read in the conditinal filter string.
      case 'f':
        conditionalFilter = true;
        filterString = string(optarg);
        break;

      case 'e':
        findHets = true;
        break;

      // Filter out reads that are not marked as 'PASS'
      case 'l':
        filterFail = true;
        break;

      // Filter on SNP quality.
      case 'q':
        filterQuality = true;
        char* value;
        value = optarg;
        filterQualityValue = atof(value);
        break;

      // Mark all records as passed filters.
      case 'm':
        markPass = true;
        break;

      // Remove genotypes from the output file.
      case 'r':
        removeGenotypes = true;
        break;
 
      // Strip out records with the specified info field.
      case 's':
        stripRecords = true;
        stripInfo = optarg;
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

  return 0;
}

// Perform the required filtering tasks on the variants.
void filterTool::filter(vcf& v) {
  information info;

  for (vector<variantDescription>::iterator iter = v.variantsIter->second.begin(); iter != v.variantsIter->second.end(); iter++) {
    string filterString = "";

// Split up the info string if the information is required.
    if (stripRecords || removeInfo) {
      v.processInfoFields(iter->info);
      info = v.getInfo(removeInfoString);
    }

// Mark the record as "PASS" if --mark-as-pass was applied.
    if (markPass) {filterString = "PASS";}

// Check for quality filtering.
    if (filterQuality && !markPass) {
      if (iter->quality < filterQualityValue) {
        ostringstream s;
        s << filterQualityValue;
        filterString = (filterString != "") ? filterString + ";" + "Q" + s.str() : "Q" + s.str();
      }
    }

// Remove info fields if necessary.
    if (removeInfo) {
      if (v.infoTags.count(removeInfoString) > 0) {
        info = v.getInfo(removeInfoString);
        string infoTag = removeInfoString;
        if (v.headerInfoFields[removeInfoString].type != "Flag") {infoTag += "=";}
        size_t found = iter->info.find(infoTag);
        size_t end   = iter->info.find(";", found + 1);
        string replaceString = "";

        // If the info tag is at the end of the line, remove the trailing ";".
        // Also, if this is the only entry in the info string, replace the
        // info string with ".".
        if (end != string::npos) {end = end - found + 1;}
        else {
        if (found != 0) {found--;}
          else {replaceString = ".";}
          end = iter->info.length();
        }
        iter->info.replace(found, end, replaceString);
      }
    }

// If genotypes need to be analysed, do that here.
    if (findHets) {
      string foundGeno = "";
      unsigned int count = 0;

      vector<string> geno = split(iter->genotypeString, "\t");
      for (vector<string>::iterator genIter = geno.begin(); genIter != geno.end(); genIter++) {
        foundGeno = (split(*genIter, ":") )[genotypePosition];
        if (foundGeno == "0/1" || foundGeno == "1/0") {
          iter->info += ";HET=" + v.samples[count];
        }
        count++;
      }
    }

    filterString = (filterString == "") ? iter->filters : filterString;
    iter->filters = filterString;

    writeRecord = true;
    if (filterFail && iter->filters != "PASS") {writeRecord = false;}
    if (stripRecords) {
      for (vector<string>::iterator stripIter = stripInfoList.begin(); stripIter != stripInfoList.end(); stripIter++) {
        if (v.infoTags.count(*stripIter) != 0) {
          writeRecord = false;
          break;
        }
      }
    }
    if (writeRecord) {
      string record = v.buildRecord(*iter);
      *output << record << endl;
    }
  }
}

// Run the tool.
int filterTool::Run(int argc, char* argv[]) {
  int getOptions = filterTool::parseCommandLine(argc, argv);

  vcf v; // Define vcf object.
  v.openVcf(vcfFile); // Open the vcf file.
  output = openOutputFile(outputFile);

// Read in the header information.
  v.parseHeader();
  string taskDescription = "##vcfCTools=filter";
  if (markPass) {taskDescription += "marked all records as PASS";}

// If find hets is specified, check that genotypes exist.
  if (findHets && !v.hasGenotypes) {
    cerr << "Input files does not contain genotype information." << endl;
    cerr << "Cannot identify sample genotypes." << endl;
    exit(1);
  }

// If the genotypes are to be removed, set the removeGenotypes value to 
// true for the vcf object.
  if (removeGenotypes) {v.removeGenotypes = true;}

// If records are to be stripped out of the vcf file, check the inputted
// IDs and populate the list of IDs to be stripped.
  if (stripRecords) {
    size_t found = stripInfo.find(",");
    stripInfoList.clear();
    if (found == string::npos) {stripInfoList.push_back(stripInfo);}
    else {stripInfoList = split(stripInfo, ",");}
    for (vector<string>::iterator iter = stripInfoList.begin(); iter != stripInfoList.end(); iter++) {
      if (v.headerInfoFields.count(*iter) == 0) {
        cerr << "WARNING: Info ID " << *iter << " is to be stripped, but does not appear in the header." << endl;
      }
    }
  }
  
// If a conditional filter string has been included, build up the methods to
// evaulate it.
  //if (conditionalFilter) {
  //}

// Write out the header.
  writeHeader(output, v, removeGenotypes, taskDescription);

// If information is being removed from the info string, ensure
// that the info string is processed.
  //if (removeInfo || stripRecords) {v.processInfo = true;}

// Read through all the entries in the file.  First construct the
// structure to contain the variants in memory and populate.
  v.success = v.getRecord();
  string currentReferenceSequence = v.variantRecord.referenceSequence;
  v.success = v.buildVariantStructure(recordsInMemory, currentReferenceSequence, false, output);

// If the genotypes need to be found, determine which position in the genotype
// format string the genotypes appear.
  if (findHets) {
    unsigned int count = 0;
    v.variantsIter = v.variants.begin();
    vector<variantDescription>::iterator i = v.variantsIter->second.begin();
    vector<string> s = split(i->genotypeFormatString, ":");
    for (vector<string>::iterator iter = s.begin(); iter != s.end(); iter++) {
      if (*iter == "GT") {break;}
      count++;
    }
    genotypePosition = count;
  }

  while (v.success) {
    if (v.variantRecord.referenceSequence == currentReferenceSequence) {

      // Add the new variant to the structure and process the first element,
      // finally erasing this element.  The structure should consequently
      // maintain the same size.
      v.addVariantToStructure();
      v.variantsIter = v.variants.begin();

// Perform all filtering tasks on this variant.
      filter(v);

      v.variants.erase(v.variantsIter);
      v.success = v.getRecord();
    } else {

    // Process all of the elements remaining in the variants structure (i.e.
    // those for the current reference sequence), deleting them as they are
    // processed, then set the currentReferenceSequence to the new value and
    // rebuild the variant structure.
      for (v.variantsIter = v.variants.begin(); v.variantsIter != v.variants.end(); v.variantsIter++) {

// Perform all filtering tasks on this variant.
        filter(v);
        v.variants.erase(v.variantsIter);
      }
      currentReferenceSequence = v.variantRecord.referenceSequence;
      v.success = v.buildVariantStructure(recordsInMemory, currentReferenceSequence, false, output);
    }
  }

// Close the vcf files.
  v.closeVcf();

  return 0;
}
