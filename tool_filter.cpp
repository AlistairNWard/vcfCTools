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
  markPass = false;
  filterFail = false;
  filterQuality = false;
  removeGenotypes = false;
  removeInfo = false;
  stripRecords = false;
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
  cout << "  -f, --fail-filter" << endl;
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
      {"fail-filter",no_argument, 0, 'f'},
      {"mark-as-pass", no_argument, 0, 'm'},
      {"quality", required_argument, 0, 'q'},
      {"remove-genotypes", required_argument, 0, 'r'},
      {"strip-records", required_argument, 0, 's'},

      {0, 0, 0, 0}
    };

    int option_index = 0;
    argument = getopt_long(argc, argv, "hi:o:d:fq:mrs:", long_options, &option_index);

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

      // Filter out reads that are not marked as 'PASS'
      case 'f':
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

// Run the tool.
int filterTool::Run(int argc, char* argv[]) {
  int getOptions = filterTool::parseCommandLine(argc, argv);

  output = openOutputFile(outputFile);

  vcf v; // Define vcf object.

// Open the vcf file.
  v.openVcf(vcfFile);

// Read in the header information.
  v.parseHeader();
  string taskDescription = "##vcfCTools=filter";
  if (markPass) {taskDescription += "marked all records as PASS";}

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
  
// Write out the header.
  writeHeader(output, v, removeGenotypes, taskDescription);

// If information is being removed from the info string, ensure
// that the info string is processed.
  if (removeInfo || stripRecords) {v.processInfo = true;}

//Parse the vcf file and check if any of the filters are failed.  If
// so, build up a string of failed filters.
  while (v.getRecord()) {
    string filterString = "";

// Mark the record as "PASS" if --mark-as-pass was applied.
    if (markPass) {filterString = "PASS";}

// Check for quality filtering.
    if (filterQuality && !markPass) {
      if (v.quality < filterQualityValue) {
        ostringstream s;
        s << filterQualityValue;
        filterString = (filterString != "") ? filterString + ";" + "Q" + s.str() : "Q" + s.str();
      }
    }

// Remove info fields if necessary.
    
    if (removeInfo && v.infoTags.count(removeInfoString) != 0) {
      string infoTag = removeInfoString;
      if (v.headerInfoFields[removeInfoString].type != "Flag") {infoTag += "=";}
      size_t found = v.info.find(infoTag);
      size_t end = v.info.find(";", found + 1);
      string replaceString = "";

      // If the info tag is at the end of the line, remove the trailing ";".
      // Also, if this is the only entry in the info string, replace the
      // info string with ".".
      if (end != string::npos) {end = end - found + 1;}
      else {
        if (found != 0) {found--;}
        else {replaceString = ".";}
        end = v.info.length();
      }
      v.info.replace(found, end, replaceString);
    }

    filterString = (filterString == "") ? v.filters : filterString;
    v.filters = filterString;

    writeRecord = true;
    if (filterFail && v.filters != "PASS") {writeRecord = false;}
    if (stripRecords) {
      for (vector<string>::iterator iter = stripInfoList.begin(); iter != stripInfoList.end(); iter++) {
        if (v.infoTags.count(*iter) != 0) {
          writeRecord = false;
          break;
        }
      }
    }
    if (writeRecord) {
      string record = v.buildRecord(removeGenotypes);
      *output << record << endl;
    }
  }

// Close the vcf files.
  v.closeVcf();

  return 0;
}
