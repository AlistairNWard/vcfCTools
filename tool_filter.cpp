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
{}

// Destructor.
filterTool::~filterTool(void) {}

// Help
int filterTool::Help(void) {
  cout << "Template help" << endl;
  cout << "Usage: ./vcfCTools template [options]." << endl;
  cout << endl;
  cout << "Options:" << endl;
  cout << "  -h, --help" << endl;
  cout << "	display intersect help." << endl;
  cout << "  -i, --in" << endl;
  cout << "	input vcf files (two, or one if intersecting with bed file)." << endl;
  cout << endl;

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
      {"quality", required_argument, 0, 'q'},
      {"mark-as-pass", no_argument, 0, 'm'},

      {0, 0, 0, 0}
    };

    int option_index = 0;
    argument = getopt_long(argc, argv, "hi:o:q:m", long_options, &option_index);

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

  return 0;
}

// Run the tool.
int filterTool::Run(int argc, char* argv[]) {
  markPass = false;
  filterQuality = false;
  int getOptions = filterTool::parseCommandLine(argc, argv);

  output = openOutputFile(outputFile);
  
  vcf v; // Define vcf object.

// Open the vcf file.
  v.openVcf(vcfFile);

// Read in the header information.
  v.parseHeader();
  string taskDescription = "##vcfCTools=";
  if (markPass) {taskDescription += "marked all records as PASS";}

  writeHeader(output, v, true, taskDescription);

//Parse the vcf file and check if any of the filters are failed.  If
// so, build up a string of failed filters.
  while (v.getRecord()) {
    string filterString = "";

// Mark the record as "PASS" if --mark-as-pass was applied.
    if (markPass) {v.filters = "PASS";}

// Check for quality filtering.
    if (filterQuality) {
      if (v.quality < filterQualityValue) {
        ostringstream s;
        s << filterQualityValue;
        filterString = (filterString != "") ? filterString + ";" + "Q" + s.str() : "Q" + s.str();
      }
    }

    filterString = (filterString == "") ? "PASS" : filterString;
    v.filters = filterString;
    string record = v.buildRecord(true);
    *output << record << endl;
  }

// Close the vcf files.
  v.closeVcf();

  return 0;
}
