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

#include <cmath>
#include "tool_temp.h"

using namespace std;
using namespace vcfCTools;

// intersectTool imlementation.
tempTool::tempTool(void)
  : AbstractTool()
{}

// Destructor.
tempTool::~tempTool(void) {}

// Help
int tempTool::Help(void) {
  cout << "Template help" << endl;
  cout << "Usage: ./vcfCTools temp [options]." << endl;
  cout << endl;
  cout << "Options:" << endl;
  cout << "  -h, --help" << endl;
  cout << "	display intersect help." << endl;
  cout << "  -i, --in" << endl;
  cout << "	input vcf files (two, or one if intersecting with bed file)." << endl;
  cout << "  -o, --out" << endl;
  cout << "	output vcf file." << endl;
  cout << "  -q, --quality" << endl;
  cout << "	filter on variant quality." << endl;
  cout << "  -r, --remove-genotypes" << endl;
  cout << "	do not include genotypes in the output vcf file." << endl;
  cout << endl;

  return 0;
}

// Parse the command line and get all required and optional arguments.
int tempTool::parseCommandLine(int argc, char* argv[]) {
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
      {"mark-as-pass", no_argument, 0, 'm'},
      {"quality", required_argument, 0, 'q'},
      {"remove-genotypes", required_argument, 0, 'r'},

      {0, 0, 0, 0}
    };

    int option_index = 0;
    argument = getopt_long(argc, argv, "hi:o:q:mr", long_options, &option_index);

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

      // Remove genotypes from the output file.
      case 'r':
        removeGenotypes = true;
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
int tempTool::Run(int argc, char* argv[]) {
  markPass = false;
  filterQuality = false;
  removeGenotypes = false;
  int getOptions = tempTool::parseCommandLine(argc, argv);

  output = openOutputFile(outputFile);
  
  vcf v; // Define vcf object.

// Open the vcf file.
  v.openVcf(vcfFile);

// Read in the header information.
  v.parseHeader();
  string taskDescription = "##vcfCTools=filter";
  if (markPass) {taskDescription += "marked all records as PASS";}

  writeHeader(output, v, removeGenotypes, taskDescription);
  v.processInfo = true;

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

    filterString = (filterString == "") ? v.filters : filterString;
    v.filters = filterString;

    //string tag = "SA";
    //information sInfo = v.getInfo(tag);
    //unsigned int ABA = atoi(sInfo.values[0].c_str());
    //tag = "ABR";
    //sInfo = v.getInfo(tag);
    //unsigned int ABR = atoi(sInfo.values[0].c_str());

    //double AB = double(ABR) / ( double(ABR) + double(ABA) );
    //long double phred;
    //if (ABR == 0 && ABA == 0) {
    //  AB = 1.0;
    //  phred = 1000;
    //}
    //else {
    //  double successes = double(ABA);
    //  double trials = double(ABA) + double(ABR);
    //  double prob = 0.5;
    //  long double ABP = log(0.5) + (-2 * pow(trials * prob - successes, 2) / trials);
    //  phred = -10 * M_LOG10E * ABP;
    //}

    //ostringstream build;
    size_t found = v.info.find(";SR=");
    if (found != string::npos) {
      size_t end = v.info.find_first_of(";",found + 1);
    //build << AB;
      string newValue = "";
      v.info.replace(found, end - found - 1, newValue);
    //build.str("");
    }

    found = v.info.find(";SA=");
    if (found != string::npos) {
      size_t end = v.info.find_first_of(";",found + 1);
    //build << phred;
      string newValue = "";
      v.info.replace(found, end - found - 1, newValue);
    //build.str("");
    }

    string record = v.buildRecord(removeGenotypes);
    *output << record << endl;
  }

// Close the vcf files.
  v.closeVcf();

  return 0;
}
