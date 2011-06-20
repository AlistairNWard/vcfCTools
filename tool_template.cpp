// ******************************************************
// vcfCTools (c) 2011 Alistair Ward
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ------------------------------------------------------
// Last modified: 18 February 2011
// ------------------------------------------------------
// Template for tool creation
// ******************************************************

#include "tool_template.h"

using namespace std;
using namespace vcfCTools;

// intersectTool imlementation.
templateTool::templateTool(void)
  : AbstractTool()
{}

// Destructor.
templateTool::~templateTool(void) {}

// Help
int templateTool::Help(void) {
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
int templateTool::parseCommandLine(int argc, char* argv[]) {
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

    {0, 0, 0, 0}
  };

    int option_index = 0;
    argument = getopt_long(argc, argv, "hi:", long_options, &option_index);

    if (argument == -1) {break;}
    switch (argument) {

      // Input vcf file - required input.
      case 'i':
        vcfFile = optarg;
        break;

      // Help.
      case 'h':
        return Help();

      //
      case '?':
        cerr << "Unknown option: " << argv[optind - 1] << endl;
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

  return 0;
}

// Run the tool.
int templateTool::Run(int argc, char* argv[]) {
  int getOptions = templateTool::parseCommandLine(argc, argv);

  return 0;
}
