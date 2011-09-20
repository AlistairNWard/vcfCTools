// ******************************************************
// vcfCTools (c) 2011 Alistair Ward
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ------------------------------------------------------
// Last modified: 18 February 2011
// ------------------------------------------------------
// Collate information.
// ******************************************************

#include "tool_distributions.h"

using namespace std;
using namespace vcfCTools;

// statsTool imlementation.
distributionsTool::distributionsTool(void)
  : AbstractTool()
{
  usePrimary               = false;
  secondaryQuality         = false;
  useDistributions         = false;
  useDistQ                 = false;
  currentReferenceSequence = "";
  processComplex           = false;
  processIndels            = false;
  processMnps              = false;
  processSnps              = false;
}

// Destructor.
distributionsTool::~distributionsTool(void) {}

// Help
int distributionsTool::Help(void) {
  cout << "Stats help" << endl;
  cout << "Usage: ./vcfCTools distributions [options]." << endl;
  cout << endl;
  cout << "Options:" << endl;
  cout << "  -h, --help" << endl;
  cout << "     display intersect help." << endl;
  cout << "  -i, --in" << endl;
  cout << "     input vcf file." << endl;
  cout << "  -o, --output" << endl;
  cout << "     output vcf file." << endl;
  cout << "  -d, --distribution" << endl;
  cout << "     generate a distibution of the supplied info fields (or QUAL)." << endl;
  cout << "  -p, --primary-information" << endl;
  cout << "     Collate statistics on this info field." << endl;
  cout << "  -s, --secondary-information" << endl;
  cout << "     store information about these fields (comma separated list)." << endl;
  cout << "  -1, --snps" << endl;
  cout << "     analyse SNPs." << endl;
  cout << "  -2, --mnps" << endl;
  cout << "     analyse MNPs." << endl;
  cout << "  -3, --indels" << endl;
  cout << "     analyse indels." << endl;
  cout << "  -4, --complex" << endl;
  cout << "     analyse complex events." << endl;
  return 0;
}

// Parse the command line and get all required and optional arguments.
int distributionsTool::parseCommandLine(int argc, char* argv[]) {
  commandLine = argv[0];
  for (int i = 2; i < argc; i++) {
    commandLine += " ";
    commandLine += argv[i];
  }

  int argument; // Counter for getopt.
  // Define the long options.
  static struct option long_options[] = {
    {"help", no_argument, 0, 'h'},
    {"in", required_argument, 0, 'i'},
    {"out", required_argument, 0, 'o'},
    {"distribution", required_argument, 0, 'd'},
    {"primary-information", required_argument, 0, 'p'},
    {"secondary-information", required_argument, 0, 's'},
    {"snps", no_argument, 0, '1'},
    {"mnps", no_argument, 0, '2'},
    {"indels", no_argument, 0, '3'},
    {"complex", no_argument, 0, '4'},

    {0, 0, 0, 0}
  };

  while (true) {
    int option_index = 0;
    argument = getopt_long(argc, argv, "hi:o:d:p:s:1234", long_options, &option_index);

    if (argument == -1)
      break;

    switch (argument) {
      // Input vcf file - required input.
      case 'i':
        vcfFile = optarg;
        break;

      // Help.
      case 'h':
        return Help();

      // Output file.
      case 'o':
        outputFile = optarg;
        break;

      // Distributions to generate.
      case 'd':
        useDistributions = true;
        distString = optarg;
        break;

      // Collate information on this info field.
      case 'p':
        usePrimary = true;
        primaryInfo = optarg;
        break;

      // Store information about these information fields.
      case 's':
        secondaryInfoString = optarg;
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

      // Analyse complex events.
      case '4':
        processComplex = true;
        break;
      
      //
      case '?':
        cerr << "Unknown option: " << argv[optind - 2] << endl;
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

// Check that a primary and at least one secondary info fields are
// specified.
  if ( !usePrimary && !useDistributions) {
    cerr << "Options --primary-information (-p) and secondary-information (-s) or " << endl;
    cerr << "--distributions (-d) must be specified." << endl;
    exit(1);
  }

  return 0;
}

// Collate the information.
void distributionsTool::distributions(vcf& v, variant& var) {

  // SNPs
  if (var.processSnps) {
    for (var.variantIter = var.vmIter->second.biSnps.begin(); var.variantIter != var.vmIter->second.biSnps.end(); var.variantIter++) {
      performCollate(v, var.vmIter->first, *var.variantIter);
    }
    for (var.variantIter = var.vmIter->second.multiSnps.begin(); var.variantIter != var.vmIter->second.multiSnps.end(); var.variantIter++) {
      performCollate(v, var.vmIter->first, *var.variantIter);
    }
  }

  // MNPs
  if (var.processMnps) {
    for (var.variantIter = var.vmIter->second.mnps.begin(); var.variantIter != var.vmIter->second.mnps.end(); var.variantIter++) {
      performCollate(v, var.vmIter->first, *var.variantIter);
    }
  }

  // Indels
  if (var.processIndels) {
    for (var.variantIter = var.vmIter->second.indels.begin(); var.variantIter != var.vmIter->second.indels.end(); var.variantIter++) {
      performCollate(v, var.vmIter->first, *var.variantIter);
    }
  }
}

void distributionsTool::performCollate(vcf& v, int position, variantDescription& varIter) {
  variantInfo info;
  string pi;

  info.processInfoFields(varIter.info);
  if (usePrimary) {
    if (info.values.size() > 1) {
      cerr << "Cannot distributions information if multiple values exist." << endl;
      exit(1);
    } else {
      if (primaryInfo == "QUAL") {
        ostringstream s;
        s << varIter.quality;
        pi = s.str();
      } else {
        info.getInfo(primaryInfo, varIter.referenceSequence, position);
        pi = info.values[0];
      }
    }
    distributionsdInfo[pi].number++;

    // Now get the secondary information.
    for (vector<string>::iterator sIter = secondaryInfo.begin(); sIter != secondaryInfo.end(); sIter++) {
      if (info.infoTags.count(*sIter) > 0) {
        info.getInfo(*sIter, varIter.referenceSequence, position);
        distributionsdInfo[pi].secondary[*sIter] += atof(info.values[0].c_str());
      }
    }

    // Get the variant quality information if it was requested.
    if (secondaryQuality) {
      distributionsdInfo[pi].secondary["QUAL"] += varIter.quality;
    }
  }

// Now build up the requested distributions.
  if (useDistributions) {
    for (vector<string>::iterator iter = distFields.begin(); iter != distFields.end(); iter++) {
      if (info.infoTags.count(*iter) > 0) {
        info.getInfo(*iter, varIter.referenceSequence, position);
        dist[*iter][info.values[0]]++;
      }
    }
    if (useDistQ) {distQ[varIter.quality]++;}
  }
}

// Write out the information.
void distributionsTool::writePrimaryInfo() {
  map<string, distributionsStruct>::iterator iter;
  map<string, double>::iterator dIter;

  // Header.
  *output << primaryInfo << "	number";
  iter = distributionsdInfo.begin();
  for (dIter = iter->second.secondary.begin(); dIter != iter->second.secondary.end(); dIter++) {
    *output << "	" << dIter->first;
  }
  *output << endl;

  // Data.
  for (iter = distributionsdInfo.begin(); iter != distributionsdInfo.end(); iter++) {
    *output << iter->first << "	" << iter->second.number;
    for (dIter = iter->second.secondary.begin(); dIter != iter->second.secondary.end(); dIter++) {
      *output << "	" << dIter->second / iter->second.number;
    }
    *output << endl;
  }
}

// Write out the distributions.
void distributionsTool::writeDistributions() {
  map<string, map<string, unsigned int> >::iterator dIter;
  map<string, unsigned int>::iterator ddIter;
  map<double, unsigned int>::iterator qIter;

  if (useDistQ) {
    *output << "Distribution of quality scores." << endl;
    *output << setw(8) << "QUAL";
    *output << setw(8) << "number";
    *output << endl;
    for (qIter = distQ.begin(); qIter != distQ.end(); qIter++) {
      *output << setw(8) << qIter->first;
      *output << setw(8) << qIter->second;
      *output << endl;
    }
    *output << endl;
  }

  if (dist.size() != 0) {*output << "Distributions of info field parameters." << endl;}
  for (dIter = dist.begin(); dIter != dist.end(); dIter++) {
    *output << setw(8) << dIter->first;
    *output << setw(8) << "number";
    *output << endl;
    for (ddIter = dIter->second.begin(); ddIter != dIter->second.end(); ddIter++) {
      *output << setw(8) << ddIter->first;
      *output << setw(8) << ddIter->second;
      *output << endl;
    }
    *output << endl;
  }
}

// Run the tool.
int distributionsTool::Run(int argc, char* argv[]) {
  int getOptions = distributionsTool::parseCommandLine(argc, argv);

  vcf v; // Create a vcf object.
  variant var; // Create a variant object;
  var.determineVariantsToProcess(processSnps, processMnps, processIndels, processComplex, false, true, false);

  v.openVcf(vcfFile);
  output = openOutputFile(outputFile);
  v.parseHeader();

// Check that the primary and all secondary information appear in the header.
// First break up the comma separate list and populate the vector 
// secondaryInfo. "Q" is allowed for variant quality.
  if (!usePrimary && secondaryInfoString != "") {
    cerr << "Cannot provide secondary information fields without a primary field." << endl;
    exit(1);
  }

  if (usePrimary) {
    secondaryInfo = split(secondaryInfoString, ",");
    if (v.headerInfoFields.count(primaryInfo) == 0) {cerr << "WARNING: No header information for " << primaryInfo << endl;}
    vector<string>::iterator iter = secondaryInfo.begin();
    while (iter != secondaryInfo.end()) {
      if (*iter == "QUAL") {
        secondaryQuality = true;
        iter = secondaryInfo.erase(iter);
      } else {
        if (v.headerInfoFields.count(*iter) == 0) {cerr << "WARNING: No header information for " << *iter << endl;}
        iter++;
      }
    }
  }

// If distributions are to be generated, break up the comma separated list of
// fields to generate distributions for and add to the vector distFields.  If
// the quality distribution is required, set distQ to true.
  if (useDistributions) {
    distFields = split(distString, ",");
    vector<string>::iterator iter = distFields.begin();
    while (iter != distFields.end()) {
      if (*iter == "QUAL") {
        useDistQ = true;
        iter = distFields.erase(iter);
      } else {
        if (v.headerInfoFields.count(*iter) == 0) {cerr << "WARNING: No header information for " << *iter << endl;}
        iter++;
      }
    }
  }

// Read through all the entries in the file.  First construct the
// structure to contain the variants in memory and populate.
  v.update = true;
  while (v.success) {
    // Build the variant structure for this reference sequence.
    if (var.variantMap.size() == 0) {
      currentReferenceSequence = v.variantRecord.referenceSequence;
      v.success = var.buildVariantStructure(v);
    }

    // Loop over the variant structure until it is empty.  While v.update is true,
    // i.e. when the reference sequence is still the current reference sequence,
    // keep adding variants to the structre.
    while (var.variantMap.size() != 0) {
      if (v.variantRecord.referenceSequence == currentReferenceSequence && v.success) {
        var.addVariantToStructure(v.position, v.variantRecord, false);
        v.success = v.getRecord();
      }
      var.vmIter = var.variantMap.begin();
      distributions(v, var);
      var.variantMap.erase(var.vmIter);
    }
  }

// Write out the distributions information.
  if (usePrimary) {writePrimaryInfo();}
  if (useDistributions) {writeDistributions();}

// Close the vcf file and return.
  v.closeVcf();

  return 0;
}
