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
  appliedFilters           = false;
  cleardbSnp               = false;
  conditionalFilter        = false;
  currentReferenceSequence = "";
  filterFail               = false;
  filterQuality            = false;
  filterString             = "";
  findHets                 = false;
  keepRecords              = false;
  markPass                 = false;
  processComplex           = false;
  processIndels            = false;
  processMnps              = false;
  processRearrangements    = false;
  processSnps              = false;
  processSvs               = false;
  removeGenotypes          = false;
  removeInfo               = false;
  splitMnps                = false;
  stripRecords             = false;
  useSampleList            = false;
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
  //cout << "  -c, --clear-dbsnp" << endl;
  //cout << "	clear the rsid field and remove dbSNP info flag." << endl;
  //cout << "  -d, --delete-info" << endl;
  //cout << "	remove the entry from the info field (comma separated list)." << endl;
  //cout << "  -e, --find-hets" << endl;
  //cout << "	output the samples that are hets to the info string." << endl;
  //cout << "  -k, --keep-records" << endl;
  //cout << "	only keep records containing the specified info field (comma separated list)." << endl;
  cout << "  -l, --fail-filter" << endl;
  cout << "	filter out all variants that are not marked as 'PASS'." << endl;
  cout << "  -m, --mark-as-pass" << endl;
  cout << "	mark all records as 'PASS'." << endl;
  cout << "  -p, --split-mnps" << endl;
  cout << "	split MNPs into SNPs." << endl;
  cout << "  -q, --quality" << endl;
  cout << "	filter on variant quality." << endl;
  cout << "  -r, --remove-genotypes" << endl;
  cout << "	do not include genotypes in the output vcf file." << endl;
  //cout << "  -s, --samples" << endl;
  //cout << "	output variants that occur in the provided list of samples.." << endl;
  //cout << "  -t, --strip-records" << endl;
  //cout << "	strip out records containing the specified info field (comma separated list)." << endl;
  cout << "  -1, --snps" << endl;
  cout << "	analyse SNPs." << endl;
  cout << "  -2, --mnps" << endl;
  cout << "	analyse MNPs." << endl;
  cout << "  -3, --indels" << endl;
  cout << "	analyse indels." << endl;
  cout << "  -4, --complex" << endl;
  cout << "	analyse complex events." << endl;
  cout << "  -5, --structural-variants" << endl;
  cout << "     analyse structural variantion events." << endl;
  cout << "  -6, --rearrangements" << endl;
  cout << "     analyse complex rearrangement events." << endl;
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
      {"split-mnps", no_argument, 0, 'p'},
      {"clear-dbsnp", no_argument, 0, 'c'},
      {"delete-info", required_argument, 0, 'd'},
      {"keep-records", required_argument, 0, 'k'},
      {"fail-filter",no_argument, 0, 'l'},
      {"find-hets",no_argument, 0, 'e'},
      {"mark-as-pass", no_argument, 0, 'm'},
      {"quality", required_argument, 0, 'q'},
      {"remove-genotypes", required_argument, 0, 'r'},
      {"samples", required_argument, 0, 's'},
      {"strip-records", required_argument, 0, 't'},
      {"snps", no_argument, 0, '1'},
      {"mnps", no_argument, 0, '2'},
      {"indels", no_argument, 0, '3'},
      {"complex", no_argument, 0, '4'},
      {"structural-variants", no_argument, 0, '5'},
      {"rearrangements", no_argument, 0, '6'},

      {0, 0, 0, 0}
    };

    int option_index = 0;
    argument = getopt_long(argc, argv, "hi:o:cd:k:elpq:mrs:t:123456", long_options, &option_index);

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

      // Clear the rsid field and the dbSNP info tag.
      case 'c':
        cleardbSnp = true;
        break;

      // Remove an entry from the info field.
      case 'd':
        removeInfo = true;
        removeInfoString = optarg;
        break;

      case 'e':
        findHets = true;
        break;

      // Only keep records containing the following info fields.
      case 'k':
        keepRecords = true;
        keepInfoFields = optarg;
        break;

      // Filter out reads that are not marked as 'PASS'
      case 'l':
        filterFail = true;
        break;

      // Mark all records as passed filters.
      case 'm':
        markPass = true;
        break;

      // When processing MNPs, split them up into
      // multiple SNPs.
      case 'p':
        splitMnps = true;
        break;

      // Filter on SNP quality.
      case 'q':
        filterQuality  = true;
        appliedFilters = true;
        char* value;
        value = optarg;
        filterQualityValue = atof(value);
        break;

      // Remove genotypes from the output file.
      case 'r':
        removeGenotypes = true;
        break;
 
      // Read in a list of samples.  Only variants that occur
      // in these samples will be considered.  This option
      // requires genotypes to be present.
      case 's':
        useSampleList = true;
        samplesListFile = optarg;
        break;

      // Strip out records with the specified info field.
      case 't':
        stripRecords = true;
        stripInfo = optarg;
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

      // Analyse complex.
      case '4':
        processComplex = true;
        break;

      // Analyse structural variants.
      case '5':
        processSvs = true;
        break;

      // Analyse complex rearrangements.
      case '6':
        processRearrangements = true;
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

// Check that a vcf file was specified.
   if (vcfFile == "") {
     cerr << "A vcf file must be specified (--in, -i)." << endl;
     exit(1);
   }

// If keepRecords and stripRecords have been simultaneously specified, terminate
// with an error.  This situation could lead to ambiguous decisions.
  if (stripRecords && keepRecords) {
    cerr << "ERROR: Cannot specify --strip-records (-t) and --keep-records (-k) simultaneously." << endl;
    exit(1);
  }

  return 0;
}

// Check if info tags appear in the header and build a list
// of info flags.
vector<string> filterTool::checkInfoFields(vcfHeader& header, vcf& v, string& infoString) {
  vector<string> infoList;

  size_t found = infoString.find(",");
  infoList.clear();
  if (found == string::npos) {infoList.push_back(infoString);}
  else {infoList = split(infoString, ",");}
  for (vector<string>::iterator iter = infoList.begin(); iter != infoList.end(); iter++) {
    if (header.infoFields.count(*iter) == 0) {
      cerr << "WARNING: Info ID " << *iter << " is to be stripped, but does not appear in the header." << endl;
    }
  }

  return infoList;
}

// Perform the required filtering tasks on the variants.
void filterTool::filter(variant& var) {

  // Loop over all records at this locus.
  var.ovIter = var.ovmIter->second.begin();
  for (; var.ovIter != var.ovmIter->second.end(); var.ovIter++) {

    // Define the filterString.  If the filter is '.', set the filterString
    // to a blank field, otherwise, set it equal to the current filter
    // string in the variant record.  This filter can be modified by the
    // various filter operations and will be used as the final value of the
    // filter field.
    string filterString = (var.ovIter->filters == ".") ? "" : var.ovIter->filters;
  
    // If the record is to be marked as having passed, change the filters
    // to PASS.
    if (markPass) {filterString = "PASS";}
  
    // Check for quality filtering.
    if (filterQuality && !markPass) {
      if (var.ovIter->quality < filterQualityValue) {
        ostringstream s;
        s << filterQualityValue;
        filterString = (filterString != "") ? filterString + ";" + "Q" + s.str() : "Q" + s.str();
      }
    }
  
    // If the filterString is blank, the record didn't fail any of the filters
    // so set the filter field to PASS, otherwise set it to the filter string.
    filterString = (filterString == "") ? "." : filterString;
    var.ovIter->filters = (filterString == "." && appliedFilters) ? "PASS" : filterString;

    // If only records marked as PASS as to be kept, check the value and update the
    // filtered vector as necessary.
    if (filterFail && var.ovIter->filters != "PASS") {
      vector<bool>::iterator fIter = var.ovIter->filtered.begin();
      for (; fIter != var.ovIter->filtered.end(); fIter++) {*fIter = true;}
    }
  }


//  // SNPs.
//  if (var.processSnps) {
//
//    // Set the reference allele as that from the first SNP in the list.  Each subsequent
//    // SNP will have the reference allele checked against this to ensure consistency.
//    // Also, an array containing all alternates is created and identical SNPs are filtered
//    // out.
//    //
//    // This checking is currently only performed for biallelic SNPs.
//    if (var.vmIter->second.biSnps.size() > 0) {
//      string masterRef = var.vmIter->second.biSnps[0].ref;
//      map<string, int> altAlleles;
//      bool duplicateAlts = false;
//  
//      for (var.variantIter = var.vmIter->second.biSnps.begin(); var.variantIter != var.vmIter->second.biSnps.end(); var.variantIter++) {
//        if (var.variantIter->ref != masterRef) {
//          cerr << "Different SNPs at the same locus have different alt alleles" << endl;
//          cerr << "Please check the contents of the vcf file" << endl;
//          cerr << var.variantIter->referenceSequence << ":" << var.vmIter->first;
//          cerr << ".  Present reference alleles: " << masterRef << " and " << var.variantIter->ref << endl;
//          exit(1);
//        }
//  
//        // Keep count of the number of occurences of each alternate allele.
//        altAlleles[var.variantIter->altString]++;
//        if (altAlleles[var.variantIter->altString] > 1) {duplicateAlts = true;}
//      }
//  
//      // If any alternate allele was observed more than once, some of the SNPs need to be
//      // removed as they are duplicated.
//      if (duplicateAlts) {
//        variantInfo info;
//  
//        for (map<string, int>::iterator iter = altAlleles.begin(); iter != altAlleles.end(); iter++) {
//  
//          // Find the number of alts for each alternate allele.  If this is greater than 1,
//          // loop through the SNPs and throw out all but 1 of these SNPs.  As an initial
//          // check, see if any of the duplicated SNPs occurred as a result of decomposing
//          // MNPs.  These should be removed first.  In the absence of this flag, the first
//          // record will be retained, all others discarded.
//          if (iter->second > 1) {
//            var.variantIter = var.vmIter->second.biSnps.begin();
//            while (var.variantIter != var.vmIter->second.biSnps.end()) {
//              if (var.variantIter->altString == iter->first) {
//                info.processInfoFields(var.variantIter->info);
//                if (info.infoTags.count("FROM_MNP") != 0) {
//                  var.variantIter = var.vmIter->second.biSnps.erase(var.variantIter);
//                  altAlleles[iter->first]--;
//                  cerr << "WARNING: Removed duplicated SNP at " << var.vmIter->second.referenceSequence << ":" << var.vmIter->first << endl;
//                } else {
//                  ++var.variantIter;
//                }
//              } else {
//                ++var.variantIter;
//              }
//            }
//          }
//        }
//  
//        // Check if there are still duplicated SNPs.  If so, keep only the first instance
//        // of each.
//        for (map<string, int>::iterator iter = altAlleles.begin(); iter != altAlleles.end(); iter++) {
//          int count = 0;
//          if (iter->second > 1) {
//            var.variantIter = var.vmIter->second.biSnps.begin();
//            while (var.variantIter != var.vmIter->second.biSnps.end()) {
//              if (var.variantIter->altString == iter->first) {
//                if (count != 0) {
//                  var.variantIter = var.vmIter->second.biSnps.erase(var.variantIter);
//                  altAlleles[iter->first]--;
//                  cerr << "WARNING: Removed duplicated SNP at " << var.vmIter->second.referenceSequence << ":" << var.vmIter->first << endl;
//                } else {
//                  ++var.variantIter;
//                }
//                count++;
//              } else {
//                ++var.variantIter;
//              }
//            }
//          }
//        }
//      }
//    }
//
//    // Perform filtering on the SNPs.
//    for (var.variantIter = var.vmIter->second.biSnps.begin(); var.variantIter != var.vmIter->second.biSnps.end(); var.variantIter++) {
//      performFilter(v, var.vmIter->first, *var.variantIter);
//    }
//    for (var.variantIter = var.vmIter->second.multiSnps.begin(); var.variantIter != var.vmIter->second.multiSnps.end(); var.variantIter++) {
//      performFilter(v, var.vmIter->first, *var.variantIter);
//    }
//  }
//
//  // MNPs.
//  if (var.processMnps) {
//    for (var.variantIter = var.vmIter->second.mnps.begin(); var.variantIter != var.vmIter->second.mnps.end(); var.variantIter++) {
//      performFilter(v, var.vmIter->first, *var.variantIter);
//    }
//  }
//
//  // Indels.
//  if (var.processIndels) {
//    for (var.variantIter = var.vmIter->second.indels.begin(); var.variantIter != var.vmIter->second.indels.end(); var.variantIter++) {
//      performFilter(v, var.vmIter->first, *var.variantIter);
//    }
//  }
}

void filterTool::performFilter(vcf& v, int position, variantDescription& varIter) {
//  variantInfo info;
//  vector<string> geno;
//
//  // Set the filter string to blank if the variant is being filtered rather than having
//  // info removed or hets found etc.  If no actual filtering is taking place, the filter
//  // field should not be modified but left as is.
//  string filterString = (!filterQuality && !conditionalFilter) ? varIter.filters : "";
//
//  // If a samples list was provided, check that one of the listed samples shows evidence
//  // for the alternate allele.  If none of the samples do, do not output this record.
//  if (useSampleList || findHets) {geno = split(varIter.genotypeString, "\t");}
//  if (useSampleList) {
//    cout << "CHECK " << geno[1] << endl;
//    exit(0);
// }
//
//  // Split up the info string if the information is required.
//  if (stripRecords || removeInfo || keepRecords) {info.processInfoFields(varIter.info);}
//
//  // Mark the record as "PASS" if --mark-as-pass was applied.
//  if (markPass) {filterString = "PASS";}
//
//  // Check for quality filtering.
//  if (filterQuality && !markPass) {
//    if (varIter.quality < filterQualityValue) {
//      ostringstream s;
//      s << filterQualityValue;
//      filterString = (filterString != "") ? filterString + ";" + "Q" + s.str() : "Q" + s.str();
//    }
//  }
//
//  // Remove info fields if necessary.
//  if (removeInfo) {
//    for (vector<string>::iterator iter = removeInfoList.begin(); iter != removeInfoList.end(); iter++) {
//      if (info.infoTags.count(*iter) > 0) {
//        info.getInfo(*iter, varIter.referenceSequence, position);
//        size_t found = varIter.info.find(*iter);
//        size_t end   = varIter.info.find(";", found + 1);
//        string replaceString = "";
//
//        // If the info tag is at the end of the line, remove the trailing ";".
//        // Also, if this is the only entry in the info string, replace the
//        // info string with ".".
//        if (end != string::npos) {
//          end = end - found + 1;
//        } else {
//          if (found != 0) {found--;}
//          else {replaceString = ".";}
//          end = varIter.info.length();
//        }
//        varIter.info.replace(found, end, replaceString);
//      }
//    }
//  }
//
//  // If genotypes need to be analysed, do that here.
//  if (findHets) {
//    string foundGeno = "";
//    string hetString = "";
//    unsigned int count = 0;
//
//    for (vector<string>::iterator genIter = geno.begin(); genIter != geno.end(); genIter++) {
//      foundGeno = (split(*genIter, ":") )[genotypePosition];
//      if (foundGeno == "0/1" || foundGeno == "1/0") {
//        hetString = (hetString == "") ? v.samples[count] : hetString + "," + v.samples[count];
//      }
//      count++;
//    }
//    if (hetString != "" ) {varIter.info += ";HET=" + hetString;}
//  }
//
//  // If the filterString is blank, the record didn't fail any of the filters
//  // so set it the variant filter field to PASS, otherwise set it to
//  // the filter string.
//  varIter.filters = (filterString == "") ? "PASS" : filterString;
//
//  writeRecord = true;
//  if (filterFail && varIter.filters != "PASS") {writeRecord = false;}
//
//  // Check if this record is to be stripped out.
//  if (stripRecords) {
//    for (vector<string>::iterator stripIter = stripInfoList.begin(); stripIter != stripInfoList.end(); stripIter++) {
//      if (info.infoTags.count(*stripIter) != 0) {
//        writeRecord = false;
//        break;
//      }
//    }
//  }
//
//  // Check if this record should be kept.
//  if (keepRecords) {
//    writeRecord = false;
//    for (vector<string>::iterator keepIter = keepInfoList.begin(); keepIter != keepInfoList.end(); keepIter++) {
//      if (info.infoTags.count(*keepIter) != 0) {
//        writeRecord = true;
//        break;
//      }
//    }
//  }
//
//  // If the dbSNP informaion is to be removed, set the rsid value to '.',
//  if (cleardbSnp) {varIter.rsid = ".";}
//
//  // If genotypes are to be removed, set the genotype strings to be blank.
//  if (removeGenotypes) {
//    varIter.genotypeFormatString = "";
//    varIter.genotypeString       = "";
//  }
//
//  if (writeRecord) {
//    buildRecord(position, varIter);
//    *output << varIter.record << endl;
//  }
}

// Run the tool.
int filterTool::Run(int argc, char* argv[]) {
  int getOptions = filterTool::parseCommandLine(argc, argv);

  // Depending on the filtering being performed, it may or may not be necessary
  // to look at each individual allele.  For example, if the only action is to
  // filter out genotypes, the different alternate alleles do not need to be
  // interrogated.
  //bool processAlleles = (stripRecords || findHets || keepRecords || splitMnps || useSampleList) ? true : false;
  bool processAlleles = true;

  // Define an output object and open the output file.
  output ofile;
  ofile.outputStream = ofile.openOutputFile(outputFile);

  // Define the vcf object.
  vcf v;
  v.openVcf(vcfFile);

  // Define the header object and read in the header information.
  vcfHeader header;
  header.parseHeader(v.input);

  // Define the variant object.
  variant var;
  var.determineVariantsToProcess(processSnps, processMnps, processIndels, processComplex, processSvs, processRearrangements, splitMnps, processAlleles, false);

  string taskDescription = "##vcfCTools=filter";
  if (markPass) {taskDescription += "marked all records as PASS";}

  // If MNPs are to be broken up, add a line to the header explaining this.
  if (splitMnps) {
    string text = "Indicates that this SNP was generated from the decomposition of an MNP.\">";
    header.infoLine["FROM_MNP"] = "##INFO=<ID=FROM_MNP,Number=0,Type=Flag,Description=" + text;
  }

  // If find hets is specified, check that genotypes exist.
  //if (findHets && !v.hasGenotypes) {
  //  cerr << "Input files does not contain genotype information." << endl;
  //  cerr << "Cannot identify sample genotypes." << endl;
  //  exit(1);
  //}

  // Add an information line to the header describing the HET info.
  //if (findHets) {v.headerInfoLine["HET"] = "##INFO=<ID=HET,Number=.,Type=String,Description=\"List of samples heterozygous at this locus.\">";}

  // If variants that are found in the list of provided samples are to be outputted,
  // check that genotypes exist and that the provided file contains at least one 
  // sample.  Also check that the samples provided exist in the vcf file.
  //if (useSampleList) {
  //  samples sl; // Create a samples list object.
  //  sl.openSamplesFile(samplesListFile);
  //  sl.getSamples(v);
  //}

  // If the genotypes are to be removed, set the removeGenotypes value to 
  // true for the vcf object.
  if (removeGenotypes) {var.removeGenotypes = true;}

  // If records are to be stripped out of the vcf file, check the inputted
  // IDs and populate the list of IDs to be stripped.
  //if (stripRecords) {stripInfoList = checkInfoFields(v, stripInfo);}

  // If records are to be kept based on the contents of the info fields,
  // check the inputted IDs and populate the list of IDs to be kept.
  //if (keepRecords) {keepInfoList = checkInfoFields(v, keepInfoFields);}

  // If info fields are to be deleted, check the inputted IDs exist in the
  // header and populate the list of IDs to be kept.
  //if (removeInfo) {removeInfoList = checkInfoFields(v, removeInfoString);}

  // If dbSNP info is to be cleared, add the dbSNP tags to the remove info
  // list if it exists.  If not cretae it.
  //if (cleardbSnp) {
  //  removeInfo = true;
  //  removeInfoList.push_back("dbSNP");
  //  removeInfoList.push_back("dbSNPX");
  //  removeInfoList.push_back("dbSNPM");
  //}

  // Write out the header.
  header.writeHeader(ofile.outputStream, removeGenotypes, taskDescription);

// Read through all the entries in the file.  First construct the
// structure to contain the variants in memory and populate.

// If the genotypes need to be found, determine which position in the genotype
// format string the genotypes appear.
//  if (findHets || useSampleList) {
//    unsigned int count = 0;
//    if (!v.hasGenotypes) {
//      cerr << "Genotypes must be present to pick out variants occuring in a list of samples." << endl;
//      exit(1);
//    }
//    vector<string> s = split(v.variantRecord.genotypeFormatString, ":");
//    for (vector<string>::iterator iter = s.begin(); iter != s.end(); iter++) {
//      if (*iter == "GT") {break;}
//      count++;
//    }
//    genotypePosition = count;
//  }

  // Get the first record from the vcf file.
  v.success = v.getRecord();
  while (v.success) {

    // Build the variant structure for this reference sequence.
    if (var.originalVariantsMap.size() == 0) {
      currentReferenceSequence = v.variantRecord.referenceSequence;
      v.success                = var.buildVariantStructure(v);
    }

    // Loop over the variant structure until it is empty.  While v.update is true,
    // i.e. when the reference sequence is still the current reference sequence,
    // keep adding variants to the structure.
    while (var.originalVariantsMap.size() != 0) {
      if (v.variantRecord.referenceSequence == currentReferenceSequence && v.success) {
        var.addVariantToStructure(v.position, v.variantRecord);
        v.success = v.getRecord();
      }
      var.ovmIter = var.originalVariantsMap.begin();

      // Perform all filtering tasks on this variant.
      filter(var);
      var.buildOutputRecord(ofile, header);
      var.originalVariantsMap.erase(var.ovmIter);
    }
  }

  // Close the vcf files.
  v.closeVcf();

  // Flush the output buffer.
  ofile.flushOutputBuffer();

  return 0;
}
