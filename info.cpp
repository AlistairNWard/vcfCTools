// ******************************************************
// vcfCTools (c) 2011 Alistair Ward
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ------------------------------------------------------
// Last modified: 18 February 2011
// ------------------------------------------------------
// Tools for manipulating the information string.
// ******************************************************

#include "info.h"
#include "split.h"

using namespace std;
using namespace vcfCTools;

// Constructor
variantInfo::variantInfo(string& info, map<string, headerInfoStruct>& header) {
  infoString = info;
  headerInfo = header;
};

// Desctructor.
variantInfo::~variantInfo(void) {};

// If alleles have been removed, modify the info field so that all
// fields have the correct number of entries.
void variantInfo::modifyInfo(vector<int>& alleleIDs) {
  string modifiedInfo = "";

  retrieveFields();

  // Loop through each of the info fields and check for any fields that should
  // contain values corresponding to the different alternate alleles.
  infoIter = infoFields.begin();
  for (; infoIter != infoFields.end(); infoIter++) {
    if (infoIter->number == "A") {
      vector<string> values            = split(infoIter->value, ",");
      vector<string>::iterator vIter   = values.begin();
      vector<int>::iterator alleleIter = alleleIDs.begin();
      alleleIter++;

      modifiedInfo += infoIter->tag + "=";
      for (; alleleIter != alleleIDs.end(); alleleIter++) {
        if (*alleleIter != -1) {modifiedInfo += *vIter + ",";}
        vIter++;
      }

      // Strip the trailing ",".
      modifiedInfo = modifiedInfo.substr(0, modifiedInfo.size() - 1);
      modifiedInfo += ";";
    } else if (infoIter->number == "0") {
      modifiedInfo += infoIter->tag + ";";
    } else {
      modifiedInfo += infoIter->tag + "=" + infoIter->value + ";";
    }
  }

  // Strip the trailing ";" and then replace the existing info string
  // with the modified one.
  modifiedInfo = modifiedInfo.substr(0, modifiedInfo.size() - 1);
  infoString   = modifiedInfo;
}

// Split the info string into its components and populate the arrays
// to store the information.
void variantInfo::retrieveFields() {
  vector<string> infoArray = split(infoString, ";");
  vector<string>::iterator infoIter = infoArray.begin();
  for (; infoIter != infoArray.end(); infoIter++) {

    // Split the entry on "=" to get the tag.
    vector<string> fields = split(*infoIter, "=");
    infoStruct info;
    info.tag    = fields[0];

    // In order to process the info field, the header information is
    // required.  If this does not exist, terminate the program.
    if (headerInfo.count(info.tag) == 0) {
      cerr << "ERROR: No header description for info tag: " << info.tag << endl;
      exit(1);
    }

    info.number = headerInfo[info.tag].number;
    info.type   = headerInfo[info.tag].type;
    if (fields.size() > 1) {
      info.value          = fields[1];
      vector<string> temp = split(fields[1], ",");
      info.numberValues   = temp.size();
    } else {
      info.value        = "";
      info.numberValues = 0;
    }
    infoFields.push_back(info);
  }
}

// Split the info entries into its constituent parts and determine
// all the tags present..
//void variantInfo::processInfoFields(string& infoString) {
//  infoTags.clear();
//  vector<string> infoEntries = split(infoString, ';');
//  for (vector<string>::iterator iter = infoEntries.begin(); iter != infoEntries.end(); iter++) {
//    size_t found = (*iter).find_first_of("=");
//    string tag = (*iter).substr(0, found);
//    if (found == string::npos) {infoTags[tag] = "flag";}
//    else {infoTags[tag] = (*iter).substr(found + 1, (*iter).size());}
//  }
//}
//
////
//void variantInfo::getInfo(string tag, string& ref, int position) {
//
//// Get the information for a specific info tag.  Also check that it contains
//// the correct number and type of entries.
//
//// If this routine has been called and processInfoFields has no, or there are
//// no fields in the info string, terminate the program.  Information can only
//// be retrieved if the info fields have been processed and so entering this
//// routine without having procesed the info fields will result in a failue toi
//// extract information.
//  if (infoTags.size() == 0) {
//    cerr << "Routine vcf::getInfo called while there is no information in the" << endl;
//    cerr << "info string, or the routine processInfoFields has not previously been" << endl;
//    cerr << "called.  Please ensure that information exists in the vcf file" << endl;
//    cerr << "and that processInfoFields has been called." << endl;
//    cerr << endl;
//    exit(1);
//  }
//
//// First check that the tag exists in the information string.  Then split
//// the entry on commas.  For flag entries, do not perform the split.
//  if (infoTags.count(tag) > 0) {
//
//    // Check if the entry is a flag.
//    if (infoTags[tag] == "flag") {
//      type = "flag";
//      values.push_back("true");
//    } else {
//
//      // Check if the info field is a comma separated list.
//      type = "value";
//      values = split(infoTags[tag], ",");
//    }
//
//  // Requested info tag is not in the info string
//  } else {
//    cerr << "Tag: " << tag << " is not present in the info string." << endl;
//    cerr << "Coordinates: " << ref << ":" << position << endl;
//    cerr << "Terminating program." << endl;
//    exit(1);
//  }
//}
//
//// Construct the info fields for each alternate.  This may require splitting comma separated
//// info entries.
//vector<string> variantInfo::buildAltInfo(string& ref, int position, int numberAlts) {
//  vector<string> altInfoStrings;
//  altInfoStrings.resize(numberAlts);
//  vector<string>::iterator altIter;
//  string entry;
//
//  for (map<string, string>::iterator iter = infoTags.begin(); iter != infoTags.end(); iter++) {
//    values.clear();
//    getInfo(iter->first, ref, position);
//
//    // If the info field has a single entry (or is a flag), add this values to each alternate
//    // info string.
//    if (values.size() == 1) {
//      if (values[0] == "true") {entry = iter->first;}
//      else {entry = iter->first + "=" + values[0];}
//
//      // Add to each info string.
//      for (altIter = altInfoStrings.begin(); altIter != altInfoStrings.end(); altIter++) {
//        *altIter = (*altIter == "") ? *altIter + entry : *altIter + ";" + entry;
//      }
//
//    // If there is a comma seperated list and the number of entries is equal to the number
//    // stated in the header line, include all these entries in the new variant string.
//    } else if (values.size() == header[iter->first].number) {
//      entry = iter->first + "=" + infoTags[iter->first];
//
//      // Add to each info string.
//      for (altIter = altInfoStrings.begin(); altIter != altInfoStrings.end(); altIter++) {
//        *altIter = (*altIter == "") ? entry : *altIter + ";" + entry;
//      }
//
//    // If there is a comma separated set of values and the same number of entries as there
//    // are alternates, include successive values in each string.
//    } else if (values.size() == numberAlts) {
//      unsigned int entryNumber = 0;
//      for (altIter = altInfoStrings.begin(); altIter != altInfoStrings.end(); altIter++) {
//        entry = iter->first + "=" + values[entryNumber];
//        *altIter = (*altIter == "") ? *altIter + entry : *altIter + ";" + entry;
//        entryNumber++;
//      }
//
//    // If one of the above are true and the number of entries is variable (i.e. in the INFO
//    // description in the header is '.'), include the whole entry for each variant record.
//    } else if (header[iter->first].number == -1) {
//
//      // Add to each info string.
//      for (altIter = altInfoStrings.begin(); altIter != altInfoStrings.end(); altIter++) {
//        *altIter = (*altIter == "") ? *altIter + entry : *altIter + ";" + entry;
//      }
//
//    // There are multiple entries in the info field, but the number of entries does not
//    // correspond to the number of alternate alleles.  Terminate the program with an 
//    // error.
//    } else {
//      cerr << "ERROR: Incorrect number of fields in the " << iter->first << " info field. ";
//      cerr << " (" << ref << ":" << position << ")" << endl;
//      exit(1);
//    }
//  }
//
//  return altInfoStrings;
//}
