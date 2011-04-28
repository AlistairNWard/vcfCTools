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
variantInfo::variantInfo(void) {};

// Desctructo.
variantInfo::~variantInfo(void) {};

// Split the info entries into its constituent parts and determine
// all the tags present..
void variantInfo::processInfoFields(string& infoString) {
  infoTags.clear();
  vector<string> infoEntries = split(infoString, ';');
  for (vector<string>::iterator iter = infoEntries.begin(); iter != infoEntries.end(); iter++) {
    size_t found = (*iter).find_first_of("=");
    string tag = (*iter).substr(0, found);
    if (found == string::npos) {infoTags[tag] = "flag";}
    else {infoTags[tag] = (*iter).substr(found + 1, (*iter).size());}
  }
}

//
void variantInfo::getInfo(string tag, string& ref, int position) {

// Get the information for a specific info tag.  Also check that it contains
// the correct number and type of entries.

// If this routine has been called and processInfoFields has no, or there are
// no fields in the info string, terminate the program.  Information can only
// be retrieved if the info fields have been processed and so entering this
// routine without having procesed the info fields will result in a failue toi
// extract information.
  if (infoTags.size() == 0) {
    cerr << "Routine vcf::getInfo called while there is no information in the" << endl;
    cerr << "info string, or the routine processInfoFields has not previously been" << endl;
    cerr << "called.  Please ensure that information exists in the vcf file" << endl;
    cerr << "and that processInfoFields has been called." << endl;
    cerr << endl;
    exit(1);
  }

// First check that the tag exists in the information string.  Then split
// the entry on commas.  For flag entries, do not perform the split.
  if (infoTags.count(tag) > 0) {

    // Check if the entry is a flag.
    if (infoTags[tag] == "flag") {
      type = "flag";
      values.push_back("true");
    } else {

      // Check if the info field is a comma separated list.
      type = "value";
      values = split(infoTags[tag], ",");
    }

  // Requested info tag is not in the info string
  } else {
    cerr << "Tag: " << tag << " is not present in the info string." << endl;
    cerr << "Coordinates: " << ref << ":" << position << endl;
    cerr << "Terminating program." << endl;
    exit(1);
  }
}
