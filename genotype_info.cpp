// ******************************************************
// vcfCTools (c) 2011 Alistair Ward
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ------------------------------------------------------
// Last modified: 18 February 2011
// ------------------------------------------------------
// Tools for manipulating the sample genotypes
// ******************************************************

#include "genotype_info.h"
#include "split.h"

using namespace std;
using namespace vcfCTools;

// Constructor
genotypeInfo::genotypeInfo(string format, string gen, map<string, headerInfoStruct>& header) {
  genotypeFormat = format;
  genotypeString = gen;
  headerInfo     = header;
};

// Desctructor.
genotypeInfo::~genotypeInfo(void) {};

// Modify the genotypes to reflect modifications to the alleles
// present in this record.
void genotypeInfo::modifyGenotypes(vector<int>& alleleIDs) {
  int gtInt;
  string modifiedGenotype;
  string modifiedGenotypeString;
  string modifiedGT;
  string modifiedField;
  vector<string> originalGenotypes;
  vector<string> fields;
  vector<string>::iterator fieldsIter;
  vector<int>::iterator alleleIter;

  // Split the genotype format and genotype string.
  genotypeStrings = split(genotypeString, "\t");
  formats         = split(genotypeFormat, ":");

  // Loop through each sample modifying the genotype as necessary.
  vector<string>::iterator iter = genotypeStrings.begin();
  for (; iter != genotypeStrings.end(); iter++) {
    modifiedGenotype = "";

    // Check that the number of fields in the genotype string matches
    // the number of fields specified in the format.  If there is
    // only one entry in the genotypes, this indicates that the
    // genotypes couldn't be determined and thus the entry is '.'.  In
    // this case, the genotypes can be left as is.
    if (*iter == ".") {
      modifiedGenotype = *iter;
    } else {
      originalGenotypes = split(*iter, ":");
      if (originalGenotypes.size() != formats.size()) {
        cerr << "ERROR: Number of entries in genotype string does not match the format." << endl;
        exit(1);
      } else {

        // Loop through the entries in the genotype format string.  If
        // the number of entries is an integer, this field does not need
        // to be modified.  However if the number of entries is A or G,
        // then these fields will require modification.  The GT (genotype)
        // field will also require modification.
        vector<string>::iterator fIter = formats.begin();
        vector<string>::iterator gIter = originalGenotypes.begin();
        for (; fIter != formats.end(); fIter++) {
          if (*fIter == "GT") {

            // Loop over the allele identifiers and replace them with the
            // contents of the imported array alleleIDs.
            vector<string> originalGT       = split(*gIter, "/");
            vector<string>::iterator gtIter = originalGT.begin();
            modifiedGT = "";
            for (; gtIter != originalGT.end(); gtIter++) {
              gtInt = atoi((*gtIter).c_str());
              ostringstream gt;
              gt << alleleIDs[gtInt];
              modifiedGT = (alleleIDs[gtInt] == -1) ? (modifiedGT + "./") : (modifiedGT + gt.str() + "/");
            }

            // Strip off the extra '/' character.
            modifiedGT = modifiedGT.substr(0, modifiedGT.size() - 1);
            modifiedGenotype += modifiedGT + ":";

          // Only retain the fields for the retained alleles.
          } else if (headerInfo[*fIter].number == "A") {
            fields = split(*gIter, ",");
            fieldsIter = fields.begin();
            modifiedField = "";

            // Loop through the alleleIDs array (skipping the first -
            // reference - entry) and if the value is 0, this alterate
            // was skipped and so the value should no be included in the
            // new genotype string.
            alleleIter = alleleIDs.begin();
            alleleIter++;
            for (; alleleIter != alleleIDs.end(); alleleIter++) {
              if (*alleleIter != -1) {
                modifiedField += *fieldsIter + ",";
              }
              fieldsIter++;
            }

            // Strip off the extra '/' character.
            modifiedField = modifiedField.substr(0, modifiedField.size() - 1);
            modifiedGenotype += modifiedField + ":";

          } else if (headerInfo[*fIter].number == "G") {
            modifiedGenotype += ".:";
          } else if (headerInfo[*fIter].number == ".") {
            modifiedGenotype += *gIter + ":";
          } else {
            modifiedGenotype += *gIter + ":";
          }
          gIter++;
        }

        // Strip off the extra ':' character.
        modifiedGenotype = modifiedGenotype.substr(0, modifiedGenotype.size() - 1);
      }
    }
    modifiedGenotypeString += modifiedGenotype + "\t";
  }
  // Strip off the extra tab character.
  modifiedGenotypeString = modifiedGenotypeString.substr(0, modifiedGenotypeString.size() - 1);

  // Replace the original genotype string with the modified version.
  genotypeString = modifiedGenotypeString;
}
