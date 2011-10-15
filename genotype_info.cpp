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

using namespace std;
using namespace vcfCTools;

// Constructor
genotypeInfo::genotypeInfo(string format, string gen) {
  genotypeFormat  = format;
  genotypeString  = gen;
  genotypeFormats = split(genotypeFormat, ":");
  genotypes       = split(gen, "\t");
  genotypeFields.clear();
};

// Desctructor.
genotypeInfo::~genotypeInfo(void) {};

// Modify the genotypes to reflect modifications to the alleles
// present in this record.
void genotypeInfo::modifyGenotypes(vcfHeader& header, vector<int>& alleleIDs) {
  bool success;
  int gtInt;
  string connector;
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
            success = getAlleles();
            if (!success) {
              cerr << "ERROR: Genotype includes both phased and unphased alleles." << endl;
              cerr << "Validate the vcf file for further information." << endl;
              exit(1);
            } else {
              connector = (phased) ? "|" : "/";
            }
            vector<string>::iterator gtIter = originalAlleles.begin();
            modifiedGT = "";
            for (; gtIter != originalAlleles.end(); gtIter++) {
              gtInt = atoi((*gtIter).c_str());
              ostringstream gt;
              gt << alleleIDs[gtInt];
              modifiedGT = (alleleIDs[gtInt] == -1) ? (modifiedGT + "." + connector) : (modifiedGT + gt.str() + connector);
            }

            // Strip off the extra '/' character.
            modifiedGT = modifiedGT.substr(0, modifiedGT.size() - 1);
            modifiedGenotype += modifiedGT + ":";

          // Only retain the fields for the retained alleles.
          } else if (header.formatFields[*fIter].number == "A") {
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

          } else if (header.formatFields[*fIter].number == "G") {
            modifiedGenotype += ".:";
          } else if (header.formatFields[*fIter].number == ".") {
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

// Parse all of the sample information and determine whether this person is
// homozygous reference, heterozygous, homozygous non-reference or unknown.
void genotypeInfo::processFormats(vcfHeader& header) {
  unsigned int formatID = 0;
  genotypeFormats = split(genotypeFormat, ":");
  vector<string>::iterator formatIter = genotypeFormats.begin();
  genoStruct geno;

  for (; formatIter != genotypeFormats.end(); formatIter++) {

    // In order to process the info field, the header information is
    // required.  If this does not exist, terminate the program.
    if (header.formatFields.count(*formatIter) == 0) {
      cerr << "ERROR: No header description for format tag: " << *formatIter << endl;
      exit(1);
    }

    geno.number = header.formatFields[*formatIter].number;
    geno.type   = header.formatFields[*formatIter].type;
    geno.ID     = formatID;
    genotypeFields[*formatIter] = geno;
    formatID++;
  }
}

// Extract information about a specified tag.
void genotypeInfo::validateGenotypes(vcfHeader& header, string& referenceSequence, int& position, unsigned int& noAlts, bool& error) {
  bool success;
  int integerValue;
  unsigned int sampleID = 0;

  // Read the format string and check that all of the entries are defined.
  processFormats(header);

  // Loop through each of the samples and check that the genotype information is
  // consistent with the alternate alleles and the format string.
  genotypeStrings               = split(genotypeString, "\t");
  vector<string>::iterator iter = genotypeStrings.begin();
  for (; iter != genotypeStrings.end(); iter++) {

    // If no information is known about the genotypes for this sample, the string 
    // will simply be '.'.  In this case, there is no need to interrogate the
    // string.
    if (*iter != ".") {
      genotypes = split(*iter, ":");

      // First check if there are the correct number of entries. There should be as
      // many fields (seperated by a ":") as there are formats.
      if (genotypeFields.size() != genotypes.size()) {
        cerr << "ERROR: Inconsistent number of fields in the genotype string for sample: " << header.samples[sampleID];
        cerr << " at " << referenceSequence << ":" << position << "." << endl;
        error = true;
      }

      // Determine the number of alleles in the genotype.  In order to do this check
      // to see if a '/' is present in the genotype.  If so, split on '/', otherwise
      // split on '|'.  For the purposes of validation, check that only '/' or '|' is
      // present and not both.
      success = getAlleles();
      if (!success) {
        cerr << "ERROR: Genotype includes both phased and unphased alleles for sample: " << header.samples[sampleID];
        cerr << " at " << referenceSequence << ":" << position << "." << endl;
        error = true;
      }

      // Parse the individual entries and ensure that they are consistent with the
      // format.
      genoIter = genotypeFields.begin();
      for (; genoIter != genotypeFields.end(); genoIter++) {

        // Determine the number of values for this format field.
        values = split(genotypes[genotypeFields[genoIter->first].ID], ",");

        // The genotypes themselves need to be treated seperately.
        if (genoIter->first == "GT") {

          vector<string>::iterator aIter = originalAlleles.begin();
          for (; aIter != originalAlleles.end(); aIter++) {

            // First check that all the entries in the genotypes are integers or '.'.
            integerValue = atoi( (*aIter).c_str() );
            if (integerValue == 0 && *aIter != "0" && *aIter != ".") {
              cerr << "ERROR: Invalid entry in genotype for sample: " << header.samples[sampleID] << " at ";
              cerr << referenceSequence << ":" << position << "." << endl;
              error = true;

            // Next check that the value is not larger than the number of alternate
            // alleles.  The numbers in the genotypes refer to the allele, 0
            // represents the reference, 1 the first alternate allele etc.  If there
            // are 3 alternate alleles, the genotype can only contain integers in the
            // range [0-3].
            } else {
              if (integerValue > noAlts) {
                cerr << "ERROR: Genotype entry for sample: " << header.samples[sampleID] << " at ";
                cerr << referenceSequence << ":" << position << " represents a non-existent allele (max value: ";
                cerr << noAlts << ")." << endl;
                error = true;
              }
            }
          }

        // Now check all the other fields in the genotype string.
        } else {
 
          // First check that the number of elements in each field is consistent
          // with the expected number.
          if (genoIter->second.number == "A") {
            if (values.size() != noAlts) {
              cerr << "ERROR: Incorrect number of entries for sample " << header.samples[sampleID];
              cerr << " at " << referenceSequence << ":" << position << " in " << genoIter->first << " field" << endl;
              error = true;
            }
  
            // Check the values conform to the expected types.
            checkTypes(referenceSequence, position, error);

          // The number of fields is expected to be equal to the number of alleles in 
          // the genotypes.
          } else if (genoIter->second.number == "G") {

            // Now check that the correct number is observed.  The number of fields
            // depends on the number of possible alleles (i.e. the number of alternate
            // alleles plus the reference allele) and the ploidy.  For N alleles and
            // a plody of k, the number of values should be (N+k-1)!/k!(N-1)!.
            unsigned int A         = noAlts + numberInGeno;
            unsigned int num       = fact(A);
            unsigned int denA      = fact(numberInGeno);
            unsigned int denB      = fact(noAlts);
            unsigned int numberInG = num / (denA * denB);
            if (values.size() != numberInG) {
              cerr << "ERROR: Incorrect number of entries for sample: " << header.samples[sampleID];
              cerr << " at " << referenceSequence << ":" << position << " in " << genoIter->first << " field" << endl;
              error = true;
            }

            // Check the values conform to the expected types.
            checkTypes(referenceSequence, position, error);

          // If the info number = '.', any number of entries is acceptable, so do not
          // check the number of values.
          } else if (genoIter->second.number == ".") {

            // Check the values conform to the expected types.
            checkTypes(referenceSequence, position, error);

          // Finally, if a specified number of entries is given, check that this number
          // is observed.  Also check that the entry is itself an integer as there are no
          // other allowed types.
          } else {

            // First check that the number is an integer.
            integerValue = atoi( genoIter->second.number.c_str() );
            if (integerValue == 0) {
              cerr << "ERROR: Invalid type in genotype field: " << genoIter->first << " at " << referenceSequence;
              cerr << ":" << position << "." << endl;
              error = true;
            }
      
            // Now check that the correct number is observed.
            if (values.size() != integerValue) {
              cerr << "ERROR: Incorrect number of entries for sample " << header.samples[sampleID];
              cerr << " at " << referenceSequence << ":" << position << " in " << genoIter->first << " field" << endl;
              error = true;
            }
      
            // Check the values conform to the expected types.
            checkTypes(referenceSequence, position, error);
          }
        }
      }
    }
    sampleID++;
  }
}

// get the alleles in the genotypes.  These should be a set of integers representing the
// alleles listed in the alt allele column.  The genotypes should be phased or unphased.
bool genotypeInfo::getAlleles() {
  size_t count;
  phased       = true;
  unphased     = true;
  bool success = true;

  // Check for unphased genotypes.
  count = genotypes[genotypeFields["GT"].ID].find("/");
  if (count == string::npos) {unphased = false;}

  // Check for phased genotypes.
  count = genotypes[genotypeFields["GT"].ID].find("|");
  if (count == string::npos) {phased = false;}

  if (phased && unphased) {
    success = false;
  } else if (phased) {
    originalAlleles = split(genotypes[genotypeFields["GT"].ID], "|");
  } else if (unphased) {
    originalAlleles = split(genotypes[genotypeFields["GT"].ID], "/");
  } else {
    originalAlleles.clear();
    originalAlleles.push_back(genotypes[genotypeFields["GT"].ID]);
  }

  // Determine the number of alleles in the genotype.
  numberInGeno = originalAlleles.size();

  return success;
}

// Check that the observed values conform to the expected type.
void genotypeInfo::checkTypes(string& referenceSequence, int& position, bool& error) {
  int integerValue;
  float floatValue;
  vector<string>::iterator gIter;

  // Now check that the values are of the correct type.
  gIter = values.begin();
  for (; gIter != values.end(); gIter++) {

    // Check that integers appear if integers are expected.
    if (genoIter->second.type == "Integer") {
      integerValue = atoi( gIter->c_str() );
      if (integerValue == 0 && *gIter != "0") {
        cerr << "ERROR: Expected integer in genotype field: " << genoIter->first << " at " << referenceSequence;
        cerr << ":" << position << "." << endl;
        error = true;
      }

    // If floating point values are expected, check that this is seen.
    } else if (genoIter->second.type == "Float") {
      floatValue = atof( gIter->c_str() );
      if (floatValue == 0. && (*gIter != "0" && *gIter != "0." && *gIter != "0.0")) {
        cerr << "ERROR: Expected floats in genotype field: " << genoIter->first << " at " << referenceSequence;
        cerr << ":" << position << "." << endl;
        error = true;
      }
    } else if (genoIter->second.type == "Character") {

      // Check that the values are of length 1.
      if (gIter->length() != 1) {
        cerr << "ERROR: Expected a single character for genotype field: " << genoIter->first << "at";
        cerr << referenceSequence << ":" << position << "." << endl;
        error = true;
      }
    } else if (genoIter->second.type == "String") {
    }
  }
}
