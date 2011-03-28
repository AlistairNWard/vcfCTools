// ******************************************************
// vcfCTools (c) 2011 Alistair Ward
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ------------------------------------------------------
// Last modified: 18 February 2011
// ------------------------------------------------------
// vcfClass describes the vcf class and all operations.
// ******************************************************

#include "vcf.h"

#include <iostream>
#include <string>
#include <fstream>
#include <cstdlib>
#include <ctype.h>
#include <vector>

using namespace std;
using namespace vcfCTools;

// Constructor.
vcf::vcf(void) {
  hasHeader = true;
  processInfo = false;
  hasGenotypes = true;
  processGenotypes = false;
  numberDataSets = 0;
  dbsnpVcf = false;
}

// Destructor.
vcf::~vcf(void)
{}

// Open a vcf file.
bool vcf::openVcf(string filename) {
  vcfFilename = filename;
  if (vcfFilename != "-") {
    file.open(vcfFilename.c_str(), ifstream::in);
    input = &file;
    if (!file.is_open()) {
      cerr << "Failed to open file: " << vcfFilename << endl;
      exit(1);
    }
  }
  else {input = &cin;}
}

// Close the vcf file.
void vcf::closeVcf() {
  if (vcfFilename != "-") {
    file.close();
    if (file.is_open()) {
      cerr << "Failed to close file: " << vcfFilename << endl;
    }
  }
}

// Parse the vcf header.
void vcf::parseHeader() {
  bool success = true;
  while(getline(*input, headerLine)) {
    if (headerLine.substr(0, 6) == "##INFO") {success = headerInfo(headerLine, 0);}
    else if (headerLine.substr(0, 8) == "##FORMAT") {success = headerInfo(headerLine, 1);}
    else if (headerLine.substr(0, 6) == "##FILE") {success = headerFiles(headerLine);}
    else if (headerLine.substr(0, 2) == "##") {success = headerAdditionalInfo(headerLine);}
    else if (headerLine.substr(0, 1) == "#") {success = headerTitles(headerLine);}
    else {success = noHeader();}
    if (!success) {break;}
  }
}

// Parse information from the info and format descriptors.
bool vcf::headerInfo(string& headerLine, unsigned int headerLineType) {

// Break up the string and fond the info or format tag.
  size_t start = headerLine.find_first_of("<");
  size_t end = headerLine.find_first_of(">");
  string line = headerLine.substr(start, end - start + 1);

  // Find the tag, type and number of expected entries.
  size_t idPosition       = line.find("ID=");
  size_t numberPosition   = line.find("Number=");
  size_t typePosition     = line.find("Type=");
  size_t descPosition     = line.find("Description=");

  headerInfoStruct infoTag;
  string tag = line.substr(idPosition + 3, numberPosition - idPosition - 4);
  if (idPosition == string::npos || numberPosition == string::npos || 
      typePosition == string::npos || descPosition == string::npos) {
    infoTag.number  = 0;
    infoTag.type    = "unknown";
    infoTag.success = false;
  }
  else {
    string number = line.substr(numberPosition + 7, typePosition - numberPosition - 8);
    string type = line.substr(typePosition + 5, descPosition - typePosition - 6);
    string description = line.substr(descPosition + 12, line.size() - descPosition - 11);
    infoTag.success = true;

// The number field can take the value ".", if the number of values is variable.  If
// the alternate alleles field is a comma separated list, this could force other
// fields to allow a comma separated list also.
    if (number == ".") {infoTag.number = 0;}
    else {
      if (isdigit(number[0])) {
        infoTag.number  = atoi(number.c_str());
      }
      else {
        infoTag.success = false;
        infoTag.number  = 0;
      }
    }
    infoTag.type = type;
  }

  if (headerLineType == 0) {
    headerInfoFields[tag] = infoTag;
    headerInfoLine[tag]   = headerLine;
  }
  else if (headerLineType == 1) {
    headerFormatFields[tag] = infoTag;
    headerFormatLine[tag] = headerLine;
  }

  return true;
}

// Check to see if the records contain information from multiple different
// sources.  If vcfPytools has been used to find the intersection or union
// of two vcf files, the records may have been merged to keep all the
// information available.  If this is the case, there will be a ##FILE line
// for each set of information in the file.  The order of these files needs
// to be maintained.
bool vcf::headerFiles(string& headerLine) {
  size_t found1 = headerLine.find("ID=");
  size_t found2 = headerLine.find(",File");
  size_t found3 = headerLine.find("\">");
  if (found1 == string::npos ||found1 == string::npos || found1 == string::npos) {
    cerr << "ERROR:" << endl;
    cerr << "Unable to resolve the header lines containing merged filenames (##FILE)." << endl;
    exit(1);
  }
  unsigned int fileID = atoi( (headerLine.substr(found1 + 3, found2 - found1 - 3)).c_str() );
  if (fileID == 0) {
    cerr << "ERROR:" << endl;
    cerr << "Header line describing files whose entries have been merged (##FILE)," << endl;
    cerr << "has a non-integer file ID." << endl;
    exit(1);
  }
  string filename = headerLine.substr(found2 + 7, found3 - found2 - 7);
  if (includedDataSets.count(fileID) != 0) {
    cerr << "ERROR: file " << vcfFilename << endl;
    cerr << "Multiple files in the ##FILE list have identical ID values." << endl;
    exit(1);
  }
  includedDataSets[fileID] = vcfFilename;

// Set the number of files with information in this vcf file.
  if (fileID > numberDataSets) {numberDataSets = fileID;}

  return true;
}

// Parse additional information from the header.
bool vcf::headerAdditionalInfo(string& headerLine) {
  if (headerText == "") {headerText = headerLine;}
  else {headerText += "\n" + headerLine;}

  return true;
}

// Parse the header titles containing sample names if genotypes are present.
bool vcf::headerTitles(string& headerLine) {
  headerTitlesText = headerLine;

  vector<string> titles = split(headerLine,'\t');
  if (titles.size() > 8) {
    hasGenotypes = true;
    for (int i = 9; i < titles.size(); i++) {samples.push_back(titles[i]);}
    numberSamples = samples.size();
  }
  else if (titles.size() == 8) {
    hasGenotypes = false;
    processGenotypes = false;
    numberSamples = 0;
  }
  else {
    cerr << "Not all vcf standard fields are available." << endl;
    cerr << endl;
    cerr << "Got:     " << headerLine << endl;
    cerr << "Require: #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO" << endl;
    exit(1);
  }

  return false;
}

// No header is present.
bool vcf::noHeader() {

// Check if any of the other header strings have been populated.  It is possible
// that there is header information, but the title line (#CHROM POS ...) is missing
// or malformed.  If this is the case, terminate with a warning.

  if (!headerText.empty() || headerInfoFields.size() != 0 || headerFormatFields.size() != 0 || includedDataSets.size() != 0) {
    cerr << "No titles line in the vcf file (#CHROM POS etc.)" << endl;
    cerr << "vcfCTools requires a properly constructed header or no header at all to ensure correct operation." << endl;
    cerr << "Please add the required title line to the vcf file:" << endl;
    cerr << "	" << vcfFilename << endl;
    exit(1);
  }

  hasHeader = false;

  closeVcf();
  openVcf(vcfFilename);

  return false;
}

// Get the next record from the vcf file.
bool vcf::getRecord() {
  bool success = true;
  success = getline(*input, record);

// Return false if no more records remain.
  if (!success) {return false;}

  vector<string> recordFields;
  if (processGenotypes) {recordFields = split(record, '\t');}
  else {recordFields = split(record, '\t', 9);}

// Populate the variant values.
  referenceSequence = recordFields[0];
  position = atoi(recordFields[1].c_str());
  rsid     = recordFields[2];
  ref      = recordFields[3];
  alt      = recordFields[4];
  sQuality = recordFields[5];
  quality  = atof(recordFields[5].c_str());
  filters  = recordFields[6];
  info     = recordFields[7];
  if (recordFields.size() < 9) {hasGenotypes = false;}
  else {genotypeFormatString = recordFields[8];}

// If the position is not an integer, the conversion to an integer will have
// failed and position = 0.  In this case, terminate with an error.
  //if (position == 0 || quality == 0) {
  if (position == 0) {
    cerr << "Error processing record." << endl;
    if (position == 0) {cerr << "Variant position is not an integer." << endl;}
    if (quality == 0) {cerr << "Variant quality is not an integer or a floating point number." << endl;}
    cerr << endl;
    cerr << "Record:" << endl << endl << record << endl;
    exit(1);
  }

// Determine the variant type and whether or not there are multiple alternate
// alleles.
  size_t found = alt.find(",");
  if (found != string::npos) {hasMultipleAlternates = true;}
  else {
    hasMultipleAlternates = false;
    isSNP = (ref.size() == 1 && (ref.size() - alt.size()) == 0) ? true : false;
    isMNP = (ref.size() != 1 && (ref.size() - alt.size()) == 0) ? true : false;
    isDeletion  = ( (ref.size() - alt.size()) > 0) ? true : false;
    isInsertion = ( (alt.size() - ref.size()) > 0) ? true : false;
  }

// Add the reference sequence to the map.  If it didn't previously
// exist append the reference sequence to the end of the list as well. 
// This ensures that the order in which the reference sequences appeared
// in the header can be preserved.
  if (referenceSequences.count(referenceSequence) == 0) {
    referenceSequences[referenceSequence] = true;
    referenceSequenceVector.push_back(referenceSequence);
  }

// If required, parse the genotype format string and create a vector
// containing all of the individual sample genotype strings.
  if (processGenotypes) {
    genotypeFormat = split(genotypeFormatString, ":");

// Check that the number of genotype fields is equal to the number of samples
    genotypes = recordFields;
    genotypes.erase(genotypes.begin(), genotypes.begin() + 9);
    if (samples.size() != genotypes.size()) {
      cerr << "Error processing genotypes." << endl;
      cerr << "The number or genotypes (" << genotypes.size() << ") is not equal to the number of samples (" << samples.size() << ")." << endl;
      exit(1);
    }
  }

// If required, process the info fields.
  if (processInfo) {processInfoFields();}

  return true;
}

bool vcf::getVariantGroup(variantGroup& vg) {
  bool success = true;
  bool inGroup = true;
  vector<string> alts;
  unsigned int compStart;
  unsigned int compEndMatch;
  string::iterator refIter;
  string::reverse_iterator refRevIter;

// Define the variant group structure.
  vg.clear();
  vg.ref = ref;

// Determine the alt allele types (SNPs, MNPs, indels) and compare to
// the reference sequence.
  if (hasMultipleAlternates) {alts = split(alt, ",");}
  else {alts.push_back(alt);}

  cout << ref << endl;
  for (vector<string>::iterator iter = alts.begin(); iter != alts.end(); iter++) {
    compStart = 0;
    compEndMatch = 0;
    if ((*iter).length() == ref.length()) {

// Insertions.
    } else if ((*iter).length() > ref.length()) {
      isInsertion = true;

// Match the alt to the reference and find where the insertion occurs.
      refIter = ref.begin();
      for (string::iterator altIter = (*iter).begin(); altIter != (*iter).end(); altIter++) {
        if ((*altIter) != (*refIter)) {
          break;
        }
        else {
          refIter++;
          compStart++;
        }
      }

// Deletions.
    } else if ((*iter).length() < ref.length()) {
      isDeletion = true;

// Match the alt to the reference and find where the deletion occurs.
      refIter = ref.begin();
      for (string::iterator altIter = (*iter).begin(); altIter != (*iter).end(); altIter++) {
        if ((*altIter) != (*refIter)) {
          refRevIter = ref.rbegin();
          for (string::reverse_iterator altRevIter = (*iter).rbegin(); altRevIter != (*iter).rend(); altRevIter++) {
            if ((*altRevIter) != (*refRevIter)) {
              break;
            } else {
              refRevIter++;
              compEndMatch++;
            }
          }
          break;
        }
        else {
          refIter++;
          compStart++;
        }
      }
    }
    if (compEndMatch == 0) {
      cout << ref.substr(0, (*iter).length()) << endl;;
    } else {
      if (isInsertion) {
      } else if (isDeletion) {
        string a = ref.substr(0, compStart);
        a.append(ref.length() - (*iter).length(), '-');
        a.append((*iter).substr((*iter).length() - compEndMatch, compEndMatch));
        cout << a << endl;
      }
    }
  }
  cout << endl;

// Define the start and end position of the variant cluster from the first
// record and clear other entries.
  vg.referenceSequence = referenceSequence;
  vg.start = position;
  vg.ref = ref;
  vg.end = position + ref.size() - 1;

// Loop over all records that contain variants that overlap the reference
// sequence.  The reference sequence can be extended with each new variant,
// but stop looping if the end of the file is reached or the start position
// of the next variant reference sequence does not overlap the cluster region.
  while (success && inGroup) {
    if (referenceSequence == vg.referenceSequence && position <= vg.end) {
      vg.noRecords++;

// Modify the reference sequence.
      if ( (position + ref.size() - 1) > vg.end) {
        vg.ref += ref.substr( (vg.end - position + 1), (ref.size() - vg.end + position - 1) );
        vg.end = position + ref.size() - 1;
      }

      vg.noAlts = (hasMultipleAlternates) ? vg.noAlts + alts.size() : vg.noAlts + 1;
      success = getRecord();
    } else {
      inGroup = false;
    }
  }
  vg.noGroups++;

  return success;
}

// Process the info entries.
void vcf::processInfoFields() {
  infoTags.clear();
  vector<string> infoEntries = split(info, ';');
  for (vector<string>::iterator iter = infoEntries.begin(); iter != infoEntries.end(); iter++) {
    size_t found = (*iter).find_first_of("=");
    string tag = (*iter).substr(0, found);
    if (found == string::npos) {infoTags[tag] = "0";}
    else {infoTags[tag] = (*iter).substr(found + 1, (*iter).size());}
  }
}

// Parse the genotype format string for this record and put the information into
// genotypeFormat.  Create a vector containing all the genotype information.
// There should be as many genotypes as samples.  Actual processing of the
// individual sample genotype strings is performed when required by
// vcf::getGenotypeInfo.
void vcf::processGenotypeFields(string& genotypeString) {
  genotypeTags.clear();
  vector<string> genotypeEntries = split(genotypeString, ':');

// Check that the genotype string has the correct number of entries.
  if (genotypeString == "./." || genotypeString == ".") {}
  else {
    if (genotypeEntries.size() != genotypeFormat.size()) {
      cerr << "Error processing genotypes." << endl;
      cerr << "The genotype string does not contain the correct number of entries:" << endl;
      cerr << "\tFormat  : " << genotypeFormatString << endl;
      cerr << "\tGenotype: " << genotypeString << endl;
      exit(1);
    }
    vector<string>::iterator formatIter = genotypeFormat.begin();
    for (vector<string>::iterator iter = genotypeEntries.begin(); iter != genotypeEntries.end(); iter++) {
      genotypeTags[*formatIter] = *iter;
      formatIter++;
    }
  }
}

// Get the information for a specific info tag.  Also check that it contains
// the correct number and type of entries.
information vcf::getInfo(string& tag) {
  information sInfo;
  sInfo.tag = tag;

// If this routine has been called and processInfo is set to false, terminate the
// program.  Information can only be retrieved if the info fields have been
// processed and so entering this routine without having procesed the info fields
// will results in a failue to extract information.
  if (!processInfo) {
    cerr << "Routine vcf::getInfo called while processInfo = false." << endl;
    cerr << "The tool calling this routine must set processInfo = true for this object." << endl;
    cerr << "If processInfo = false, the info fields are not interrogated (in order to save time)," << endl;
    cerr << "but the getInfo routine is useless in this instance as required data structures" << endl;
    cerr << "have not been populated." << endl;
    cerr << endl;
    cerr << "Please check the logic of the called tool." << endl;
    exit(1);
  }

// Check if the tag exists in the header information.  If so,
// determine the number and type of entries asscoiated with this
// tag.
  if (headerInfoFields.count(tag) > 0) {
    sInfo.number = headerInfoFields[tag].number;
    sInfo.type   = headerInfoFields[tag].type;

// First check that the tag exists in the information string.  Then split
// the entry on commas.  For flag entries, do not perform the split.
    if (infoTags.count(tag) > 0) {
      if (sInfo.number == 0 && sInfo.type == "Flag") {sInfo.values.push_back("true");}
      else if (sInfo.number != 0 && sInfo.type == "Flag") {
        cerr << "Error processing info string." << endl;
        cerr << "Header inforamtion for entry: " << tag << " lists a flag with a non-zero number of entries." << endl;
        exit(1);
      }
      else {
        sInfo.values = split(infoTags[tag],",");
        if (sInfo.values.size() != sInfo.number) {
          cerr << "Error processing info string." << endl;
          cerr << "Unexpected number of entries for info field " << tag << " at " << referenceSequence << ":" << position << endl;
          exit(1);
        }
      }
    }
    else {sInfo.number = 0;}
  }
  else {
    cerr << "Error processing info string." << endl;
    cerr << "No information in the header for info entry: " << tag << endl;
    exit(1);
  }

  return sInfo;
}

// Get the genotype information.
information vcf::getGenotypeInfo(string& tag) {
  information gInfo;
  if (headerFormatFields.count(tag) > 0) {
    gInfo.number = headerFormatFields[tag].number;
    gInfo.type   = headerFormatFields[tag].type;

    gInfo.values = split(genotypeTags[tag], ",");
    if (gInfo.values.size() != gInfo.number) {
      cerr << "Error processing info string." << endl;
      cerr << "Unexpected number of entries for genotype entry " << tag << " at " << referenceSequence << ":" << position << endl;
      exit(1);
    }
  }
  else {
    cerr << "Error processing info string." << endl;
    cerr << "No information in the header for genotype format entry: " << tag << endl;
    exit(1);
  }

  return gInfo;
}

// Parse through the vcf file until the correct reference sequence is
// encountered and the position is greater than or equal to that requested.
bool vcf::parseVcf(string& compReferenceSequence, unsigned int compPosition, bool write, ostream* output) {
  bool success = true;
  if (referenceSequence != compReferenceSequence) {
    while (referenceSequence != compReferenceSequence and success) {
      if (write) {*output << record << endl;}
      success = getRecord();
    }
  }
  while ( (referenceSequence == compReferenceSequence) && (position < compPosition) && success) {
    if (write) {*output << record << endl;}
    success = getRecord();
  }

  return success;
}

// Parse through the vcf file until the correct reference sequence is
// encountered and then construct groups of variants occupying
// overlapping reference sequence and parse through these until the
// start position of the cluster is greater than or equal to the
// requested value.
bool vcf::parseVcfGroups(variantGroup& vc, string& compReferenceSequence, unsigned int compPosition, bool write, ostream* output) {
  bool success = true;
  if (vc.referenceSequence != compReferenceSequence) {
    while (referenceSequence != compReferenceSequence && success) {
      //if (write) {*output << record << endl;}
      success = getRecord();
    }
  }
  if (success) {success = getVariantGroup(vc);}
  while ( (vc.referenceSequence == compReferenceSequence) && (vc.end < compPosition) && success) {
    //if (write) {*output << record << endl;}
    success = getVariantGroup(vc);
  }

  return success;
}

// Construct a vcf record.
string vcf::buildRecord(bool removeGenotypes) {
  ostringstream sPosition;
  sPosition << position;
  string build= referenceSequence + "\t" + \
                sPosition.str() + "\t" + \
                rsid + "\t" + \
                ref + "\t" + \
                alt + "\t" + \
                sQuality + "\t" + \
                filters + "\t" + \
                info;

  if (hasGenotypes && !removeGenotypes) {
    size_t found = record.find(genotypeFormatString);
    build += "\t" + record.substr(found);
  }

  return build;
}
