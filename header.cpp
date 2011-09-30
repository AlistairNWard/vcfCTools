// ******************************************************
// vcfCTools (c) 2011 Alistair Ward
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ------------------------------------------------------
// Last modified: 18 February 2011
// ------------------------------------------------------
// Tools for manipulating the vcf header.
// ******************************************************

#include "header.h"

#define ALT 0
#define FILTER 1
#define FORMAT 2
#define INFO 3

using namespace std;
using namespace vcfCTools;

// Constructor
vcfHeader::vcfHeader(void) {
  hasAssemblyInfo   = false;
  hasContigInfo     = false;
  hasPedigreeInfo   = false;
  hasPedigreeDbInfo = false;
  hasSampleInfo     = false;
  assembly          = "";
  fileFormat        = "";
};

// Destructor.
vcfHeader::~vcfHeader(void) {};

// Parse the vcf header.
void vcfHeader::parseHeader(istream* input) {
  while(getline(*input, line)) {
    if (line.substr(0, 5) == "##ALT")              {parseInfo(ALT);}
    else if (line.substr(0, 10) == "##assembly")   {parseAssembly();}
    else if (line.substr(0, 8)  == "##contig")     {parseContig();}
    else if (line.substr(0, 12) == "##fileformat") {parseFileFormat();}
    else if (line.substr(0, 8)  == "##FILTER")     {parseInfo(FILTER);}
    else if (line.substr(0, 8)  == "##FORMAT")     {parseInfo(FORMAT);}
    else if (line.substr(0, 6)  == "##INFO")       {parseInfo(INFO);}
    else if (line.substr(0, 10) == "##PEDIGREE")   {parsePedigree();}
    else if (line.substr(0, 12) == "##pedigreeDB") {parsePedigreeDb();}
    else if (line.substr(0, 8)  == "##SAMPLE")     {parseSample();}
    else if (line.substr(0, 2)  == "##")           {parseAdditionalInfo();}
    else if (line.substr(0, 1)  == "#") {
      parseTitles();
      break;
    }

    // A header is required, so if no header lines were found, terminate
    // with an error.
    else {
      cerr << "ERROR: No header lines present. At a minimum, the column descriptions are required." << endl;
      exit(1);
    }
  }
}

// Parse assembly information if present.
void vcfHeader::parseAssembly() {
  hasAssemblyInfo = true;
}

// Parse file format information.
void vcfHeader::parseFileFormat() {
}

// Parse contig information.
void vcfHeader::parseContig() {
  hasContigInfo = true;
}

// Parse information from the info and format descriptors.
void vcfHeader::parseInfo(unsigned int lineType) {
  string description;
  string number;
  string tag;
  string type;

// Break up the string and find the info or format tag.
  size_t start = line.find_first_of("<");
  size_t end   = line.find_first_of(">");
  string sLine = line.substr(start, end - start + 1);

  // Find the tag, type and number of expected entries.
  size_t idPosition       = sLine.find("ID=");
  size_t numberPosition   = sLine.find("Number=");
  size_t typePosition     = sLine.find("Type=");
  size_t descPosition     = sLine.find("Description=");

  headerInfo infoTag;
  infoTag.success = true;

  // The info and format fields include the type and number, while the
  // alt an filter fields do not.  Treat the two types separately.
  if (lineType == ALT || lineType == FILTER) {
    tag                = sLine.substr(idPosition + 3, descPosition - idPosition - 4);
    description        = sLine.substr(descPosition + 12, sLine.size() - descPosition - 13);
    infoTag.type       = ".";
    infoTag.number     = ".";
  } else if (lineType == INFO || lineType == FORMAT) {
    tag         = sLine.substr(idPosition + 3, numberPosition - idPosition - 4);
    number      = sLine.substr(numberPosition + 7, typePosition - numberPosition - 8);
    type        = sLine.substr(typePosition + 5, descPosition - typePosition - 6);
    description = sLine.substr(descPosition + 12, sLine.size() - descPosition - 13);

    // The number field can take an integer value if the number of values for this is
    // fixed.  If there is an entry per alternate, this should take the value 'A', one
    // entry per genotype takes the value 'G' and if the number is unknown or variable
    // then this should be '.'.
    if (number == "A") {
      infoTag.number = "A";
    } else if (number == "G") {
      infoTag.number = "G";
    } else if (number == ".") {
      infoTag.number = ".";
    } else if (isdigit(number[0])) {
      infoTag.number = number;
    } else {
      infoTag.success = false;
      infoTag.number  = ".";
      cerr << "WARNING: Malformed info/format in header for tag: ";
      cerr << tag << ".  Unknown value in Number field." << endl;
    }
    infoTag.type = type;
  }
  infoTag.description = description;

  // Depending on the header line type (for example, an info or format decriptor),
  // populate the correct header structure.
  if (lineType == ALT) {
    altFields[tag] = infoTag;
    altLine[tag]   = line;

  } else if (lineType == INFO) {
    infoFields[tag] = infoTag;
    infoLine[tag]   = line;

  } else if (lineType == FILTER) {
    filterFields[tag] = infoTag;
    filterLine[tag]   = line;

  } else if (lineType == FORMAT) {
    formatFields[tag] = infoTag;
    formatLine[tag]   = line;
  }
}

// Parse additional information from the header.
void vcfHeader::parseAdditionalInfo() {
  if (text == "") {text = line;}
  else {text += "\n" + line;}
}

// Parse pedigree information.
void vcfHeader::parsePedigree() {
  hasPedigreeInfo = true;
}

// Parse pedigree database information.
void vcfHeader::parsePedigreeDb() {
  hasPedigreeDbInfo = true;
}

// Parse sample information.
void vcfHeader::parseSample() {
  hasSampleInfo = true;
}

// Parse the header titles containing sample names if genotypes are present.
void vcfHeader::parseTitles() {
  vector<string> titles = split(line,'\t');
  numberSamples         = 0;
  if (titles.size() > 8) {
    for (int i = 9; i < titles.size(); i++) {samples.push_back(titles[i]);}
    numberSamples = samples.size();
  } else if (titles.size() < 8) {
    cerr << "ERROR: Not all vcf standard fields are present in the header.  Please check the header." << endl;
    exit(1);
  }
}

// Write out the header to the output.
void vcfHeader::writeHeader(ostream* output, bool removeGenotypes, string& description) {
  map<string, string>::iterator iter;

  // If the original vcf file contained header text and/or more text has
  // been added, write this out to the output.
  if (text != "") {*output << text << endl;}

  // Write out a description of the tasks performed by vcfCTools.
  *output << description << endl;

  // Write out all of the info, filter, format and symbolic alternate
  // allele  descriptions.
  for (iter = infoLine.begin(); iter != infoLine.end(); iter++) {*output << (*iter).second << endl;}
  for (iter = filterLine.begin(); iter != filterLine.end(); iter++) {*output << (*iter).second << endl;}
  for (iter = formatLine.begin(); iter != formatLine.end(); iter++) {*output << (*iter).second << endl;}
  for (iter = altLine.begin(); iter != altLine.end(); iter++) {*output << (*iter).second << endl;}

  // Write out the standard column headers.
  *output << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";

  // If the genotypes are not filtered out (assuming they exist in the input
  // vcf file), include the format header and all of the samples.
  if (!removeGenotypes && numberSamples != 0) {
    *output << "\tFORMAT";
    vector<string>::iterator sIter = samples.begin();
    for (; sIter != samples.end(); sIter++) {*output << "\t" << *sIter;}
  }
  *output << endl;
}

