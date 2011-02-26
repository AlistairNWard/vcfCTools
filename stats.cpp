// ******************************************************
// vcfCTools (c) 2011 Alistair Ward
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ------------------------------------------------------
// Last modified: 18 February 2011
// ------------------------------------------------------
// Tools for the statistics class.
// ******************************************************

#include "stats.h"

using namespace std;
using namespace vcfCTools;

// Constructor.
statistics::statistics(void) {
  lastSnpPosition = -1;
  currentReferenceSequence = "";
}

// Destructor.
statistics::~statistics(void) {}

void statistics::generateStatistics(vcf& v) {

// Initialise some variables.
  isTransition   = false;
  isTransversion = false;
  isInsertion    = false;
  isDeletion     = false;
  variantStruct variant;

// Check if this variant is annotated as being in dbsnp and/or hapmap.
  inDbsnp = (v.rsid == ".") ? false : true;
  inHapmap = (v.infoTags.count("HM3") == 0 && v.infoTags.count("HM3A") == 0) ? false : true;

// Determine if the event is an indel/SV or a SNP.  If the variant is
// an indel, keep track of whether it is an insertion or a deletion.  If
// the variant is a SNP, determine if it is a transition of a transversion.
  if (v.ref.size() - v.alt.size() == 0 && v.ref.size() == 1) {

// Generate a string as a pair the pair of alleles, in lower case and in alphabetical
// order.  A simple comparison can then be made to determine if the SNP is a 
// transition or a transversion.
    string alleles = v.ref + v.alt;
    for (int i = 0; i < 2; i++) {alleles[i] = tolower(alleles[i]);}
    sort(alleles.begin(), alleles.end());

    // Populate hapmap values.
    if (inHapmap) {variants[v.referenceSequence][v.filters].hapmap += 1;}

    // Transition:   A <-> G or C <-> T.
    if (alleles == "ag" || alleles == "ct") {
      isTransition = true;
      if (inDbsnp) {variants[v.referenceSequence][v.filters].knownTransitions += 1;}
      else {variants[v.referenceSequence][v.filters].novelTransitions += 1;}
    }

    // Transversion: A <-> C, A <-> T, C <-> G or G <-> T.
    if (alleles == "ac" || alleles == "at" || alleles == "cg" || alleles == "gt") {
      isTransversion = true;
      if (inDbsnp) {variants[v.referenceSequence][v.filters].knownTransversions += 1;}
      else {variants[v.referenceSequence][v.filters].novelTransversions += 1;}
    }

// Keep track of the last SNP position, so that the distribution of 
// the distance between SNPs can be maintained.
    if (v.referenceSequence != currentReferenceSequence) {
      currentReferenceSequence = v.referenceSequence;
      lastSnpPosition = -1;
    }
    if (lastSnpPosition != -1) {
      unsigned int distance = v.position - lastSnpPosition;
      snpDistribution[distance] += 1;
    }
  }
  else if (v.ref.size() - v.alt.size() == 0) {
    cout << "Haven't handled MNPs yet." << endl;
    exit(0);
  }
  // Deletions.
  else if (v.ref.size() > v.alt.size()) {
    isDeletion = true;
    variants[v.referenceSequence][v.filters].deletions[v.ref.size() - v.alt.size()] += 1;
  }
  // Insertions.
  else if (v.ref.size() < v.alt.size()) {
    isInsertion = true;
    variants[v.referenceSequence][v.filters].insertions[v.alt.size() - v.ref.size()] += 1;
  }
}

// Print out the statistics to the output file.
void statistics::printStatistics(ostream* output) {
  countByFilter();
  *output << setw(22) << "";
  *output << setw(60) << "--------------------------# SNPs--------------------------";
  *output << setw(9) << "";
  *output << setw(24) << "------ts/tv ratio-----";
  *output << endl;
  *output << setw(22) << "filter";
  *output << setw(12) << "total";
  *output << setw(12) << "novel ts";
  *output << setw(12) << "novel tv";
  *output << setw(12) << "known ts";
  *output << setw(12) << "known tv";
  *output << setw(9) << setprecision(6) << "% dbsnp";
  *output << setw(8) << setprecision(6) << "total";
  *output << setw(8) << setprecision(6) << "novel";
  *output << setw(8) << setprecision(6) << "known";
  *output << setw(12) << "hapmap";
  *output << endl;

// Print the total number of variants over all reference sequences.
  for (map<string, variantStruct>::iterator iter = totalVariants["total"].begin(); iter != totalVariants["total"].end(); iter++) {
    if ((*iter).first != "all" && (*iter).first != "PASS") {
      string filter = (*iter).first;
      printVariantStruct(output, filter, (*iter).second);
    }
  }
  *output << setw(22) << "";
  *output << "---------------------------------------------------------------------------------------------------------";
  *output << endl;
  string filter = "PASS";
  printVariantStruct(output, filter, totalVariants["total"]["PASS"]);
  filter = "Total";
  printVariantStruct(output, filter, totalVariants["total"]["all"]);
  *output << setw(22) << "";
  *output << "---------------------------------------------------------------------------------------------------------";
  *output << endl;
}

// The structure containing the numbers of the different variant types is
// currently organised by reference sequence and by filter.  Some of the
// filters, however, are combinations of multiple filters (e.g. Q10,DP100
// could represent a variant filtered out as it is below a quality threshold
// of 10 and a depth threshold of 100).  Calculate the total number of variants
// of each kind (e.g novel transition), for each filter and count the total
// over all reference sequences.
void statistics::countByFilter() {
  for (map<string, map<string, variantStruct> >::iterator iter = variants.begin(); iter != variants.end(); iter++) {
    for (map<string, variantStruct>::iterator filterIter = (*iter).second.begin(); filterIter != (*iter).second.end(); filterIter++) {
      vector<string> filters = split((*filterIter).first, ";");
      totalVariants[(*iter).first]["all"] = totalVariants[(*iter).first]["all"] + variants[(*iter).first][(*filterIter).first];
      totalVariants["total"]["all"] = totalVariants["total"]["all"] + variants[(*iter).first][(*filterIter).first];
      for (vector<string>::iterator vIter = filters.begin(); vIter != filters.end(); vIter++) {
        totalVariants[(*iter).first][(*vIter)] = totalVariants[(*iter).first][(*vIter)] + variants[(*iter).first][(*filterIter).first];
        totalVariants["total"][(*vIter)] = totalVariants["total"][(*vIter)] + variants[(*iter).first][(*filterIter).first];
      }
    }
  }
}

// Print the contents of the structure variantStruct to screen in a standard format.
void statistics::printVariantStruct(ostream* output, string& filter, variantStruct& var) {
  int novel         = var.novelTransitions + var.novelTransversions;
  int known         = var.knownTransitions + var.knownTransversions;
  int transitions   = var.novelTransitions + var.knownTransitions;
  int transversions = var.novelTransversions + var.knownTransversions;
  int totalSnp      = novel + known;

  float dbsnp     = (totalSnp == 0) ? 0. : (100. * float(known) / float(totalSnp));
  float tstv      = (transversions == 0) ? 0. : (float(transitions) / float(transversions));
  float noveltstv = (var.novelTransversions == 0) ? 0. : (float(var.novelTransitions) / float(var.novelTransversions));
  float knowntstv = (var.knownTransversions == 0) ? 0. : (float(var.knownTransitions) / float(var.knownTransversions));

  *output << setw(22) << filter;
  *output << setw(12) << setprecision(10) << totalSnp;
  *output << setw(12) << var.novelTransitions,
  *output << setw(12) << var.novelTransversions;
  *output << setw(12) << var.knownTransitions;
  *output << setw(12) << var.knownTransversions;
  *output << setw(9) << setprecision(4) << dbsnp;
  *output << setw(8) << setprecision(3) << tstv;
  *output << setw(8) << setprecision(3) << noveltstv;
  *output << setw(8) << setprecision(3) << knowntstv;
  *output << setw(12) << var.hapmap;
  *output << endl;
}
