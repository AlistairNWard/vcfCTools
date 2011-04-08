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
  hasIndel = false;
  hasMnp = false;
  hasSnp = false;
  hasAnnotations = false;
}

// Destructor.
statistics::~statistics(void) {}

void statistics::generateStatistics(vcf& v, bool generateAfs) {

// Initialise some variables.
  isTransition   = false;
  isTransversion = false;
  variantStruct variant;
  information info;
  string tag = "AC";

// Check if this variant is annotated as being in dbsnp.
  inDbsnp = (v.rsid == ".") ? false : true;

// Deal with multi-allelic variants.
  if (v.hasMultipleAlternates) {
    variants[v.referenceSequence][v.filters].multiAllelic++;

// Determine if the event is an indel/SV or a SNP.  If the variant is
// an indel, keep track of whether it is an insertion or a deletion.  If
// the variant is a SNP, determine if it is a transition of a transversion.
  } else {
    if (generateAfs) {info = v.getInfo(tag);}
    if (v.isSNP[0]) {
      hasSnp = true;

// Generate a string as a pair the pair of alleles, in lower case and in alphabetical
// order.  A simple comparison can then be made to determine if the SNP is a 
// transition or a transversion.
      string alleles = v.ref + v.alt[0];
      for (int i = 0; i < 2; i++) {alleles[i] = tolower(alleles[i]);}
      sort(alleles.begin(), alleles.end());

    // Transition:   A <-> G or C <-> T.
      if (alleles == "ag" || alleles == "ct") {
        isTransition = true;
        if (inDbsnp) {
          variants[v.referenceSequence][v.filters].knownTransitions += 1;
          if (v.infoTags.count("dbSNPX") != 0) {variants[v.referenceSequence][v.filters].diffKnownTransitions += 1;}
          if (generateAfs) {variants[v.referenceSequence][v.filters].afs[atoi(info.values[0].c_str())].knownTransitions++;}
        }
        else {
          variants[v.referenceSequence][v.filters].novelTransitions += 1;
          if (generateAfs) {variants[v.referenceSequence][v.filters].afs[atoi(info.values[0].c_str())].novelTransitions++;}
        }
      }

    // Transversion: A <-> C, A <-> T, C <-> G or G <-> T.
      if (alleles == "ac" || alleles == "at" || alleles == "cg" || alleles == "gt") {
        isTransversion = true;
        if (inDbsnp) {
          variants[v.referenceSequence][v.filters].knownTransversions += 1;
          if (v.infoTags.count("dbSNPX") != 0) {variants[v.referenceSequence][v.filters].diffKnownTransversions += 1;}
          if (generateAfs) {variants[v.referenceSequence][v.filters].afs[atoi(info.values[0].c_str())].knownTransversions++;}
        }
        else {
          variants[v.referenceSequence][v.filters].novelTransversions += 1;
          if (generateAfs) {variants[v.referenceSequence][v.filters].afs[atoi(info.values[0].c_str())].novelTransversions++;}
        }
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
    else if (v.isMNP[0]) {
      hasMnp = true;
      variants[v.referenceSequence][v.filters].mnps[v.alt[0].size()]++;
    // Deletions.
    } else if (v.isDeletion[0]) {
      hasIndel = true;
      variants[v.referenceSequence][v.filters].indels[v.ref.size() - v.alt[0].size()].deletions++;
    // Insertions.
    } else if (v.isInsertion[0]) {
      hasIndel = true;
      variants[v.referenceSequence][v.filters].indels[v.alt[0].size() - v.ref.size()].insertions++;
    }

// If the vcf file has been annotated and stats on these annotations are requested,
// search for the ANN tag and build stats for the following string.
    if (v.infoTags.count("ANN") != 0) {
      string tag = "ANN";
      information sInfo = v.getInfo(tag);
      hasAnnotations = true;

    // Update annotation statistics.
      for (vector<string>::iterator iter = sInfo.values.begin(); iter != sInfo.values.end(); iter++) {
        if (annotationNames.count((*iter)) == 0) {annotationNames[(*iter)] = 1;}
        if (v.isSNP[0] && isTransition) {variants[v.referenceSequence][v.filters].annotationsTs[(*iter)] += 1;}
        else if (v.isSNP[0] && isTransversion) {variants[v.referenceSequence][v.filters].annotationsTv[(*iter)] += 1;}
        else if (v.isDeletion[0]) {variants[v.referenceSequence][v.filters].annotationsDel[(*iter)] += 1;}
        else if (v.isInsertion[0]) {variants[v.referenceSequence][v.filters].annotationsIns[(*iter)] += 1;}
      }
    }
  }
}

// Print out the statistics to the output file.
void statistics::printSnpStatistics(ostream* output) {
  *output << setw(22) << "";
  *output << setw(60) << "--------------------------# SNPs--------------------------";
  *output << setw(18) << "";
  *output << setw(24) << "------ts/tv ratio-----";
  *output << endl;
  *output << setw(22) << "filter";
  *output << setw(12) << "total";
  *output << setw(12) << "novel ts";
  *output << setw(12) << "novel tv";
  *output << setw(12) << "known ts";
  *output << setw(12) << "known tv";
  *output << setw(18) << setprecision(6) << "% dbsnp (% diff)";
  *output << setw(8) << setprecision(6) << "total";
  *output << setw(8) << setprecision(6) << "novel";
  *output << setw(8) << setprecision(6) << "known";
  *output << setw(8) << "multi";
  *output << endl;

// Print the total number of variants over all reference sequences.
  for (map<string, variantStruct>::iterator iter = totalVariants["total"].begin(); iter != totalVariants["total"].end(); iter++) {
    if (iter->first != "all" && iter->first != "PASS") {
      string filter = iter->first;
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
  map<string, map<string, variantStruct> >::iterator variantIter;
  map<string, variantStruct>::iterator filterIter;
  vector<string>::iterator fIter;

  // Iterate over whole structure.
  for (variantIter = variants.begin(); variantIter != variants.end(); variantIter++) {

    // Iterate over the different filters.
    for (filterIter = variantIter->second.begin(); filterIter != variantIter->second.end(); filterIter++) {

      // Break up the filter string into all the individual filters.
      vector<string> filters = split(filterIter->first, ";");

      // Populate the structure containing information across all reference sequences
      // and filters.
      totalVariants[variantIter->first]["all"] = totalVariants[variantIter->first]["all"] + variants[variantIter->first][filterIter->first];
      totalVariants["total"]["all"] = totalVariants["total"]["all"] + variants[variantIter->first][filterIter->first];

      // Iterate over the individual filters.
      for (fIter = filters.begin(); fIter != filters.end(); fIter++) {
        totalVariants[variantIter->first][(*fIter)] = totalVariants[variantIter->first][(*fIter)] + variants[variantIter->first][filterIter->first];
        totalVariants["total"][(*fIter)] = totalVariants["total"][(*fIter)] + variants[variantIter->first][filterIter->first];
      }
    }
  }
}

// Print the contents of the structure variantStruct to screen in a standard format.
void statistics::printVariantStruct(ostream* output, string& filter, variantStruct& var) {
  int novel         = var.novelTransitions + var.novelTransversions;
  int known         = var.knownTransitions + var.knownTransversions;
  int diffKnown     = var.diffKnownTransitions + var.diffKnownTransversions;
  int transitions   = var.novelTransitions + var.knownTransitions;
  int transversions = var.novelTransversions + var.knownTransversions;
  int multiAllelic  = var.multiAllelic;
  int totalSnp      = novel + known + multiAllelic;

  float dbsnp     = (totalSnp == 0) ? 0. : (100. * float(known) / float(totalSnp));
  float dbsnpX    = (totalSnp == 0) ? 0. : (100. * float(diffKnown) / float(totalSnp));
  float tstv      = (transversions == 0) ? 0. : (float(transitions) / float(transversions));
  float noveltstv = (var.novelTransversions == 0) ? 0. : (float(var.novelTransitions) / float(var.novelTransversions));
  float knowntstv = (var.knownTransversions == 0) ? 0. : (float(var.knownTransitions) / float(var.knownTransversions));

  *output << setw(22) << filter;
  *output << setw(12) << setprecision(10) << totalSnp;
  *output << setw(12) << var.novelTransitions,
  *output << setw(12) << var.novelTransversions;
  *output << setw(12) << var.knownTransitions;
  *output << setw(12) << var.knownTransversions;
  *output << setw(9) << setprecision(3) << dbsnp;
  *output << setw(2) << " (";
  *output << setw(5) << setprecision(3) << dbsnpX;
  *output << setw(2) << ") ";
  *output << setw(8) << setprecision(3) << tstv;
  *output << setw(8) << setprecision(3) << noveltstv;
  *output << setw(8) << setprecision(3) << knowntstv;
  *output << setw(8) << multiAllelic;
  *output << endl;
}

// Print out annotation information for SNPs.
void statistics::printSnpAnnotations(ostream* output) {
  for (map<string, unsigned int>::iterator annIter = annotationNames.begin(); annIter != annotationNames.end(); annIter++) {
    *output << endl;
    string annotationName = annIter->first;
    *output << "SNP annotation information for: " << annotationName << endl;
    *output << endl;
    *output << setw(22) << "filter";
    *output << setw(16) << "total";
    *output << setw(16) << "transitions";
    *output << setw(16) << "transversions";
    *output << setw(16) << "ts/tv";
    *output << endl;
    for (map<string, variantStruct>::iterator iter = totalVariants["total"].begin(); iter != totalVariants["total"].end(); iter++) {
      if (iter->first != "all" && iter->first != "PASS") {
        string filter = iter->first;
        printSnpAnnotationStruct(output, filter, (*iter).second, annotationName);
      }
    }
    *output << setw(22) << "";
    *output << "--------------------------------------------------------------------";
    *output << endl;
    string filter = "PASS";
    printSnpAnnotationStruct(output, filter, totalVariants["total"]["PASS"], annotationName);
    filter = "Total";
    printSnpAnnotationStruct(output, filter, totalVariants["total"]["all"], annotationName);
    *output << setw(22) << "";
    *output << "--------------------------------------------------------------------";
    *output << endl;
  }
  *output << endl;
}

// Print out allele frequency information.
void statistics::printAfs(ostream* output) {
  map<unsigned int, snpTypes>::iterator iter;
  unsigned int transitions, transversions, novel, known;
  float dbsnp, tstv, novelTstv, knownTstv;

  *output << "Statistics by allele frequency:";
  *output << endl;
  *output << endl;
  *output << setw(16) << "";
  *output << setw(60) << "--------------------------# SNPs--------------------------";
  *output << setw(16) << "";
  *output << setw(24) << "-----ts/tv ratio-----";
  *output << endl;
  *output << setw(14) << "allele count";
  *output << setw(12) << "total";
  *output << setw(12) << "novel ts";
  *output << setw(12) << "novel tv";
  *output << setw(12) << "known ts";
  *output << setw(12) << "known tv";
  *output << setw(18) << "dbsnp";
  *output << setw(8) << "total";
  *output << setw(8) << "novel";
  *output << setw(8) << "known";
  *output << endl;
  for (iter = totalVariants["total"]["PASS"].afs.begin(); iter != totalVariants["total"]["PASS"].afs.end(); iter++) {
    novel = iter->second.novelTransitions + iter->second.novelTransversions;
    known = iter->second.knownTransitions + iter->second.knownTransversions;
    transitions = iter->second.novelTransitions + iter->second.knownTransitions;
    transversions = iter->second.novelTransversions + iter->second.knownTransversions;

    dbsnp = ( (known + novel) == 0 ) ? 0. : float(known) / ( float(novel) + float(known) );
    tstv = (transversions == 0) ? 0 : float(transitions) / float(transversions);
    novelTstv = (iter->second.novelTransversions == 0) ? 0 : float(iter->second.novelTransitions) / float(iter->second.novelTransversions);
    knownTstv = (iter->second.knownTransversions == 0) ? 0 : float(iter->second.knownTransitions) / float(iter->second.knownTransversions);
    *output << setw(12) << iter->first;
    *output << setw(12) << novel + known;
    *output << setw(12) << iter->second.novelTransitions;
    *output << setw(12) << iter->second.novelTransversions;
    *output << setw(12) << iter->second.knownTransitions;
    *output << setw(12) << iter->second.knownTransversions;
    *output << setw(18) << setprecision(3) << dbsnp;
    *output << setw(8) << setprecision(3) << tstv;
    *output << setw(8) << setprecision(3) << novelTstv;
    *output << setw(8) << setprecision(3) << knownTstv;
    *output << endl;
  }
}

// Print out the information structure for annotated SNPs.
void statistics::printSnpAnnotationStruct(ostream* output, string& filter, variantStruct& var, string& ann) {
  float tstv = (var.annotationsTv[ann] == 0) ? 0. : (float(var.annotationsTs[ann]) / float(var.annotationsTv[ann]));

  *output << setw(22) << filter;
  *output << setw(16) << setprecision(10) << var.annotationsTs[ann] + var.annotationsTv[ann];
  *output << setw(16) << var.annotationsTs[ann];
  *output << setw(16) << var.annotationsTv[ann];
  *output << setw(16) << setprecision(3) << tstv;
  *output << endl;
}

// Print out statistics on MNPs.
void statistics::printMnpStatistics(ostream* output) {
  string filterTag;

  *output << "Statistics on MNPs:" << endl;
  *output << endl;

  // First print out the total number of found MNPs (regardless of the filter).
  filterTag = "all";
  printMnpFilter(filterTag, output);

  // Now print out the total number of "PASS" MNPs.
  filterTag = "PASS";
  printMnpFilter(filterTag, output);
}

// Print out particular MNP statistics.
void statistics::printMnpFilter(string& tag, ostream* output) {
  map<unsigned int, unsigned int>::iterator iter;
  unsigned int totalMnps = 0;

  for (iter = totalVariants["total"][tag].mnps.begin(); iter != totalVariants["total"][tag].mnps.end(); iter++) {
    totalMnps += iter->second;
  }
  if (totalMnps != 0) {
    *output << "MNPs with filter field: " << tag << endl;
    *output << "Total number = " << totalMnps << endl;
    *output << endl;
    *output << left;
    *output << setw(15) << "length (bp)";
    *output << setw(10) << "number";
    *output << endl;
    for (iter = totalVariants["total"][tag].mnps.begin(); iter != totalVariants["total"][tag].mnps.end(); iter++) {
      *output << left;
      *output << setw(15) << iter->first;
      *output << setw(10) << iter->second;
      *output << endl;
    }
    *output << endl;
  }
}

// Print out statistics on indels.
void statistics::printIndelStatistics(ostream* output) {
  map<unsigned int, indel>::iterator iter;
  unsigned int totalInsertions = 0;
  unsigned int totalDeletions = 0;
  float ratio;

  *output << "Indel statistics:";
  *output << endl;
  *output << endl;
  *output << setw(12) << "length";
  *output << setw(12) << "insertions";
  *output << setw(12) << "deletions";
  *output << setw(12) << "ins/del";
  *output << endl;
  for (iter = totalVariants["total"]["PASS"].indels.begin(); iter != totalVariants["total"]["PASS"].indels.end(); iter++) {
    ratio = (iter->second.deletions == 0) ? 0 : float(iter->second.insertions) / float(iter->second.deletions);
    *output << setw(12) << iter->first;
    *output << setw(12) << iter->second.insertions;
    *output << setw(12) << iter->second.deletions;
    *output << setw(12) << setprecision(3) << ratio;
    *output << endl;
    totalInsertions += iter->second.insertions;
    totalDeletions += iter->second.deletions;
  }
  ratio = (totalDeletions == 0) ? 0 : float(totalInsertions) / float(totalDeletions);
  *output << endl;
  *output << "Total indels:     " << totalInsertions + totalDeletions << endl;
  *output << "Total insertions: " << totalInsertions << endl;
  *output << "Total deletions:  " << totalDeletions << endl;
  *output << "Ratio:            " << setprecision(3) << ratio;
  *output << endl;
}
