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
  lastSnpPosition = -1;
  lastIndelPosition = -1;
  currentReferenceSequence = "";
  hasIndel = false;
  hasMnp = false;
  hasSnp = false;

// Initialise the arrays.
  variants.clear();
}

// Destructor.
statistics::~statistics(void) {}

// Parse the variants at this locus and generate statistics.
void statistics::generateStatistics(variant& var, vcf& v, int position, bool useAnnotations, vector<string>& annFlags, bool generateAfs) {
  unsigned int ac;
  double af;
  variantInfo info;

// Deal with each alternate allele in turn.  Only generate statistics on the requested
// variant classes.  By default, this is all classes.

  // SNPs.
  if (var.processSnps) {
    unsigned int iterationNumber = 0;
    // Biallelic SNPs.
    for (var.variantIter = var.vmIter->second.biSnps.begin(); var.variantIter != var.vmIter->second.biSnps.end(); var.variantIter++) {
      iterationNumber++;
      info.processInfoFields(var.variantIter->info);

      if (generateAfs) {
        info.getInfo(string("AC"), var.variantIter->referenceSequence, var.vmIter->first);
        ac = atoi(info.values[0].c_str());
        info.getInfo(string("AF"), var.variantIter->referenceSequence, var.vmIter->first);
        af = atof(info.values[0].c_str());
      }
      // Check if this variant is annotated as being in dbsnp.
      inDbsnp = (var.variantIter->rsid == ".") ? false : true;
      hasSnp = true;

      // If this locus contains multiple SNPs, then the vcf file must contain at least two
      // SNPs at this locus.  In this case, treat the variant as a triallelic SNP.  If three
      // alternates come up, treat it as a quad-allelic.
      if (var.vmIter->second.biSnps.size() > 1) {
        if (iterationNumber == 1) {variants[var.variantIter->referenceSequence][var.variantIter->filters].multiAllelic++;}
      }
  
      // If this is the only biallelic SNP at this locus, determine if it is a transition or a transversion
      // and update the statistics accordingly.
      else {
        isTransition   = false;
        isTransversion = false;
  
        // Generate a string as a pair the pair of alleles, in lower case and in alphabetical
        // order.  A simple comparison can then be made to determine if the SNP is a 
        // transition or a transversion.
        string alleles = var.variantIter->ref + var.variantIter->altString;
        for (int i = 0; i < 2; i++) {alleles[i] = tolower(alleles[i]);}
        sort(alleles.begin(), alleles.end());
  
        // Transition:   A <-> G or C <-> T.
        if (alleles == "ag" || alleles == "ct") {
          isTransition = true;
          if (inDbsnp) {
            variants[var.variantIter->referenceSequence][var.variantIter->filters].knownTransitions++;
            if (info.infoTags.count("dbSNPX") != 0) {variants[var.variantIter->referenceSequence][var.variantIter->filters].diffKnownTransitions++;}
            if (generateAfs) {
              variants[var.variantIter->referenceSequence][var.variantIter->filters].acs[ac].knownTransitions++;
              variants[var.variantIter->referenceSequence][var.variantIter->filters].afs[af].knownTransitions++;
            }
          } else {
            variants[var.variantIter->referenceSequence][var.variantIter->filters].novelTransitions++;
            if (generateAfs) {
              variants[var.variantIter->referenceSequence][var.variantIter->filters].acs[ac].novelTransitions++;
              variants[var.variantIter->referenceSequence][var.variantIter->filters].afs[af].novelTransitions++;
            }
          }

          // Annotations.
          if (useAnnotations) {
            getAnnotations(annFlags, info, variants[var.variantIter->referenceSequence][var.variantIter->filters].annotationsTs);}
          
          //  for (vector<string>::iterator annIter = annFlags.begin(); annIter != annFlags.end(); annIter++) {
          //    cout << *annIter << endl;
          //  }

            //if (info.values.size() != 0) {
            //  hasAnnotations = true;
            //  for (vector<string>::iterator annIter = info.values.begin(); annIter != info.values.end(); annIter++) {
             //   if (annotationNames.count(*annIter) == 0) {annotationNames[*annIter] = 1;}
            //   variants[var.variantIter->referenceSequence][var.variantIter->filters].annotationsTs[*annIter]++;
            //  }
         
          //}
  
          // Transversion: A <-> C, A <-> T, C <-> G or G <-> T.
        } else if (alleles == "ac" || alleles == "at" || alleles == "cg" || alleles == "gt") {
          isTransversion = true;
          if (inDbsnp) {
            variants[var.variantIter->referenceSequence][var.variantIter->filters].knownTransversions++;
            if (info.infoTags.count("dbSNPX") != 0) {variants[var.variantIter->referenceSequence][var.variantIter->filters].diffKnownTransversions++;}
            if (generateAfs) {
              variants[var.variantIter->referenceSequence][var.variantIter->filters].acs[ac].knownTransversions++;
              variants[var.variantIter->referenceSequence][var.variantIter->filters].afs[af].knownTransversions++;
            }
          } else {
            variants[var.variantIter->referenceSequence][var.variantIter->filters].novelTransversions++;
            if (generateAfs) {
              variants[var.variantIter->referenceSequence][var.variantIter->filters].acs[ac].novelTransversions++; 
              variants[var.variantIter->referenceSequence][var.variantIter->filters].afs[af].novelTransversions++; 
            }
          }

          // Annotations.
          if (info.values.size() != 0) {
            hasAnnotations = true;
            for (vector<string>::iterator annIter = info.values.begin(); annIter != info.values.end(); annIter++) {
              if (annotationNames.count(*annIter) == 0) {annotationNames[*annIter] = 1;}
              variants[var.variantIter->referenceSequence][var.variantIter->filters].annotationsTv[*annIter]++;
            }
          }
        }
      }
  
      // Keep track of the last SNP position, so that the distribution of 
      // the distance between SNPs can be maintained.
      if (var.variantIter->referenceSequence != currentReferenceSequence) {
        currentReferenceSequence = var.variantIter->referenceSequence;
        lastSnpPosition = -1;
      }
      if (lastSnpPosition != -1) {
        unsigned int distance = position - lastSnpPosition;
        snpDistribution[distance]++;
      }
    }

    // Multiallelic SNPs.
    for (var.variantIter = var.vmIter->second.multiSnps.begin(); var.variantIter != var.vmIter->second.multiSnps.end(); var.variantIter++) {
      info.processInfoFields(var.variantIter->info);
      
      // Check for annotations.
      if (info.infoTags.count("ANN") != 0) {
        hasAnnotations = true;
        info.getInfo(string("ANN"), var.variantIter->referenceSequence, var.vmIter->first);
        if (info.values.size() != 0) {
          for (vector<string>::iterator annIter = info.values.begin(); annIter != info.values.end(); annIter++) {
            if (annotationNames.count(*annIter) == 0) {annotationNames[*annIter] = 1;}
            if (var.variantIter->isTriallelicSnp) {
              variants[var.variantIter->referenceSequence][var.variantIter->filters].annotationsTriallelicSnp[*annIter]++;
            } else if (var.variantIter->isQuadallelicSnp) {
              variants[var.variantIter->referenceSequence][var.variantIter->filters].annotationsQuadallelicSnp[*annIter]++;
            }
          }
        }
      }
      variants[var.variantIter->referenceSequence][var.variantIter->filters].multiAllelic++;
    }
  }

  // MNPs.
  if (var.processMnps) {
    for (var.variantIter = var.vmIter->second.mnps.begin(); var.variantIter != var.vmIter->second.mnps.end(); var.variantIter++) {
      info.processInfoFields(var.variantIter->info);

      // Check for annotations.
      if (info.infoTags.count("ANN") != 0) {
        hasAnnotations = true;
        info.getInfo(string("ANN"), var.variantIter->referenceSequence, var.vmIter->first);
        if (info.values.size() != 0) {
          for (vector<string>::iterator annIter = info.values.begin(); annIter != info.values.end(); annIter++) {
            if (annotationNames.count(*annIter) == 0) {annotationNames[*annIter] = 1;}
            variants[var.variantIter->referenceSequence][var.variantIter->filters].annotationsMnp[*annIter]++;
          }
        }
      }
      hasMnp = true;
      variants[var.variantIter->referenceSequence][var.variantIter->filters].mnps[var.variantIter->altString.size()]++;
    }
  }

  //Indels.
  if (var.processIndels) {
    for (var.variantIter = var.vmIter->second.indels.begin(); var.variantIter != var.vmIter->second.indels.end(); var.variantIter++) {
      info.processInfoFields(var.variantIter->info);

      // Check for annotations.
      if (info.infoTags.count("ANN") != 0) {
        hasAnnotations = true;
        info.getInfo(string("ANN"), var.variantIter->referenceSequence, var.vmIter->first);
        if (info.values.size() != 0) {
          for (vector<string>::iterator annIter = info.values.begin(); annIter != info.values.end(); annIter++) {
            if (annotationNames.count(*annIter) == 0) {annotationNames[*annIter] = 1;}
            if (var.variantIter->isInsertion) {variants[var.variantIter->referenceSequence][var.variantIter->filters].annotationsIns[*annIter]++;}
            if (var.variantIter->isDeletion) {variants[var.variantIter->referenceSequence][var.variantIter->filters].annotationsDel[*annIter]++;}
          }
        }
      }

      hasIndel = true;
      int indelSize = var.variantIter->ref.size() - var.variantIter->altString.size();
      if (var.variantIter->isDeletion) {variants[var.variantIter->referenceSequence][var.variantIter->filters].indels[indelSize].deletions++;}
      if (var.variantIter->isInsertion) {variants[var.variantIter->referenceSequence][var.variantIter->filters].indels[abs(indelSize)].insertions++;}
    }
  }
}

// Search for annotations in the info string.  This will either involve searching
// for flags from a given vector or searching for all flags in the info field.
void statistics::getAnnotations(vector<string>& annotationFlags, variantInfo& info, map<string, unsigned int>& annotation) {

  // Search for all flags in the info field.
  if (annotationFlags.size() == 1 && annotationFlags[0] == "all") {
    cout << "LOOKING FOR ALL" << endl;

  // Search for a specified set of flags.
  } else {
    for (vector<string>::iterator iter = annotationFlags.begin(); iter != annotationFlags.end(); iter++) {
      cout << *iter << endl;
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

  double dbsnp     = (totalSnp == 0) ? 0. : (100. * double(known) / double(totalSnp));
  double dbsnpX    = (totalSnp == 0) ? 0. : (100. * double(diffKnown) / double(totalSnp));
  double tstv      = (transversions == 0) ? 0. : (double(transitions) / double(transversions));
  double noveltstv = (var.novelTransversions == 0) ? 0. : (double(var.novelTransitions) / double(var.novelTransversions));
  double knowntstv = (var.knownTransversions == 0) ? 0. : (double(var.knownTransitions) / double(var.knownTransversions));

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
    //*output << endl;
    string annotationName = annIter->first;
    //*output << "SNP annotation information for: " << annotationName << endl;
    //*output << endl;
    //*output << setw(22) << "filter";
    //*output << setw(16) << "total";
    //*output << setw(16) << "transitions";
    //*output << setw(16) << "transversions";
    //*output << setw(16) << "ts/tv";
    //*output << endl;
    //for (map<string, variantStruct>::iterator iter = totalVariants["total"].begin(); iter != totalVariants["total"].end(); iter++) {
    //  if (iter->first != "all" && iter->first != "PASS") {
    //    string filter = iter->first;
    //    printSnpAnnotationStruct(output, filter, (*iter).second, annotationName);
    //  }
    //}
    //*output << setw(22) << "";
    //*output << "--------------------------------------------------------------------";
    //*output << endl;
    //string filter = "PASS";
    //printSnpAnnotationStruct(output, filter, totalVariants["total"]["PASS"], annotationName);
    //filter = "Total";
    //printSnpAnnotationStruct(output, filter, totalVariants["total"]["all"], annotationName);
    printSnpAnnotationStruct(output, annotationName, totalVariants["total"]["all"], annotationName);
    //*output << setw(22) << "";
    //*output << "--------------------------------------------------------------------";
    //*output << endl;
  }
  *output << endl;
}

// Print out the information structure for annotated SNPs.
void statistics::printSnpAnnotationStruct(ostream* output, string& filter, variantStruct& var, string& ann) {
  double tstv = (var.annotationsTv[ann] == 0) ? 0. : (double(var.annotationsTs[ann]) / double(var.annotationsTv[ann]));

  *output << setw(22) << filter;
  *output << setw(16) << setprecision(10) << var.annotationsTs[ann] + var.annotationsTv[ann];
  *output << setw(16) << var.annotationsTs[ann];
  *output << setw(16) << var.annotationsTv[ann];
  *output << setw(16) << setprecision(3) << tstv;
  *output << endl;
}

// Print out allele count information.
void statistics::printAcs(ostream* output) {
  map<unsigned int, snpTypes>::iterator acsIter;
  unsigned int transitions, transversions, novel, known;
  double dbsnp, tstv, novelTstv, knownTstv;

  *output << "Statistics by allele count:";
  *output << endl;
  *output << endl;
  *output << setw(18) << "";
  *output << setw(60) << "--------------------------# SNPs--------------------------";
  *output << setw(10) << "";
  *output << setw(24) << "-----ts/tv ratio-----";
  *output << endl;
  *output << setw(16) << "allele count";
  *output << setw(12) << "total";
  *output << setw(12) << "novel ts";
  *output << setw(12) << "novel tv";
  *output << setw(12) << "known ts";
  *output << setw(12) << "known tv";
  *output << setw(12) << "% dbsnp";
  *output << setw(8) << "total";
  *output << setw(8) << "novel";
  *output << setw(8) << "known";
  *output << endl;
  for (acsIter = totalVariants["total"]["PASS"].acs.begin(); acsIter != totalVariants["total"]["PASS"].acs.end(); acsIter++) {
    novel = acsIter->second.novelTransitions + acsIter->second.novelTransversions;
    known = acsIter->second.knownTransitions + acsIter->second.knownTransversions;
    transitions = acsIter->second.novelTransitions + acsIter->second.knownTransitions;
    transversions = acsIter->second.novelTransversions + acsIter->second.knownTransversions;

    dbsnp = ( (known + novel) == 0 ) ? 0. : 100. * (double(known) / ( double(novel) + double(known)) );
    tstv = (transversions == 0) ? 0 : double(transitions) / double(transversions);
    novelTstv = (acsIter->second.novelTransversions == 0) ? 0 : double(acsIter->second.novelTransitions) / double(acsIter->second.novelTransversions);
    knownTstv = (acsIter->second.knownTransversions == 0) ? 0 : double(acsIter->second.knownTransitions) / double(acsIter->second.knownTransversions);
    *output << setw(12) << acsIter->first;
    *output << setw(12) << novel + known;
    *output << setw(12) << acsIter->second.novelTransitions;
    *output << setw(12) << acsIter->second.novelTransversions;
    *output << setw(12) << acsIter->second.knownTransitions;
    *output << setw(12) << acsIter->second.knownTransversions;
    *output << setw(12) << setprecision(3) << dbsnp;
    *output << setw(8) << setprecision(3) << tstv;
    *output << setw(8) << setprecision(3) << novelTstv;
    *output << setw(8) << setprecision(3) << knownTstv;
    *output << endl;
  }
  *output << endl;
}

// Print out allele frequency information.
void statistics::printAfs(ostream* output) {
  map<double, snpTypes>::iterator afsIter;
  unsigned int transitions, transversions, novel, known;
  double dbsnp, tstv, novelTstv, knownTstv;

  *output << "Statistics by allele frequency:";
  *output << endl;
  *output << endl;
  *output << setw(18) << "";
  *output << setw(60) << "--------------------------# SNPs--------------------------";
  *output << setw(10) << "";
  *output << setw(24) << "-----ts/tv ratio-----";
  *output << endl;
  *output << setw(16) << "allele frequency";
  *output << setw(12) << "total";
  *output << setw(12) << "novel ts";
  *output << setw(12) << "novel tv";
  *output << setw(12) << "known ts";
  *output << setw(12) << "known tv";
  *output << setw(12) << "% dbsnp";
  *output << setw(8) << "total";
  *output << setw(8) << "novel";
  *output << setw(8) << "known";
  *output << endl;
  for (afsIter = totalVariants["total"]["PASS"].afs.begin(); afsIter != totalVariants["total"]["PASS"].afs.end(); afsIter++) {
    novel = afsIter->second.novelTransitions + afsIter->second.novelTransversions;
    known = afsIter->second.knownTransitions + afsIter->second.knownTransversions;
    transitions = afsIter->second.novelTransitions + afsIter->second.knownTransitions;
    transversions = afsIter->second.novelTransversions + afsIter->second.knownTransversions;

    dbsnp = ( (known + novel) == 0 ) ? 0. : 100. * (double(known) / ( double(novel) + double(known)) );
    tstv = (transversions == 0) ? 0 : double(transitions) / double(transversions);
    novelTstv = (afsIter->second.novelTransversions == 0) ? 0 : double(afsIter->second.novelTransitions) / double(afsIter->second.novelTransversions);
    knownTstv = (afsIter->second.knownTransversions == 0) ? 0 : double(afsIter->second.knownTransitions) / double(afsIter->second.knownTransversions);
    *output << setw(12) << afsIter->first;
    *output << setw(12) << novel + known;
    *output << setw(12) << afsIter->second.novelTransitions;
    *output << setw(12) << afsIter->second.novelTransversions;
    *output << setw(12) << afsIter->second.knownTransitions;
    *output << setw(12) << afsIter->second.knownTransversions;
    *output << setw(12) << setprecision(6) << dbsnp;
    *output << setw(8) << setprecision(6) << tstv;
    *output << setw(8) << setprecision(6) << novelTstv;
    *output << setw(8) << setprecision(6) << knownTstv;
    *output << endl;
  }
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
  double ratio;

  *output << "Indel statistics:";
  *output << endl;
  *output << endl;
  *output << setw(12) << "length";
  *output << setw(12) << "insertions";
  *output << setw(12) << "deletions";
  *output << setw(12) << "ins/del";
  *output << endl;
  for (iter = totalVariants["total"]["PASS"].indels.begin(); iter != totalVariants["total"]["PASS"].indels.end(); iter++) {
    ratio = (iter->second.deletions == 0) ? 0 : double(iter->second.insertions) / double(iter->second.deletions);
    *output << setw(12) << iter->first;
    *output << setw(12) << iter->second.insertions;
    *output << setw(12) << iter->second.deletions;
    *output << setw(12) << setprecision(3) << ratio;
    *output << endl;
    totalInsertions += iter->second.insertions;
    totalDeletions += iter->second.deletions;
  }
  ratio = (totalDeletions == 0) ? 0 : double(totalInsertions) / double(totalDeletions);
  *output << endl;
  *output << "Total indels:     " << totalInsertions + totalDeletions << endl;
  *output << "Total insertions: " << totalInsertions << endl;
  *output << "Total deletions:  " << totalDeletions << endl;
  *output << "Ratio:            " << setprecision(3) << ratio;
  *output << endl;
}
