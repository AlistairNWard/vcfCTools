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
  lastSnpPosition          = -1;
  lastSnpPosition          = -1;
  lastIndelPosition        = -1;
  currentReferenceSequence = "";
  hasInsertion             = false;
  hasDeletion              = false;
  hasMnp                   = false;
  hasSnp                   = false;
  splitMnps                = false;

  // Initialise the arrays.
  variants.clear();

  // Initialise sample level statistics.
  processSampleSnps = false;
  generateDetailed  = false;
  sampleSnps.clear();
}

// Destructor.
statistics::~statistics(void) {}

// Parse the variants at this locus and generate statistics.
//void statistics::generateStatistics(variant& var, vcf& v, int position, bool useAnnotations, vector<string>& annFlags, bool generateAfs, ostream* output) {
void statistics::generateStatistics(variant& var, bool useAnnotations, vector<string>& annFlags, bool generateAfs, output& ofile) {
  unsigned int ac;
  unsigned int variantID;
  string alleles;

  // Loop over all records at this locus.
  var.ovIter = var.ovmIter->second.begin();
  for (; var.ovIter != var.ovmIter->second.end(); var.ovIter++) {

    // Keep track of the different variants at this locus.  This is used to ensure
    // that multiallelic sites are correctly handled.
    bool locusHasSnp       = false;
    bool locusHasMultiSnp  = false;
    bool locusHasMnp       = false;
    bool locusHasInsertion = false;
    bool locusHasDeletion  = false;

    // Loop over all of the variants at this position and add to
    // the stats structure if required.
    vector<string>::iterator refIter       = var.ovIter->reducedRef.begin();
    vector<string>::iterator altIter       = var.ovIter->reducedAlts.begin();;
    vector<variantType>::iterator typeIter = var.ovIter->type.begin();;
    variantID = 0;
    for (; refIter != var.ovIter->reducedRef.end(); refIter++) {
   
      // Biallelic SNPs.
      if (typeIter->isBiallelicSnp) {

        // Check if this variant is annotated as being in dbsnp.
        inDbsnp     = (var.ovIter->rsid == ".") ? false : true;
        hasSnp      = true;
        locusHasSnp = true;

        // Determine if the SNP is a transition or transversion.
        isTransition   = false;
        isTransversion = false;
  
        // Generate a string as a pair the pair of alleles, in lower case and in alphabetical
        // order.  A simple comparison can then be made to determine if the SNP is a 
        // transition or a transversion.
        alleles = *refIter + *altIter;
        for (int i = 0; i < 2; i++) {alleles[i] = tolower(alleles[i]);}
        sort(alleles.begin(), alleles.end());
  
        // Retrieve information from the info fields if necessary.
        ac = 0;
        if (generateAfs || processSampleSnps) {
          variantInfo info(var.ovIter->info, var.headerInfoFields);
          info.retrieveFields();
          ac = atoi(info.infoFields["AC"].values[variantID].c_str());
          //info.getInfo(string("AF"), var.variantIter->referenceSequence, var.vmIter->first);
          //af = atof(info.values[0].c_str());
        }

        // Determine if the SNP is a transition or a transversion and update the relevant
        // statistics.
        determineSnpType(var, alleles, ac);

      // Triallelic SNPs.  If hasSnp is true, another SNP allele has already been observed
      // and counted at this locus, so do not double count.
      } else if (typeIter->isTriallelicSnp) {
        if (!locusHasSnp) {variants[var.ovIter->referenceSequence][var.ovIter->filters].multiAllelic++;}
        hasSnp           = true;
        hasMultiSnp      = true;
        locusHasSnp      = true;
        locusHasMultiSnp = true;
  
      // Quadallelic SNPs.
      } else if (typeIter->isQuadallelicSnp) {
        if (!locusHasSnp) {variants[var.ovIter->referenceSequence][var.ovIter->filters].multiAllelic++;}
        hasSnp           = true;
        hasMultiSnp      = true;
        locusHasSnp      = true;
        locusHasMultiSnp = true;
  
      // MNPs.
      } else if (typeIter->isMnp) {
        variants[var.ovIter->referenceSequence][var.ovIter->filters].mnps[altIter->size()]++;
        hasMnp      = true;
        locusHasMnp = true;
        if (splitMnps) {
          for (int i = 0; i < refIter->size(); i++) {
            if ((*refIter).substr(i, 1) == (*altIter).substr(i, 1)) {continue;}
            alleles = (*refIter).substr(i, 1) + (*altIter).substr(i, 1);
            for (int i = 0; i < 2; i++) {alleles[i] = tolower(alleles[i]);}
            sort(alleles.begin(), alleles.end());
            determineSnpType(var, alleles, ac);
          }
        }

      // Insertions.
      } else if (typeIter->isInsertion) {
        int insertionSize = altIter->size() - refIter->size();
        variants[var.ovIter->referenceSequence][var.ovIter->filters].indels[insertionSize].insertions++;
        hasInsertion      = true;
        locusHasInsertion = true;
  
      // Deletions.
      } else if (typeIter->isDeletion) {
        int deletionSize = refIter->size() - altIter->size();
        variants[var.ovIter->referenceSequence][var.ovIter->filters].indels[deletionSize].deletions++;
        hasDeletion      = true;
        locusHasDeletion = true;
  
      // Unknown variant type.
      } else {
        cerr << "ERROR: Unknown variant type." << endl;
        cerr << "Variant at " << var.ovIter->referenceSequence;
        cerr << ":" << var.ovmIter->first << endl;
        exit(1);
      }
  
      // Iterate the alt allele and allele types.
      altIter++;
      typeIter++;
      variantID++;
    }

//  unsigned int ac;
//  double af;
//  variantInfo info;

// Deal with each alternate allele in turn.  Only generate statistics on the requested
// variant classes.  By default, this is all classes.

  // SNPs.
//  if (var.processSnps) {
//    unsigned int iterationNumber = 0;
//    // Biallelic SNPs.
//    for (var.variantIter = var.vmIter->second.biSnps.begin(); var.variantIter != var.vmIter->second.biSnps.end(); var.variantIter++) {
//      iterationNumber++;
//      info.processInfoFields(var.variantIter->info);
//
//
//      // If this locus contains multiple SNPs, then the vcf file must contain at least two
//      // SNPs at this locus.  In this case, treat the variant as a triallelic SNP.  If three
//      // alternates come up, treat it as a quad-allelic.
//      if (var.vmIter->second.biSnps.size() > 1) {
//        if (iterationNumber == 1) {variants[var.variantIter->referenceSequence][var.variantIter->filters].multiAllelic++;}
//  
//      // If this is the only biallelic SNP at this locus, determine if it is a transition or a transversion
//      // and update the statistics accordingly.
//      } else {
//        isTransition   = false;
//        isTransversion = false;
//  
//        // Generate a string as a pair the pair of alleles, in lower case and in alphabetical
//        // order.  A simple comparison can then be made to determine if the SNP is a 
//        // transition or a transversion.
//        string alleles = var.variantIter->ref + var.variantIter->altString;
//        for (int i = 0; i < 2; i++) {alleles[i] = tolower(alleles[i]);}
//        sort(alleles.begin(), alleles.end());
//  
//        // Transition:   A <-> G or C <-> T.
//        if (alleles == "ag" || alleles == "ct") {
//          isTransition = true;
//
//          // Sample level and detailed stats.
//          if (processSampleSnps) {updateSampleSnps(var, v, ac);}
//          if (generateDetailed) {updateDetailedSnps(var, v, ac, output);}
//
//          if (inDbsnp) {
//            variants[var.variantIter->referenceSequence][var.variantIter->filters].knownTransitions++;
//            if (info.infoTags.count("dbSNPX") != 0) {variants[var.variantIter->referenceSequence][var.variantIter->filters].diffKnownTransitions++;}
//            if (generateAfs) {
//              variants[var.variantIter->referenceSequence][var.variantIter->filters].acs[ac].knownTransitions++;
//              variants[var.variantIter->referenceSequence][var.variantIter->filters].afs[af].knownTransitions++;
//            }
//          } else {
//            variants[var.variantIter->referenceSequence][var.variantIter->filters].novelTransitions++;
//            if (generateAfs) {
//              variants[var.variantIter->referenceSequence][var.variantIter->filters].acs[ac].novelTransitions++;
//              variants[var.variantIter->referenceSequence][var.variantIter->filters].afs[af].novelTransitions++;
//            }
//          }
//
//          // Annotations.
//          if (useAnnotations) {getAnnotations(annFlags, info, variants[var.variantIter->referenceSequence][var.variantIter->filters].annotationsTs);}
//          
//          // Transversion: A <-> C, A <-> T, C <-> G or G <-> T.
//        } else if (alleles == "ac" || alleles == "at" || alleles == "cg" || alleles == "gt") {
//          isTransversion = true;
//
//          // Sample level stats.
//          if (processSampleSnps) {updateSampleSnps(var, v, ac);}
//          if (generateDetailed) {updateDetailedSnps(var, v, ac, output);}
//
//          if (inDbsnp) {
//            variants[var.variantIter->referenceSequence][var.variantIter->filters].knownTransversions++;
//            if (info.infoTags.count("dbSNPX") != 0) {variants[var.variantIter->referenceSequence][var.variantIter->filters].diffKnownTransversions++;}
//            if (generateAfs) {
//              variants[var.variantIter->referenceSequence][var.variantIter->filters].acs[ac].knownTransversions++;
//              variants[var.variantIter->referenceSequence][var.variantIter->filters].afs[af].knownTransversions++;
//            }
//          } else {
//            variants[var.variantIter->referenceSequence][var.variantIter->filters].novelTransversions++;
//            if (generateAfs) {
//              variants[var.variantIter->referenceSequence][var.variantIter->filters].acs[ac].novelTransversions++; 
//              variants[var.variantIter->referenceSequence][var.variantIter->filters].afs[af].novelTransversions++; 
//            }
//          }
//
//          // Annotations.
//          if (useAnnotations) {getAnnotations(annFlags, info, variants[var.variantIter->referenceSequence][var.variantIter->filters].annotationsTv);}
//        }
//      }
//  
//      // Keep track of the last SNP position, so that the distribution of 
//      // the distance between SNPs can be maintained.
//      if (var.variantIter->referenceSequence != currentReferenceSequence) {
//        currentReferenceSequence = var.variantIter->referenceSequence;
//        lastSnpPosition = -1;
//      }
//      if (lastSnpPosition != -1) {
//        unsigned int distance = position - lastSnpPosition;
//        snpDistribution[distance]++;
//      }
//    }
//
//    // Multiallelic SNPs.
//    for (var.variantIter = var.vmIter->second.multiSnps.begin(); var.variantIter != var.vmIter->second.multiSnps.end(); var.variantIter++) {
//      info.processInfoFields(var.variantIter->info);
//      
//      // Check for annotations.
//      if (useAnnotations) {
//        if (var.variantIter->isTriallelicSnp) {
//          getAnnotations(annFlags, info, variants[var.variantIter->referenceSequence][var.variantIter->filters].annotationsTriallelicSnp);
//        } else if (var.variantIter->isQuadallelicSnp) {
//          getAnnotations(annFlags, info, variants[var.variantIter->referenceSequence][var.variantIter->filters].annotationsQuadallelicSnp);
//        }
//      }
//      variants[var.variantIter->referenceSequence][var.variantIter->filters].multiAllelic++;
//    }
//  }
//
//  // MNPs.
//  if (var.processMnps) {
//    for (var.variantIter = var.vmIter->second.mnps.begin(); var.variantIter != var.vmIter->second.mnps.end(); var.variantIter++) {
//      info.processInfoFields(var.variantIter->info);
//      hasMnp = true;
//
//      // Check for annotations.
//      if (useAnnotations) {getAnnotations(annFlags, info, variants[var.variantIter->referenceSequence][var.variantIter->filters].annotationsMnp);}
//      variants[var.variantIter->referenceSequence][var.variantIter->filters].mnps[var.variantIter->altString.size()]++;
//    }
//  }
//
//  //Indels.
//  if (var.processIndels) {
//    for (var.variantIter = var.vmIter->second.indels.begin(); var.variantIter != var.vmIter->second.indels.end(); var.variantIter++) {
//      info.processInfoFields(var.variantIter->info);
//
//      // Check for annotations.
//      if (useAnnotations) {
//        if (var.variantIter->isInsertion) {
//          getAnnotations(annFlags, info, variants[var.variantIter->referenceSequence][var.variantIter->filters].annotationsIns);
//        } else if (var.variantIter->isDeletion) {
//          getAnnotations(annFlags, info, variants[var.variantIter->referenceSequence][var.variantIter->filters].annotationsDel);
//        }
//      }
//
//      hasIndel = true;
//      int indelSize = var.variantIter->ref.size() - var.variantIter->altString.size();
//      if (var.variantIter->isDeletion) {variants[var.variantIter->referenceSequence][var.variantIter->filters].indels[indelSize].deletions++;}
//      if (var.variantIter->isInsertion) {variants[var.variantIter->referenceSequence][var.variantIter->filters].indels[abs(indelSize)].insertions++;}
//    }
//  }
  }
}

// Given the SNP alleles, determine if it is a transition/transversion and update
// all necessary statistics.
void statistics::determineSnpType(variant& var, string& alleles, unsigned int ac) {

  // Transition:   A <-> G or C <-> T.
  if (alleles == "ag" || alleles == "ct") {
    isTransition = true;

    // Update the structures for sample level statistics.
    if (processSampleSnps) {updateSampleSnps(var, ac);}

    if (inDbsnp) {
      variants[var.ovIter->referenceSequence][var.ovIter->filters].knownTransitions++;
      //if (info.infoTags.count("dbSNPX") != 0) {variants[var.ovIter->referenceSequence][var.ovIter->filters].diffKnownTransitions++;}
    } else {
      variants[var.ovIter->referenceSequence][var.ovIter->filters].novelTransitions++;
    }

  // Transversion: A <-> C, A <-> T, C <-> G or G <-> T.
  } else if (alleles == "ac" || alleles == "at" || alleles == "cg" || alleles == "gt") {
    isTransversion = true;

    // Update the structures for sample level statistics.
    if (processSampleSnps) {updateSampleSnps(var, ac);}

    if (inDbsnp) {
      variants[var.ovIter->referenceSequence][var.ovIter->filters].knownTransversions++;
      //if (info.infoTags.count("dbSNPX") != 0) {variants[var.ovIter->referenceSequence][var.ovIter->filters].diffKnownTransversions++;}
    } else {
      variants[var.ovIter->referenceSequence][var.ovIter->filters].novelTransversions++;
    }
  }
}

// Update the map containing sample level statistics for SNPs.
void statistics::updateSampleSnps(variant& var, unsigned int ac) {
  unsigned int depth;
  unsigned int sampleID = 0;
  string entry;
  string geno;
  double quality;
  vector<string> sampleEntries;
  genotypeInfo gen(var.ovIter->genotypeFormat, var.ovIter->genotypes, var.headerFormatFields);
  gen.processFormats();

  // Parse each sample in turn.
  vector<string>::iterator genoIter = gen.genotypes.begin();
  for (; genoIter != gen.genotypes.end(); genoIter++) {
    sampleEntries = split(*genoIter, ":");

    // Check that the number of entries is consistent with the format string.
    if (sampleEntries.size() != gen.genotypeFormats.size() && sampleEntries.size() != 1) {
      cerr << "ERROR: Number of fields in the genotype string is inconsistent with the format." << endl;
      cerr << "Reference sequence: " << var.ovIter->referenceSequence << ", position: " << var.ovIter->position << "." << endl;
      exit(1);
    } else if (sampleEntries.size() == 1 && sampleEntries[0] == ".") {
      sampleSnps[var.samples[sampleID]].unknown++;
    } else {
      quality = (gen.genotypeFields.count("GQ") == 0) ? 0. : atof( (sampleEntries[gen.genotypeFields["GQ"].ID]).c_str() );
      if (quality >= minGenotypeQuality || gen.genotypeFields.count("GQ") == 0) {
        geno  = sampleEntries[gen.genotypeFields["GT"].ID];
        depth = (gen.genotypeFields.count("GQ") == 0) ? 0 : atoi( (sampleEntries[gen.genotypeFields["DP"].ID]).c_str() );

        // Update the total depth for the sample.
        sampleSnps[var.samples[sampleID]].totalDepth += depth;

        // Homozygous alternate SNPs.
        if (geno == "1/1" || geno == "1|1") {
          sampleSnps[var.samples[sampleID]].homAlt++;
          sampleSnps[var.samples[sampleID]].totalAltDepth += depth;
          if (inDbsnp) {
            if (isTransition) {sampleSnps[var.samples[sampleID]].knownTransitions++;}
            if (isTransversion) {sampleSnps[var.samples[sampleID]].knownTransversions++;}
          } else {
            if (isTransition) {sampleSnps[var.samples[sampleID]].novelTransitions++;}
            if (isTransversion) {sampleSnps[var.samples[sampleID]].novelTransversions++;}
          }
  
          // If the SNP is a singleton, update the singletons stat.
          if (ac == 1) {sampleSnps[var.samples[sampleID]].singletons++;}

        // Heterozygous SNPs.
        } else if (geno == "0/1" || geno == "1/0" || geno == "0|1" || geno == "1|0") {
          sampleSnps[var.samples[sampleID]].het++;
          sampleSnps[var.samples[sampleID]].totalAltDepth += depth;
          if (inDbsnp) {
            if (isTransition) {sampleSnps[var.samples[sampleID]].knownTransitions++;}
            if (isTransversion) {sampleSnps[var.samples[sampleID]].knownTransversions++;}
          } else {
            if (isTransition) {sampleSnps[var.samples[sampleID]].novelTransitions++;}
            if (isTransversion) {sampleSnps[var.samples[sampleID]].novelTransversions++;}
          }
  
          // If the SNP is a singleton, update the singletons stat.
          if (ac == 1) {sampleSnps[var.samples[sampleID]].singletons++;}
        } else if (geno == "0/0" || geno == "0|0") {
          sampleSnps[var.samples[sampleID]].homRef++;
        } else {
          sampleSnps[var.samples[sampleID]].unknown++;
        }
      }
    }
    sampleID++;
  }
}

// Print the header for detailed statistics.
void statistics::printDetailedHeader(output& ofile) {
  *ofile.outputStream << "Detailed statistics for each variant position:" << endl;
  *ofile.outputStream << setw(24) << "";
  *ofile.outputStream << setw(48) << "    ----------Homozygous reference------------";
  *ofile.outputStream << setw(48) << "    ---------------Heterozygous---------------";
  *ofile.outputStream << setw(48) << "    --------Homozygous non-reference----------";
  *ofile.outputStream << setw(12) << " -No genotype-";
  *ofile.outputStream << endl;
  *ofile.outputStream << setw(12) << "Ref._seq.";
  *ofile.outputStream << setw(12) << "Position";
  *ofile.outputStream << setw(12) << "Number";
  *ofile.outputStream << setw(12) << "Mean_depth";
  *ofile.outputStream << setw(12) << "Min_depth";
  *ofile.outputStream << setw(12) << "Max_depth";
  *ofile.outputStream << setw(12) << "Number";
  *ofile.outputStream << setw(12) << "Mean_depth";
  *ofile.outputStream << setw(12) << "Min_depth";
  *ofile.outputStream << setw(12) << "Max_depth";
  *ofile.outputStream << setw(12) << "Number";
  *ofile.outputStream << setw(12) << "Mean_depth";
  *ofile.outputStream << setw(12) << "Min_depth";
  *ofile.outputStream << setw(12) << "Max_depth";
  *ofile.outputStream << setw(12) << "Number";
  *ofile.outputStream << setw(20) << "SNP_type";
  *ofile.outputStream << setw(20) << "Het_samples";
  *ofile.outputStream << setw(20) << "Hom_alt_samples";
  *ofile.outputStream << endl;
}

// Search for annotations in the info string.  This will either involve searching
// for flags from a given vector or searching for all flags in the info field.
void statistics::getAnnotations(vector<string>& annotationFlags, variantInfo& info, map<string, unsigned int>& annotation) {
//  hasAnnotations = true;
//
//  // Search for all flags in the info field.
//  if (annotationFlags.size() == 1 && annotationFlags[0] == "all") {
//    for (map<string, string>::iterator iter = info.infoTags.begin(); iter != info.infoTags.end(); iter++) {
//      if (iter->second == "flag") {
//        annotation[iter->first]++;
//        annotationNames[iter->first] = 1;
//      }
//    }
//
//  // Search for a specified set of flags.
//  } else {
//    for (vector<string>::iterator iter = annotationFlags.begin(); iter != annotationFlags.end(); iter++) {
//      annotationNames[*iter] = 1;
//      if (info.infoTags.count(*iter) != 0) {annotation[*iter]++;}
//    }
//  }
}


// Update the map containing detailed statistics for SNPs.
void statistics::updateDetailedSnps(variant& var, vcf& v, unsigned int ac, output& ofile) {
//  vector<string> genotypes = var.extractGenotypeField( string("GT") );
//  vector<string> genotypeQualities = var.extractGenotypeField( string("GQ") );
//  vector<string> genotypeDepth = var.extractGenotypeField( string("DP") );
//
//  vector<string>::iterator qIter = genotypeQualities.begin();
//  vector<string>::iterator dIter = genotypeDepth.begin();
//  vector<string>::iterator sIter = v.samples.begin();
//
//  string hetSamples    = "";
//  string homAltSamples = "";
//
//  unsigned int unknown     = 0;
//  unsigned int minHomAlt   = 0;
//  unsigned int maxHomAlt   = 0;
//  unsigned int homAlt      = 0;
//  unsigned int homAltDepth = 0;
//  unsigned int minHet      = 0;
//  unsigned int maxHet      = 0;
//  unsigned int het         = 0;
//  unsigned int hetDepth    = 0;
//  unsigned int minHomRef   = 0;
//  unsigned int maxHomRef   = 0;
//  unsigned int homRef      = 0;
//  unsigned int homRefDepth = 0;
//
//  for (vector<string>::iterator gIter = genotypes.begin(); gIter != genotypes.end(); gIter++) {
//    double genotypeQuality = atof( (*qIter).c_str() );
//    double depth = atof( (*dIter).c_str() );
//    if (genotypeQuality >= minDetailedGenotypeQuality) {
//
//      // Homozygous alternate SNPs.
//      if (*gIter == "1/1") {
//        homAlt++;
//        homAltDepth += depth;
//        if (depth < minHomAlt || minHomAlt == 0) {minHomAlt = depth;}
//        if (depth > maxHomAlt) {maxHomAlt = depth;}
//        homAltSamples = (homAltSamples == "") ? *sIter : homAltSamples + "," + *sIter;
//
//      // Heterozygous SNPs.
//      } else if (*gIter == "0/1" || *gIter == "1/0") {
//        het++;
//        hetDepth += depth;
//        if (depth < minHet || minHet == 0) {minHet = depth;}
//        if (depth > maxHet) {maxHet = depth;}
//        hetSamples = (hetSamples == "") ? *sIter : hetSamples + "," + *sIter;
//
//      // Homozygous reference.
//      } else if (*gIter == "0/0") {
//        homRef++;
//        homRefDepth += depth;
//        if (depth < minHomRef || minHomRef == 0) {minHomRef = depth;}
//        if (depth > maxHomRef) {maxHomRef = depth;}
//
//      // Uncalled genotypes.
//      } else {
//        unknown++;
//      }
//    }
//    qIter++; // Iterate the genotype quality.
//    dIter++; // Iterate the genotype depth.
//    sIter++; // Iterate the sample name.
//  }
//
//  double aveHetDepth    = (het == 0) ? 0. : double(hetDepth) / double(het);
//  double aveHomRefDepth = (homRef == 0) ? 0. : double(homRefDepth) / double(homRef);
//  double aveHomAltDepth = (homAlt == 0) ? 0. : double(homAltDepth) / double(homAlt);
//  *ofile.outputStream << setw(12) << var.variantIter->referenceSequence;
//  *ofile.outputStream << setw(12) << var.vmIter->first;
//  *ofile.outputStream << setw(12) << homRef;
//  *ofile.outputStream << setw(12) << aveHomRefDepth;
//  *ofile.outputStream << setw(12) << minHomRef;
//  *ofile.outputStream << setw(12) << maxHomRef;
//  *ofile.outputStream << setw(12) << het;
//  *ofile.outputStream << setw(12) << aveHetDepth;
//  *ofile.outputStream << setw(12) << minHet;
//  *ofile.outputStream << setw(12) << maxHet;
//  *ofile.outputStream << setw(12) << homAlt;
//  *ofile.outputStream << setw(12) << aveHomAltDepth;
//  *ofile.outputStream << setw(12) << minHomAlt;
//  *ofile.outputStream << setw(12) << maxHomAlt;
//  *ofile.outputStream << setw(12) << unknown;
//  if (isTransition) {*ofile.outputStream << setw(20) << "transition";}
//  else if (isTransversion) {*ofile.outputStream << setw(20) << "transversion";}
//  else {*ofile.outputStream << setw(20) << "other";}
//  
//  // Output the lists of het and hom alt samples.  If there are none, set
//  // the string to "no-hets" or "no-hom-alts".
//  if (hetSamples == "") {hetSamples = "no-hets";}
//  if (homAltSamples == "") {homAltSamples = "no-hom-alts";}
//  *ofile.outputStream << "  " << hetSamples << "  ";
//  *ofile.outputStream << "  " << homAltSamples << "  ";
//  *ofile.outputStream << endl;
}

// The structure containing the numbers of the different variant types is
// currently organised by reference sequence and by filter.  Some of the
// filters, however, are combinations of multiple filters (e.g. Q10,DP100
// could represent a variant filtered out as it is below a quality threshold
// of 10 and a depth threshold of 100).  Calculate the total number of variants
// of each kind (e.g novel transition), for each filter and count the total
// over all reference sequences.
void statistics::countByFilter() {
  map<string, map<string, variantStruct> >::iterator variantIter = variants.begin();
  map<string, variantStruct>::iterator filterIter;
  vector<string>::iterator fIter;

  // Iterate over whole structure.
  for (; variantIter != variants.end(); variantIter++) {

    // Iterate over the different filters.
    for (filterIter = variantIter->second.begin(); filterIter != variantIter->second.end(); filterIter++) {

      // Break up the filter string into all the individual filters.
      vector<string> filters = split(filterIter->first, ";");

      // Populate the structure containing information across all reference sequences
      // and filters.
      totalVariants[variantIter->first]["all"] = totalVariants[variantIter->first]["all"] + filterIter->second;
      totalVariants["total"]["all"] = totalVariants["total"]["all"] + filterIter->second;

      // Iterate over the individual filters.
      for (fIter = filters.begin(); fIter != filters.end(); fIter++) {
        totalVariants[variantIter->first][(*fIter)] = totalVariants[variantIter->first][(*fIter)] + filterIter->second;
        totalVariants["total"][(*fIter)] = totalVariants["total"][(*fIter)] + filterIter->second;
      }
    }
  }
}

// Print out the statistics to the output file.
void statistics::printSnpStatistics(output& ofile) {

// Print the total number of variants over all reference sequences.

  bool writtenHeader = false;
  for (map<string, map< string, variantStruct> >::iterator iter = totalVariants.begin(); iter != totalVariants.end(); iter++) {
    for (map<string, variantStruct>::iterator vIter = iter->second.begin(); vIter != iter->second.end(); vIter++) {
      if (vIter->first == "PASS") {
        if (!writtenHeader) {
          *ofile.outputStream << "Statistics on SNPs that pass filters (marked as PASS)." << endl;
          *ofile.outputStream << endl;
          printHeader(ofile, string("reference_sequence"), true, true);
          *ofile.outputStream << endl;
          writtenHeader = true;
        }
        printVariantStruct(ofile, iter->first, vIter->second);
      }
    }
  }
 
  *ofile.outputStream << endl;
  *ofile.outputStream << "Total_statistics." << endl;
  *ofile.outputStream << endl;
  printHeader(ofile, string("filter"), true, true);
  *ofile.outputStream << endl;
  for (map<string, variantStruct>::iterator iter = totalVariants["total"].begin(); iter != totalVariants["total"].end(); iter++) {
    if (iter->first != "all" && iter->first != "PASS") {printVariantStruct(ofile, string(iter->first), iter->second);}
  }
  *ofile.outputStream << setw(22) << "";
  *ofile.outputStream << "---------------------------------------------------------------------------------------------------------";
  *ofile.outputStream << endl;
  string filter = "PASS";
  printVariantStruct(ofile, filter, totalVariants["total"]["PASS"]);
  filter = "Total";
  printVariantStruct(ofile, filter, totalVariants["total"]["all"]);
  *ofile.outputStream << setw(22) << "";
  *ofile.outputStream << "---------------------------------------------------------------------------------------------------------";
  *ofile.outputStream << endl;
}

// Print out a header line.
void statistics::printHeader(output& ofile, string text, bool dbsnpDiff, bool multi) {
  *ofile.outputStream << setw(22) << "";
  *ofile.outputStream << setw(60) << "--------------------------#SNPs---------------------------";
  if (dbsnpDiff) {*ofile.outputStream << setw(18) << "";}
  else {*ofile.outputStream << setw(12) << "";}
  *ofile.outputStream << setw(24) << "------ts/tv_ratio-----";
  *ofile.outputStream << endl;
  *ofile.outputStream << setw(22) << text;
  *ofile.outputStream << setw(12) << "total";
  *ofile.outputStream << setw(12) << "novel_ts";
  *ofile.outputStream << setw(12) << "novel_tv";
  *ofile.outputStream << setw(12) << "known_ts";
  *ofile.outputStream << setw(12) << "known_tv";
  if (dbsnpDiff) {*ofile.outputStream << setw(18) << setprecision(6) << "%dbsnp (%diff)";}
  else {*ofile.outputStream << setw(12) << "%dbsnp";}
  *ofile.outputStream << setw(8) << setprecision(6) << "total";
  *ofile.outputStream << setw(8) << setprecision(6) << "novel";
  *ofile.outputStream << setw(8) << setprecision(6) << "known";
  if (multi) {*ofile.outputStream << setw(8) << "multi";}
}

// Print the contents of the structure variantStruct to screen in a standard format.
void statistics::printVariantStruct(output& ofile, string filter, variantStruct& var) {
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

  *ofile.outputStream << setw(22) << filter;
  *ofile.outputStream << setw(12) << setprecision(10) << totalSnp;
  *ofile.outputStream << setw(12) << var.novelTransitions,
  *ofile.outputStream << setw(12) << var.novelTransversions;
  *ofile.outputStream << setw(12) << var.knownTransitions;
  *ofile.outputStream << setw(12) << var.knownTransversions;
  *ofile.outputStream << setw(9) << setprecision(3) << dbsnp;
  *ofile.outputStream << setw(2) << " (";
  *ofile.outputStream << setw(5) << setprecision(3) << dbsnpX;
  *ofile.outputStream << setw(2) << ") ";
  *ofile.outputStream << setw(8) << setprecision(3) << tstv;
  *ofile.outputStream << setw(8) << setprecision(3) << noveltstv;
  *ofile.outputStream << setw(8) << setprecision(3) << knowntstv;
  *ofile.outputStream << setw(8) << multiAllelic;
  *ofile.outputStream << endl;
}

// Print out annotation information for SNPs.
void statistics::printSnpAnnotations(output& ofile) {
  for (map<string, unsigned int>::iterator annIter = annotationNames.begin(); annIter != annotationNames.end(); annIter++) {
    //*ofile.outputStream << endl;
    string annotationName = annIter->first;
    //*ofile.outputStream << "SNP annotation information for: " << annotationName << endl;
    //*ofile.outputStream << endl;
    //*ofile.outputStream << setw(22) << "filter";
    //*ofile.outputStream << setw(16) << "total";
    //*ofile.outputStream << setw(16) << "transitions";
    //*ofile.outputStream << setw(16) << "transversions";
    //*ofile.outputStream << setw(16) << "ts/tv";
    //*ofile.outputStream << endl;
    //for (map<string, variantStruct>::iterator iter = totalVariants["total"].begin(); iter != totalVariants["total"].end(); iter++) {
    //  if (iter->first != "all" && iter->first != "PASS") {
    //    string filter = iter->first;
    //    printSnpAnnotationStruct(output, filter, (*iter).second, annotationName);
    //  }
    //}
    //*ofile.outputStream << setw(22) << "";
    //*ofile.outputStream << "--------------------------------------------------------------------";
    //*ofile.outputStream << endl;
    //string filter = "PASS";
    //printSnpAnnotationStruct(output, filter, totalVariants["total"]["PASS"], annotationName);
    //filter = "Total";
    //printSnpAnnotationStruct(output, filter, totalVariants["total"]["all"], annotationName);
    printSnpAnnotationStruct(ofile, annotationName, totalVariants["total"]["all"], annotationName);
    //*ofile.outputStream << setw(22) << "";
    //*ofile.outputStream << "--------------------------------------------------------------------";
    //*ofile.outputStream << endl;
  }
  *ofile.outputStream << endl;
}

// Print out the information structure for annotated SNPs.
void statistics::printSnpAnnotationStruct(output& ofile, string& filter, variantStruct& var, string& ann) {
  double tstv = (var.annotationsTv[ann] == 0) ? 0. : (double(var.annotationsTs[ann]) / double(var.annotationsTv[ann]));

  *ofile.outputStream << setw(22) << filter;
  *ofile.outputStream << setw(16) << setprecision(10) << var.annotationsTs[ann] + var.annotationsTv[ann];
  *ofile.outputStream << setw(16) << var.annotationsTs[ann];
  *ofile.outputStream << setw(16) << var.annotationsTv[ann];
  *ofile.outputStream << setw(16) << setprecision(3) << tstv;
  *ofile.outputStream << endl;
}

// Print out allele count information.
void statistics::printAcs(output& ofile) {
  map<unsigned int, snpTypes>::iterator acsIter;
  unsigned int transitions, transversions, novel, known;
  double dbsnp, tstv, novelTstv, knownTstv;

  *ofile.outputStream << "Statistics_by_allele_count:";
  *ofile.outputStream << endl;
  *ofile.outputStream << endl;
  *ofile.outputStream << setw(18) << "";
  *ofile.outputStream << setw(60) << "--------------------------#SNPs---------------------------";
  *ofile.outputStream << setw(10) << "";
  *ofile.outputStream << setw(24) << "-----ts/tv_ratio-----";
  *ofile.outputStream << endl;
  *ofile.outputStream << setw(16) << "allele_count";
  *ofile.outputStream << setw(12) << "total";
  *ofile.outputStream << setw(12) << "novel_ts";
  *ofile.outputStream << setw(12) << "novel_tv";
  *ofile.outputStream << setw(12) << "known_ts";
  *ofile.outputStream << setw(12) << "known_tv";
  *ofile.outputStream << setw(12) << "%dbsnp";
  *ofile.outputStream << setw(8) << "total";
  *ofile.outputStream << setw(8) << "novel";
  *ofile.outputStream << setw(8) << "known";
  *ofile.outputStream << endl;
  for (acsIter = totalVariants["total"]["PASS"].acs.begin(); acsIter != totalVariants["total"]["PASS"].acs.end(); acsIter++) {
    novel = acsIter->second.novelTransitions + acsIter->second.novelTransversions;
    known = acsIter->second.knownTransitions + acsIter->second.knownTransversions;
    transitions = acsIter->second.novelTransitions + acsIter->second.knownTransitions;
    transversions = acsIter->second.novelTransversions + acsIter->second.knownTransversions;

    dbsnp = ( (known + novel) == 0 ) ? 0. : 100. * (double(known) / ( double(novel) + double(known)) );
    tstv = (transversions == 0) ? 0 : double(transitions) / double(transversions);
    novelTstv = (acsIter->second.novelTransversions == 0) ? 0 : double(acsIter->second.novelTransitions) / double(acsIter->second.novelTransversions);
    knownTstv = (acsIter->second.knownTransversions == 0) ? 0 : double(acsIter->second.knownTransitions) / double(acsIter->second.knownTransversions);
    *ofile.outputStream << setw(12) << acsIter->first;
    *ofile.outputStream << setw(12) << novel + known;
    *ofile.outputStream << setw(12) << acsIter->second.novelTransitions;
    *ofile.outputStream << setw(12) << acsIter->second.novelTransversions;
    *ofile.outputStream << setw(12) << acsIter->second.knownTransitions;
    *ofile.outputStream << setw(12) << acsIter->second.knownTransversions;
    *ofile.outputStream << setw(12) << setprecision(3) << dbsnp;
    *ofile.outputStream << setw(8) << setprecision(3) << tstv;
    *ofile.outputStream << setw(8) << setprecision(3) << novelTstv;
    *ofile.outputStream << setw(8) << setprecision(3) << knownTstv;
    *ofile.outputStream << endl;
  }
  *ofile.outputStream << endl;
}

// Print out allele frequency information.
void statistics::printAfs(output& ofile) {
  map<double, snpTypes>::iterator afsIter;
  unsigned int transitions, transversions, novel, known;
  double dbsnp, tstv, novelTstv, knownTstv;

  *ofile.outputStream << "Statistics by allele frequency:";
  *ofile.outputStream << endl;
  *ofile.outputStream << endl;
  *ofile.outputStream << setw(18) << "";
  *ofile.outputStream << setw(60) << "--------------------------#SNPs---------------------------";
  *ofile.outputStream << setw(10) << "";
  *ofile.outputStream << setw(24) << "-----ts/tv_ratio-----";
  *ofile.outputStream << endl;
  *ofile.outputStream << setw(16) << "allele_frequency";
  *ofile.outputStream << setw(12) << "total";
  *ofile.outputStream << setw(12) << "novel_ts";
  *ofile.outputStream << setw(12) << "novel_tv";
  *ofile.outputStream << setw(12) << "known_ts";
  *ofile.outputStream << setw(12) << "known_tv";
  *ofile.outputStream << setw(12) << "%dbsnp";
  *ofile.outputStream << setw(8) << "total";
  *ofile.outputStream << setw(8) << "novel";
  *ofile.outputStream << setw(8) << "known";
  *ofile.outputStream << endl;
  for (afsIter = totalVariants["total"]["PASS"].afs.begin(); afsIter != totalVariants["total"]["PASS"].afs.end(); afsIter++) {
    novel = afsIter->second.novelTransitions + afsIter->second.novelTransversions;
    known = afsIter->second.knownTransitions + afsIter->second.knownTransversions;
    transitions = afsIter->second.novelTransitions + afsIter->second.knownTransitions;
    transversions = afsIter->second.novelTransversions + afsIter->second.knownTransversions;

    dbsnp = ( (known + novel) == 0 ) ? 0. : 100. * (double(known) / ( double(novel) + double(known)) );
    tstv = (transversions == 0) ? 0 : double(transitions) / double(transversions);
    novelTstv = (afsIter->second.novelTransversions == 0) ? 0 : double(afsIter->second.novelTransitions) / double(afsIter->second.novelTransversions);
    knownTstv = (afsIter->second.knownTransversions == 0) ? 0 : double(afsIter->second.knownTransitions) / double(afsIter->second.knownTransversions);
    *ofile.outputStream << setw(12) << afsIter->first;
    *ofile.outputStream << setw(12) << novel + known;
    *ofile.outputStream << setw(12) << afsIter->second.novelTransitions;
    *ofile.outputStream << setw(12) << afsIter->second.novelTransversions;
    *ofile.outputStream << setw(12) << afsIter->second.knownTransitions;
    *ofile.outputStream << setw(12) << afsIter->second.knownTransversions;
    *ofile.outputStream << setw(12) << setprecision(3) << dbsnp;
    *ofile.outputStream << setw(8) << setprecision(3) << tstv;
    *ofile.outputStream << setw(8) << setprecision(3) << novelTstv;
    *ofile.outputStream << setw(8) << setprecision(3) << knownTstv;
    *ofile.outputStream << endl;
  }
  *ofile.outputStream << endl;
}

// Print out statistics on MNPs.
void statistics::printMnpStatistics(output& ofile) {
  string filterTag;

  *ofile.outputStream << "Statistics_on_MNPs:" << endl;
  *ofile.outputStream << endl;

  // First print out the total number of found MNPs (regardless of the filter).
  //filterTag = "all";
  //printMnpFilter(filterTag, ofile);

  // Now print out the total number of "PASS" MNPs.
  filterTag = "PASS";
  printMnpFilter(filterTag, ofile);
}

// Print out particular MNP statistics.
void statistics::printMnpFilter(string& tag, output& ofile) {
  map<unsigned int, unsigned int>::iterator iter;
  unsigned int totalMnps = 0;

  for (iter = totalVariants["total"][tag].mnps.begin(); iter != totalVariants["total"][tag].mnps.end(); iter++) {
    totalMnps += iter->second;
  }
  if (totalMnps != 0) {
    *ofile.outputStream << "MNPs_with_filter_field: " << tag << endl;
    *ofile.outputStream << "Total_number=" << totalMnps << endl;
    *ofile.outputStream << endl;
    *ofile.outputStream << left;
    *ofile.outputStream << setw(15) << "length_(bp)";
    *ofile.outputStream << setw(10) << "number";
    *ofile.outputStream << endl;
    for (iter = totalVariants["total"][tag].mnps.begin(); iter != totalVariants["total"][tag].mnps.end(); iter++) {
      *ofile.outputStream << left;
      *ofile.outputStream << setw(15) << iter->first;
      *ofile.outputStream << setw(10) << iter->second;
      *ofile.outputStream << endl;
    }
    *ofile.outputStream << endl;
  }
}

// Print out statistics on indels.
void statistics::printIndelStatistics(output& ofile) {
  map<unsigned int, indel>::iterator iter;
  unsigned int totalInsertions = 0;
  unsigned int totalDeletions = 0;
  double ratio;

  *ofile.outputStream << "Indel_statistics:";
  *ofile.outputStream << endl;
  *ofile.outputStream << endl;
  *ofile.outputStream << setw(12) << "length";
  *ofile.outputStream << setw(12) << "insertions";
  *ofile.outputStream << setw(12) << "deletions";
  *ofile.outputStream << setw(12) << "ins/del";
  *ofile.outputStream << endl;
  for (iter = totalVariants["total"]["PASS"].indels.begin(); iter != totalVariants["total"]["PASS"].indels.end(); iter++) {
    ratio = (iter->second.deletions == 0) ? 0 : double(iter->second.insertions) / double(iter->second.deletions);
    *ofile.outputStream << setw(12) << iter->first;
    *ofile.outputStream << setw(12) << iter->second.insertions;
    *ofile.outputStream << setw(12) << iter->second.deletions;
    *ofile.outputStream << setw(12) << setprecision(3) << ratio;
    *ofile.outputStream << endl;
    totalInsertions += iter->second.insertions;
    totalDeletions += iter->second.deletions;
  }
  ratio = (totalDeletions == 0) ? 0 : double(totalInsertions) / double(totalDeletions);
  *ofile.outputStream << endl;
  *ofile.outputStream << "Total_indels:     " << totalInsertions + totalDeletions << endl;
  *ofile.outputStream << "Total_insertions: " << totalInsertions << endl;
  *ofile.outputStream << "Total_deletions:  " << totalDeletions << endl;
  *ofile.outputStream << "Ratio:            " << setprecision(3) << ratio;
  *ofile.outputStream << endl;
}

// Print out the sample level SNP statistics.
void statistics::printSampleSnps(vcf& v, output& ofile) {
  unsigned int totalTransitions = totalVariants["total"]["all"].novelTransitions + totalVariants["total"]["all"].knownTransitions;
  unsigned int totalTransversions = totalVariants["total"]["all"].novelTransversions + totalVariants["total"]["all"].knownTransversions;
  unsigned int totalSnps = totalTransitions + totalTransversions;

  *ofile.outputStream << "SNP_statistics_by_sample:" << endl;
  *ofile.outputStream << endl;
  printHeader(ofile, string("sample"), false, false);
  *ofile.outputStream << setw(12) << "Hom_ref";
  *ofile.outputStream << setw(12) << "Hets";
  *ofile.outputStream << setw(12) << "Hom_Alt";
  *ofile.outputStream << setw(12) << "Unknown";
  *ofile.outputStream << setw(12) << "Singletons";
  *ofile.outputStream << setw(12) << "Depth";
  *ofile.outputStream << setw(12) << "Alt_depth";
  *ofile.outputStream << endl;
  for (vector<string>::iterator sample = v.samples.begin(); sample != v.samples.end(); sample++) {
    int novel         = sampleSnps[*sample].novelTransitions + sampleSnps[*sample].novelTransversions;
    int known         = sampleSnps[*sample].knownTransitions + sampleSnps[*sample].knownTransversions;
    int transitions   = sampleSnps[*sample].novelTransitions + sampleSnps[*sample].knownTransitions;
    int transversions = sampleSnps[*sample].novelTransversions + sampleSnps[*sample].knownTransversions;
    int totalSnp      = novel + known;
    double depth      = (totalSnp == 0) ? 0. : double(sampleSnps[*sample].totalDepth) / double(totalSnps);
    double altDepth   = (totalSnp == 0) ? 0. : double(sampleSnps[*sample].totalAltDepth) / double(totalSnp);
    double dbsnp      = (totalSnp == 0) ? 0. : (100. * double(known) / double(totalSnp));
    double tstv       = (transversions == 0) ? 0. : (double(transitions) / double(transversions));
    double noveltstv  = (sampleSnps[*sample].novelTransversions == 0) ? 0. : (double(sampleSnps[*sample].novelTransitions) / double(sampleSnps[*sample].novelTransversions));
    double knowntstv  = (sampleSnps[*sample].knownTransversions == 0) ? 0. : (double(sampleSnps[*sample].knownTransitions) / double(sampleSnps[*sample].knownTransversions));    
    *ofile.outputStream << setw(22) << *sample;
    *ofile.outputStream << setw(12) << setprecision(10) << totalSnp;
    *ofile.outputStream << setw(12) << sampleSnps[*sample].novelTransitions,
    *ofile.outputStream << setw(12) << sampleSnps[*sample].novelTransversions;
    *ofile.outputStream << setw(12) << sampleSnps[*sample].knownTransitions;
    *ofile.outputStream << setw(12) << sampleSnps[*sample].knownTransversions;
    *ofile.outputStream << setw(12) << setprecision(3) << dbsnp;
    *ofile.outputStream << setw(8) << setprecision(3) << tstv;
    *ofile.outputStream << setw(8) << setprecision(3) << noveltstv;
    *ofile.outputStream << setw(8) << setprecision(3) << knowntstv;
    *ofile.outputStream << setw(12) << sampleSnps[*sample].homRef;
    *ofile.outputStream << setw(12) << sampleSnps[*sample].het;
    *ofile.outputStream << setw(12) << sampleSnps[*sample].homAlt;
    *ofile.outputStream << setw(12) << sampleSnps[*sample].unknown;
    *ofile.outputStream << setw(12) << sampleSnps[*sample].singletons;
    *ofile.outputStream << setw(12) << setprecision(3) << depth;
    *ofile.outputStream << setw(12) << setprecision(3) << altDepth;
    *ofile.outputStream << endl;
  }
  *ofile.outputStream << endl;
}
