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

#define SNP 1
#define TRISNP 2
#define QUADSNP 3
#define MNP 4
#define INSERTION 5
#define DELETION 6
#define COMPLEX 7

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
  generateSampleStats = false;
  generateDetailed    = false;
  sampleLevelStats.clear();
}

// Destructor.
statistics::~statistics(void) {}

// Parse the variants at this locus and generate statistics.
//void statistics::generateStatistics(variant& var, vcf& v, int position, bool useAnnotations, vector<string>& annFlags, bool generateAfs, ostream* output) {
void statistics::generateStatistics(vcfHeader& header, variant& var, bool useAnnotations, vector<string>& annFlags, bool generateAfs, output& ofile) {
  unsigned int ac;
  unsigned int variantID;
  string alleles;
  vector<unsigned int> variantIDs;

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
    variantIDs.clear();
    variantIDs.push_back(0);
    for (; refIter != var.ovIter->reducedRef.end(); refIter++) {

      // Reset some values.
      inDbsnp        = false;
      isAmination    = false;
      isDeamination  = false;
      isTransition   = false;
      isTransversion = false;

      // Biallelic SNPs.
      if (typeIter->isBiallelicSnp) {

        // Check if this variant is annotated as being in dbsnp.
        inDbsnp     = (var.ovIter->rsid == ".") ? false : true;
        hasSnp      = true;
        locusHasSnp = true;
  
        // Generate a string as a pair the pair of alleles, in lower case and in alphabetical
        // order.  A simple comparison can then be made to determine if the SNP is a 
        // transition or a transversion.
        alleles = *refIter + *altIter;
        for (int i = 0; i < 2; i++) {alleles[i] = tolower(alleles[i]);}
  
        // Retrieve information from the info fields if necessary.
        ac = 0;
        //if (generateAfs || generateSampleStats) {
          //variantInfo info(var.ovIter->info);
          //info.retrieveFields(header, true);
          //ac = atoi(info.infoFields["AC"].values[variantID].c_str());
          //info.getInfo(string("AF"), var.variantIter->referenceSequence, var.vmIter->first);
          //af = atof(info.values[0].c_str());
        //}

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
            determineSnpType(var, alleles, ac);
          }
        }

      // Insertions.
      } else if (typeIter->isInsertion) {
        insertionSize = altIter->size() - refIter->size();
        variants[var.ovIter->referenceSequence][var.ovIter->filters].indels[insertionSize].insertions++;
        hasInsertion      = true;
        locusHasInsertion = true;

      // Deletions.
      } else if (typeIter->isDeletion) {
        deletionSize = refIter->size() - altIter->size();
        variants[var.ovIter->referenceSequence][var.ovIter->filters].indels[deletionSize].deletions++;
        hasDeletion      = true;
        locusHasDeletion = true;
  
      // Complex variants. - NOT YET HANDLED.
      } else if (typeIter->isComplex) {

      // Unknown variant type.
      } else {
        cerr << "ERROR: Unknown variant type." << endl;
        cerr << "Variant at " << var.ovIter->referenceSequence;
        cerr << ":" << var.ovmIter->first << endl;
        exit(1);
      }

      // If sample level statistics are required, associate the variant ID (e.g.
      // the position in the alternate allele string), with the variant type.
      if (generateSampleStats) {
        if (typeIter->isBiallelicSnp) {
          variantIDs.push_back(SNP);
        } else if (typeIter->isTriallelicSnp) {
          variantIDs.push_back(TRISNP);
        } else if (typeIter->isQuadallelicSnp) {
          variantIDs.push_back(QUADSNP);
        } else if (typeIter->isMnp) {
          variantIDs.push_back(MNP);
        } else if (typeIter->isInsertion) {
          variantIDs.push_back(INSERTION);
        } else if (typeIter->isDeletion) {
          variantIDs.push_back(DELETION);
        } else if (typeIter->isComplex) {
          variantIDs.push_back(COMPLEX);
        }
      }

      // Iterate the alt allele and allele types.
      altIter++;
      typeIter++;
      variantID++;
    }
  
    // If sample level statistics are required, parse the genotype fields and
    // populate the relevant structures.  Parsing the genotypes should only
    // happen once per record as this can be very expensive.  If there are
    // multiple alleles in the record, each allele needs to be dealt with
    // in the genotypes.
    if (generateSampleStats) {parseGenotypes(header, var, variantIDs);}
  }
}

// Given the SNP alleles, determine if it is a transition/transversion and update
// all necessary statistics.
void statistics::determineSnpType(variant& var, string& alleles, unsigned int ac) {
  string sortedAlleles = alleles;

  // Sort the alleles in alphabetical order to aid in determining if the SNP is
  // an transition or a transversion.
  sort(sortedAlleles.begin(), sortedAlleles.end());

  // Transition:   A <-> G or C <-> T.
  if (sortedAlleles == "ag" || sortedAlleles == "ct") {
    isTransition = true;

    // From the unsorted alleles, determine if the event is an amination (T->C or
    // A->G) or deamination (C->T or G->A).
    if (alleles == "ct" || alleles == "ga") {
      isDeamination = true;
      if (inDbsnp) {
        variants[var.ovIter->referenceSequence][var.ovIter->filters].knownDeaminations++;
      } else {
        variants[var.ovIter->referenceSequence][var.ovIter->filters].novelDeaminations++;
      }
    } else {
      isAmination = true;
      if (inDbsnp) {
        variants[var.ovIter->referenceSequence][var.ovIter->filters].knownAminations++;
      } else {
        variants[var.ovIter->referenceSequence][var.ovIter->filters].novelAminations++;
      }
    }

    if (inDbsnp) {
      variants[var.ovIter->referenceSequence][var.ovIter->filters].knownTransitions++;
      //if (info.infoTags.count("dbSNPX") != 0) {variants[var.ovIter->referenceSequence][var.ovIter->filters].diffKnownTransitions++;}
    } else {
      variants[var.ovIter->referenceSequence][var.ovIter->filters].novelTransitions++;
    }

  // Transversion: A <-> C, A <-> T, C <-> G or G <-> T.
  } else if (sortedAlleles == "ac" || sortedAlleles == "at" || sortedAlleles == "cg" || sortedAlleles == "gt") {
    isTransversion = true;

    if (inDbsnp) {
      variants[var.ovIter->referenceSequence][var.ovIter->filters].knownTransversions++;
      //if (info.infoTags.count("dbSNPX") != 0) {variants[var.ovIter->referenceSequence][var.ovIter->filters].diffKnownTransversions++;}
    } else {
      variants[var.ovIter->referenceSequence][var.ovIter->filters].novelTransversions++;
    }
  }
}

// Parse the genotypes for each sample and update the relevant
// statistics structures.
void statistics::parseGenotypes(vcfHeader& header, variant& var, vector<unsigned int> variantIDs) {
  unsigned int depth;
  unsigned int sampleID = 0;
  string entry;
  string geno;
  double quality;
  vector<string> sampleEntries;
  genotypeInfo gen(var.ovIter->genotypeFormat, var.ovIter->genotypes);
  gen.processFormats(header);

  // Parse each sample in turn.
  vector<string>::iterator genoIter = gen.genotypes.begin();
  for (; genoIter != gen.genotypes.end(); genoIter++) {
    sampleEntries = split(*genoIter, ":");

    // Check that the number of entries is consistent with the format string.
    if (sampleEntries.size() != gen.genotypeFormats.size() && sampleEntries.size() != 1) {
      cerr << "ERROR: Number of fields in the genotype string is inconsistent with the format for sample ";
      cerr << header.samples[sampleID] << " at " << var.ovIter->referenceSequence << ":" << var.ovIter->position << "." << endl;
      exit(1);
    } else if (sampleEntries.size() == 1 && sampleEntries[0] == ".") {
      sampleLevelStats[header.samples[sampleID]].unknown++;
    } else {
      quality = (gen.genotypeFields.count("GQ") == 0) ? 0. : atof( (sampleEntries[gen.genotypeFields["GQ"].ID]).c_str() );
      if (quality >= minGenotypeQuality || gen.genotypeFields.count("GQ") == 0) {
        geno  = sampleEntries[gen.genotypeFields["GT"].ID];
        //depth = (gen.genotypeFields.count("GQ") == 0) ? 0 : atoi( (sampleEntries[gen.genotypeFields["DP"].ID]).c_str() );

        // Find the allele IDs from the genotype string.  This can be 0 for
        // reference and then any number up to the number of alternate alleles.
        // Only do this if the genotype isn't '.'.
        if (geno == ".") {
          sampleLevelStats[header.samples[sampleID]].unknown++;
        } else {
          size_t separator = geno.find('/');
          if (separator == string::npos) {
            separator = geno.find("|");
            if (separator == string::npos) {
              cerr << "ERROR: Unknown genotype separator in sample " << header.samples[sampleID];
              cerr << " at " << var.ovIter->referenceSequence << ":" << var.ovIter->position << "." << endl;
              exit(1);
            }
          }
          unsigned int idA = atoi((geno.substr(0, separator)).c_str());
          unsigned int idB = atoi((geno.substr(separator + 1, geno.length())).c_str());

          // Define a structure to hold boolean flags.  This will be used in other called
          // routines.
          statsFlags flags;
          flags.het            = false;
          flags.isAmination    = isAmination;
          flags.isDeamination  = isDeamination;
          flags.inDbsnp        = inDbsnp;
          flags.isTransition   = isTransition;
          flags.isTransversion = isTransversion;
  
          // If the allele IDs are the same, the sample is homozygous.  If they are
          // both '0', then they are homozygous reference and can be dealt with immediately.
          // For non-reference alleles, the array variantIDs can be used to determine
          // which alternate allele the ID corresponds to.
          size_t dotInGeno = geno.find('.');
          if (dotInGeno != string::npos) {
            sampleLevelStats[header.samples[sampleID]].unknown++;
          } else if (idA == idB && idA == 0) {
            sampleLevelStats[header.samples[sampleID]].homRef++;
  
          // Homozygous non-reference.
          } else if (idA == idB) {
            updateSampleLevelStats(flags, variantIDs[idA], header.samples[sampleID]);
  
          // Heterozygous with a ref allele.
          } else if (idA == 0 | idB == 0) {
            flags.het = true;
            unsigned int id = (idA == 0) ? variantIDs[idB] : variantIDs[idA];
            updateSampleLevelStats(flags, id, header.samples[sampleID]);
  
          // Heterozygous with two non-reference alleles.
          } else {
          }
  
          // Update the total depth for the sample.
          //sampleSnps[header.samples[sampleID]].totalDepth += depth;
        }
      }
    }
    sampleID++;
  }
}

// Update the sample level statistics.
void statistics::updateSampleLevelStats(statsFlags& flags, unsigned int id, string& sample) {

  // Update the statistics for the correct variant type.
  if (id == SNP) {
 
    // Heterozygous SNPs.
    if (flags.het) {

      // Known SNPs.
      if (flags.inDbsnp) {
        if (flags.isTransition) {sampleLevelStats[sample].knownHetTransitions++;}
        else if (flags.isTransversion) {sampleLevelStats[sample].knownHetTransversions++;}

        // De/aminations.
        if (flags.isAmination) {sampleLevelStats[sample].knownAminations++;}
        else if (flags.isDeamination) {sampleLevelStats[sample].knownDeaminations++;}

      // Novel SNPs.
      } else {
        if (flags.isTransition) {sampleLevelStats[sample].novelHetTransitions++;}
        else if (flags.isTransversion) {sampleLevelStats[sample].novelHetTransversions++;}

        // De/aminations.
        if (flags.isAmination) {sampleLevelStats[sample].novelAminations++;}
        else if (flags.isDeamination) {sampleLevelStats[sample].novelDeaminations++;}
      }
    } else {

      // Known SNPs.
      if (flags.inDbsnp) {
        if (flags.isTransition) {sampleLevelStats[sample].knownHomTransitions++;}
        else if (flags.isTransversion) {sampleLevelStats[sample].knownHomTransversions++;}

        // De/aminations.
        if (flags.isAmination) {sampleLevelStats[sample].knownAminations++;}
        else if (flags.isDeamination) {sampleLevelStats[sample].knownDeaminations++;}

      // Novel SNPs.
      } else {
        if (flags.isTransition) {sampleLevelStats[sample].novelHomTransitions++;}
        else if (flags.isTransversion) {sampleLevelStats[sample].novelHomTransversions++;}

        // De/aminations.
        if (flags.isAmination) {sampleLevelStats[sample].novelAminations++;}
        else if (flags.isDeamination) {sampleLevelStats[sample].novelDeaminations++;}
      }
    }
  } else if (id == TRISNP) {
  } else if (id == QUADSNP) {
  } else if (id == MNP) {
    if (flags.het) {sampleLevelStats[sample].hetMnps++;}
    else {sampleLevelStats[sample].homMnps++;}
  } else if (id == INSERTION) {
    if (flags.het) {
      sampleLevelStats[sample].hetInsertions++;
    }
    else {
      sampleLevelStats[sample].homInsertions++;
    }
  } else if (id == DELETION) {
    if (flags.het) {
      sampleLevelStats[sample].hetDeletions++;
    }
    else {
      sampleLevelStats[sample].homDeletions++;
    }
  } else if (id == COMPLEX) {
    if (flags.het) {sampleLevelStats[sample].hetComplex++;}
    else {sampleLevelStats[sample].homComplex++;}
  }
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
          printHeader(ofile, string("reference_sequence"), true, true, false);
          *ofile.outputStream << endl;
          writtenHeader = true;
        }
        printVariantStruct(ofile, iter->first, vIter->second);
      }
    }
  }
 
  *ofile.outputStream << endl;
  *ofile.outputStream << "Total_statistics (de/am_ratio is the deamination (C->T/G->A) amination (T->C/A->G) ratio)." << endl;
  *ofile.outputStream << endl;
  printHeader(ofile, string("filter"), true, true, false);
  *ofile.outputStream << endl;
  for (map<string, variantStruct>::iterator iter = totalVariants["total"].begin(); iter != totalVariants["total"].end(); iter++) {
    if (iter->first != "all" && iter->first != "PASS") {printVariantStruct(ofile, string(iter->first), iter->second);}
  }
  *ofile.outputStream << setw(22) << "";
  string buf;
  buf.assign(152, '-');
  *ofile.outputStream << buf;
  *ofile.outputStream << endl;
  string filter = "PASS";
  printVariantStruct(ofile, filter, totalVariants["total"]["PASS"]);
  filter = "Total";
  printVariantStruct(ofile, filter, totalVariants["total"]["all"]);
  *ofile.outputStream << setw(22) << "";
  *ofile.outputStream << buf;
  *ofile.outputStream << endl;
}

// Print out a header line.
void statistics::printHeader(output& ofile, string text, bool dbsnpDiff, bool multi, bool sampleLevel) {
  *ofile.outputStream << setw(22) << "";
  *ofile.outputStream << setw(84) << "--------------------------------------#SNPs---------------------------------------";
  //if (dbsnpDiff) {*ofile.outputStream << setw(18) << "";}
  //else {*ofile.outputStream << setw(12) << "";}
  *ofile.outputStream << setw(12) << "";
  *ofile.outputStream << setw(24) << "------ts/tv_ratio-----";
  *ofile.outputStream << setw(24) << "------de/am_ratio-----";
  if (sampleLevel) {
    *ofile.outputStream << setw(12) << "";
    *ofile.outputStream << setw(48) << "--------------------SNPs----------------------";
    *ofile.outputStream << setw(24) << "---------MNPs---------";
    *ofile.outputStream << setw(24) << "------Insertions------";
    *ofile.outputStream << setw(24) << "------Deletions-------";
    *ofile.outputStream << setw(24) << "----Complex Events----";
  }
  *ofile.outputStream << endl;
  *ofile.outputStream << setw(22) << text;
  *ofile.outputStream << setw(12) << "total";
  *ofile.outputStream << setw(12) << "novel_ts";
  *ofile.outputStream << setw(12) << "novel_tv";
  *ofile.outputStream << setw(12) << "known_ts";
  *ofile.outputStream << setw(12) << "known_tv";
  *ofile.outputStream << setw(12) << "amination";
  *ofile.outputStream << setw(12) << "deamination";
  //if (dbsnpDiff) {*ofile.outputStream << setw(18) << setprecision(6) << "%dbsnp (%diff)";}
  //else {*ofile.outputStream << setw(12) << "%dbsnp";}
  *ofile.outputStream << setw(12) << "%dbsnp";
  *ofile.outputStream << setw(8) << setprecision(6) << "total";
  *ofile.outputStream << setw(8) << setprecision(6) << "novel";
  *ofile.outputStream << setw(8) << setprecision(6) << "known";
  *ofile.outputStream << setw(8) << setprecision(6) << "total";
  *ofile.outputStream << setw(8) << setprecision(6) << "novel";
  *ofile.outputStream << setw(8) << setprecision(6) << "known";
  if (multi) {*ofile.outputStream << setw(8) << "multi";}
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

// Print the contents of the structure variantStruct to screen in a standard format.
void statistics::printVariantStruct(output& ofile, string filter, variantStruct& var) {
  int novel         = var.novelTransitions + var.novelTransversions;
  int known         = var.knownTransitions + var.knownTransversions;
  int diffKnown     = var.diffKnownTransitions + var.diffKnownTransversions;
  int transitions   = var.novelTransitions + var.knownTransitions;
  int transversions = var.novelTransversions + var.knownTransversions;
  int multiAllelic  = var.multiAllelic;
  int totalSnp      = novel + known + multiAllelic;
  int aminations    = var.novelAminations + var.knownAminations;
  int deaminations  = var.novelDeaminations + var.knownDeaminations;

  double dbsnp     = (totalSnp == 0) ? 0. : (100. * double(known) / double(totalSnp));
  //double dbsnpX    = (totalSnp == 0) ? 0. : (100. * double(diffKnown) / double(totalSnp));
  double tstv      = (transversions == 0) ? 0. : (double(transitions) / double(transversions));
  double noveltstv = (var.novelTransversions == 0) ? 0. : (double(var.novelTransitions) / double(var.novelTransversions));
  double knowntstv = (var.knownTransversions == 0) ? 0. : (double(var.knownTransitions) / double(var.knownTransversions));
  double totalad   = (aminations == 0) ? 0. : (double(deaminations) / double(aminations));
  double novelad   = (var.novelAminations == 0) ? 0. : (double(var.novelDeaminations) / double(var.novelAminations));
  double knownad   = (var.knownAminations == 0) ? 0. : (double(var.knownDeaminations) / double(var.knownAminations));

  *ofile.outputStream << setw(22) << filter;
  *ofile.outputStream << setw(12) << setprecision(10) << totalSnp;
  *ofile.outputStream << setw(12) << var.novelTransitions,
  *ofile.outputStream << setw(12) << var.novelTransversions;
  *ofile.outputStream << setw(12) << var.knownTransitions;
  *ofile.outputStream << setw(12) << var.knownTransversions;
  *ofile.outputStream << setw(12) << var.novelAminations + var.knownAminations;;
  *ofile.outputStream << setw(12) << var.novelDeaminations + var.knownDeaminations;;
  *ofile.outputStream << setw(12) << setprecision(3) << dbsnp;
  //*ofile.outputStream << setw(2) << " (";
  //*ofile.outputStream << setw(5) << setprecision(3) << dbsnpX;
  //*ofile.outputStream << setw(2) << ") ";
  *ofile.outputStream << setw(8) << setprecision(3) << tstv;
  *ofile.outputStream << setw(8) << setprecision(3) << noveltstv;
  *ofile.outputStream << setw(8) << setprecision(3) << knowntstv;
  *ofile.outputStream << setw(8) << setprecision(3) << totalad;
  *ofile.outputStream << setw(8) << setprecision(3) << novelad;
  *ofile.outputStream << setw(8) << setprecision(3) << knownad;
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
//  map<unsigned int, snpTypes>::iterator acsIter;
//  unsigned int transitions, transversions, novel, known;
//  double dbsnp, tstv, novelTstv, knownTstv;
//
//  *ofile.outputStream << "Statistics_by_allele_count:";
//  *ofile.outputStream << endl;
//  *ofile.outputStream << endl;
//  *ofile.outputStream << setw(18) << "";
//  *ofile.outputStream << setw(60) << "--------------------------#SNPs---------------------------";
//  *ofile.outputStream << setw(10) << "";
//  *ofile.outputStream << setw(24) << "-----ts/tv_ratio-----";
//  *ofile.outputStream << endl;
//  *ofile.outputStream << setw(16) << "allele_count";
//  *ofile.outputStream << setw(12) << "total";
//  *ofile.outputStream << setw(12) << "novel_ts";
//  *ofile.outputStream << setw(12) << "novel_tv";
//  *ofile.outputStream << setw(12) << "known_ts";
//  *ofile.outputStream << setw(12) << "known_tv";
//  *ofile.outputStream << setw(12) << "%dbsnp";
//  *ofile.outputStream << setw(8) << "total";
//  *ofile.outputStream << setw(8) << "novel";
//  *ofile.outputStream << setw(8) << "known";
//  *ofile.outputStream << endl;
//  for (acsIter = totalVariants["total"]["PASS"].acs.begin(); acsIter != totalVariants["total"]["PASS"].acs.end(); acsIter++) {
//    novel = acsIter->second.novelTransitions + acsIter->second.novelTransversions;
//    known = acsIter->second.knownTransitions + acsIter->second.knownTransversions;
//    transitions = acsIter->second.novelTransitions + acsIter->second.knownTransitions;
//    transversions = acsIter->second.novelTransversions + acsIter->second.knownTransversions;
//
//    dbsnp = ( (known + novel) == 0 ) ? 0. : 100. * (double(known) / ( double(novel) + double(known)) );
//    tstv = (transversions == 0) ? 0 : double(transitions) / double(transversions);
//    novelTstv = (acsIter->second.novelTransversions == 0) ? 0 : double(acsIter->second.novelTransitions) / double(acsIter->second.novelTransversions);
//    knownTstv = (acsIter->second.knownTransversions == 0) ? 0 : double(acsIter->second.knownTransitions) / double(acsIter->second.knownTransversions);
//    *ofile.outputStream << setw(12) << acsIter->first;
//    *ofile.outputStream << setw(12) << novel + known;
//    *ofile.outputStream << setw(12) << acsIter->second.novelTransitions;
//    *ofile.outputStream << setw(12) << acsIter->second.novelTransversions;
//    *ofile.outputStream << setw(12) << acsIter->second.knownTransitions;
//    *ofile.outputStream << setw(12) << acsIter->second.knownTransversions;
//    *ofile.outputStream << setw(12) << setprecision(3) << dbsnp;
//    *ofile.outputStream << setw(8) << setprecision(3) << tstv;
//    *ofile.outputStream << setw(8) << setprecision(3) << novelTstv;
//    *ofile.outputStream << setw(8) << setprecision(3) << knownTstv;
//    *ofile.outputStream << endl;
//  }
//  *ofile.outputStream << endl;
}

// Print out allele frequency information.
void statistics::printAfs(output& ofile) {
//  map<double, snpTypes>::iterator afsIter;
//  unsigned int transitions, transversions, novel, known;
//  double dbsnp, tstv, novelTstv, knownTstv;
//
//  *ofile.outputStream << "Statistics by allele frequency:";
//  *ofile.outputStream << endl;
//  *ofile.outputStream << endl;
//  *ofile.outputStream << setw(18) << "";
//  *ofile.outputStream << setw(60) << "--------------------------#SNPs---------------------------";
//  *ofile.outputStream << setw(10) << "";
//  *ofile.outputStream << setw(24) << "-----ts/tv_ratio-----";
//  *ofile.outputStream << endl;
//  *ofile.outputStream << setw(16) << "allele_frequency";
//  *ofile.outputStream << setw(12) << "total";
//  *ofile.outputStream << setw(12) << "novel_ts";
//  *ofile.outputStream << setw(12) << "novel_tv";
//  *ofile.outputStream << setw(12) << "known_ts";
//  *ofile.outputStream << setw(12) << "known_tv";
//  *ofile.outputStream << setw(12) << "%dbsnp";
//  *ofile.outputStream << setw(8) << "total";
//  *ofile.outputStream << setw(8) << "novel";
//  *ofile.outputStream << setw(8) << "known";
//  *ofile.outputStream << endl;
//  for (afsIter = totalVariants["total"]["PASS"].afs.begin(); afsIter != totalVariants["total"]["PASS"].afs.end(); afsIter++) {
//    novel = afsIter->second.novelTransitions + afsIter->second.novelTransversions;
//    known = afsIter->second.knownTransitions + afsIter->second.knownTransversions;
//    transitions = afsIter->second.novelTransitions + afsIter->second.knownTransitions;
//    transversions = afsIter->second.novelTransversions + afsIter->second.knownTransversions;
//
//    dbsnp = ( (known + novel) == 0 ) ? 0. : 100. * (double(known) / ( double(novel) + double(known)) );
//    tstv = (transversions == 0) ? 0 : double(transitions) / double(transversions);
//    novelTstv = (afsIter->second.novelTransversions == 0) ? 0 : double(afsIter->second.novelTransitions) / double(afsIter->second.novelTransversions);
//    knownTstv = (afsIter->second.knownTransversions == 0) ? 0 : double(afsIter->second.knownTransitions) / double(afsIter->second.knownTransversions);
//    *ofile.outputStream << setw(12) << afsIter->first;
//    *ofile.outputStream << setw(12) << novel + known;
//    *ofile.outputStream << setw(12) << afsIter->second.novelTransitions;
//    *ofile.outputStream << setw(12) << afsIter->second.novelTransversions;
//    *ofile.outputStream << setw(12) << afsIter->second.knownTransitions;
//    *ofile.outputStream << setw(12) << afsIter->second.knownTransversions;
//    *ofile.outputStream << setw(12) << setprecision(3) << dbsnp;
//    *ofile.outputStream << setw(8) << setprecision(3) << tstv;
//    *ofile.outputStream << setw(8) << setprecision(3) << novelTstv;
//    *ofile.outputStream << setw(8) << setprecision(3) << knownTstv;
//    *ofile.outputStream << endl;
//  }
//  *ofile.outputStream << endl;
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
void statistics::printSampleSnps(vcfHeader& header, vcf& v, output& ofile) {
  unsigned int totalTransitions = totalVariants["total"]["all"].novelTransitions + totalVariants["total"]["all"].knownTransitions;
  unsigned int totalTransversions = totalVariants["total"]["all"].novelTransversions + totalVariants["total"]["all"].knownTransversions;
  unsigned int totalSnps = totalTransitions + totalTransversions;

  *ofile.outputStream << endl << "variant_statistics_by_sample:" << endl;
  *ofile.outputStream << endl;
  printHeader(ofile, string("sample"), false, false, true);

  // SNPs
  *ofile.outputStream << setw(12) << "hom_ref";
  *ofile.outputStream << setw(12) << "het_ts";
  *ofile.outputStream << setw(12) << "het_tv";
  *ofile.outputStream << setw(12) << "hom_alt_ts";
  *ofile.outputStream << setw(12) << "hom_alt_tv";

  // MNPs.
  *ofile.outputStream << setw(12) << "het";
  *ofile.outputStream << setw(12) << "hom_alt";

  // Insertions.
  *ofile.outputStream << setw(12) << "het";
  *ofile.outputStream << setw(12) << "hom_alt";

  // Deletions.
  *ofile.outputStream << setw(12) << "het";
  *ofile.outputStream << setw(12) << "hom_alt";

  // Complex variants.
  *ofile.outputStream << setw(12) << "het";
  *ofile.outputStream << setw(12) << "hom_alt";

  // Unknown genotypes.
  *ofile.outputStream << setw(12) << "Unknown";
  //*ofile.outputStream << setw(12) << "Singletons";
  //*ofile.outputStream << setw(12) << "Depth";
  //*ofile.outputStream << setw(12) << "Alt_depth";
  *ofile.outputStream << endl;
  for (vector<string>::iterator sample = header.samples.begin(); sample != header.samples.end(); sample++) {

    // Calculate transition/transversion ratios etc for each sample.
    unsigned int aminations    = sampleLevelStats[*sample].novelAminations + sampleLevelStats[*sample].knownAminations;
    unsigned int deaminations  = sampleLevelStats[*sample].novelDeaminations + sampleLevelStats[*sample].knownDeaminations;
    unsigned int novelAm       = sampleLevelStats[*sample].novelAminations;
    unsigned int novelDe       = sampleLevelStats[*sample].novelDeaminations;
    unsigned int novelTs       = sampleLevelStats[*sample].novelHomTransitions + sampleLevelStats[*sample].novelHetTransitions;
    unsigned int novelTv       = sampleLevelStats[*sample].novelHomTransversions + sampleLevelStats[*sample].novelHetTransversions;
    unsigned int knownAm       = sampleLevelStats[*sample].knownAminations;
    unsigned int knownDe       = sampleLevelStats[*sample].knownDeaminations;
    unsigned int knownTs       = sampleLevelStats[*sample].knownHomTransitions + sampleLevelStats[*sample].knownHetTransitions;
    unsigned int knownTv       = sampleLevelStats[*sample].knownHomTransversions + sampleLevelStats[*sample].knownHetTransversions;
    unsigned int known         = knownTs + knownTv;
    unsigned int novel         = novelTs + novelTv;
    unsigned int transitions   = novelTs + knownTs;
    unsigned int transversions = novelTv + knownTv;
    unsigned int totalSnp      = transitions + transversions;
    //double depth      = (totalSnp == 0) ? 0. : double(sampleLevelStats[*sample].totalDepth) / double(totalSnps);
    //double altDepth   = (totalSnp == 0) ? 0. : double(sampleLevelStats[*sample].totalAltDepth) / double(totalSnp);
    double deam       = (aminations == 0) ? 0. : (double(aminations) / double(deaminations));
    double dbsnp      = (totalSnp == 0) ? 0. : (100. * double(known) / double(totalSnp));
    double tstv       = (transversions == 0) ? 0. : (double(transitions) / double(transversions));
    double noveldeam  = (novelAm == 0) ? 0. : (double(novelAm) / double(novelDe));
    double noveltstv  = (novelTv == 0) ? 0. : (double(novelTs) / double(novelTv));
    double knowndeam  = (knownAm == 0) ? 0. : (double(knownAm) / double(knownDe));
    double knowntstv  = (knownTv == 0) ? 0. : (double(knownTs) / double(knownTv));

    // Write out to the output stream.
    // Sample name.
    *ofile.outputStream << setw(22) << *sample;

    // Overall SNP statistics.
    *ofile.outputStream << setw(12) << setprecision(10) << totalSnp;
    *ofile.outputStream << setw(12) << novelTs;
    *ofile.outputStream << setw(12) << novelTv;
    *ofile.outputStream << setw(12) << knownTs;
    *ofile.outputStream << setw(12) << knownTv;
    *ofile.outputStream << setw(12) << aminations;
    *ofile.outputStream << setw(12) << deaminations;
    *ofile.outputStream << setw(12) << setprecision(3) << dbsnp;
    *ofile.outputStream << setw(8) << setprecision(3) << tstv;
    *ofile.outputStream << setw(8) << setprecision(3) << noveltstv;
    *ofile.outputStream << setw(8) << setprecision(3) << knowntstv;
    *ofile.outputStream << setw(8) << setprecision(3) << deam;
    *ofile.outputStream << setw(8) << setprecision(3) << noveldeam;
    *ofile.outputStream << setw(8) << setprecision(3) << knowndeam;

    // SNP genotype information.
    *ofile.outputStream << setw(12) << sampleLevelStats[*sample].homRef;
    *ofile.outputStream << setw(12) << sampleLevelStats[*sample].knownHetTransitions + sampleLevelStats[*sample].novelHetTransitions;
    *ofile.outputStream << setw(12) << sampleLevelStats[*sample].knownHetTransversions + sampleLevelStats[*sample].novelHetTransversions;
    *ofile.outputStream << setw(12) << sampleLevelStats[*sample].knownHomTransitions + sampleLevelStats[*sample].novelHomTransitions;
    *ofile.outputStream << setw(12) << sampleLevelStats[*sample].knownHomTransversions + sampleLevelStats[*sample].novelHomTransversions;

    // MNPs.
    *ofile.outputStream << setw(12) << sampleLevelStats[*sample].hetMnps;
    *ofile.outputStream << setw(12) << sampleLevelStats[*sample].homMnps;

    // Insertions.
    *ofile.outputStream << setw(12) << sampleLevelStats[*sample].hetInsertions;
    *ofile.outputStream << setw(12) << sampleLevelStats[*sample].homInsertions;

    // Deletions.
    *ofile.outputStream << setw(12) << sampleLevelStats[*sample].hetDeletions;
    *ofile.outputStream << setw(12) << sampleLevelStats[*sample].homDeletions;

    // Complex variants.
    *ofile.outputStream << setw(12) << sampleLevelStats[*sample].hetComplex;
    *ofile.outputStream << setw(12) << sampleLevelStats[*sample].homComplex;

    // Unknown genotypes (these contain at least one '.').
    *ofile.outputStream << setw(12) << sampleLevelStats[*sample].unknown;
//    *ofile.outputStream << setw(12) << sampleSnps[*sample].unknown;
//    *ofile.outputStream << setw(12) << sampleSnps[*sample].singletons;
//    *ofile.outputStream << setw(12) << setprecision(3) << depth;
//    *ofile.outputStream << setw(12) << setprecision(3) << altDepth;
    *ofile.outputStream << endl;
  }
  *ofile.outputStream << endl;
}
