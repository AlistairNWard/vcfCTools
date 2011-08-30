// ******************************************************
// vcfCTools (c) 2011 Alistair Ward
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ------------------------------------------------------
// Last modified: 18 February 2011
// ------------------------------------------------------
// vcfClass describes the vcf class and all operations.
// ******************************************************

#include "variant.h"

using namespace std;
using namespace vcfCTools;

//Constructor.
variant::variant(void) {
  isDbsnp         = false;
  processSnps     = false;
  processMnps     = false;
  processIndels   = false;
  processAll      = false;
  recordsInMemory = 1000;
  removeGenotypes = false;
  splitMnps       = false;
};

// Destructor.
variant::~variant(void) {};

// Determine which variants to process.  This must be all variant types
// or only one, no other combination.  If all variants are processed,
// the output will contain the variants written as they were in the input
// vcf.  If a single variant type is being considered, the shortest alleles
// will be written to the output.  For example, if a SNP and MNP are present
// at the same locus, AG -> TG,TC, choosing SNPs only will cause the SNP to
// be written as A -> T.  This could also cause the coordinate to change.
void variant::determineVariantsToProcess(bool snps, bool mnps, bool indels, bool isSplitMnps, bool alleles, bool store) {
  if (snps) {processSnps = true;}
  if (mnps) {processMnps = true;}
  if (indels) {processIndels = true;}
  if (!snps && !mnps && !indels) {
    processSnps   = true;
    processMnps   = true;
    processIndels = true;
    processAll    = true;
  }

  // If the MNPs are to be broken into SNPs, set the flag
  if (isSplitMnps) {splitMnps = true;}

  // Set the assessAltAlleles flag.  This determines whether the alt alleles
  // should be reduced to their minimum size (e.g. CC -> CG reduced to C -> G).
  // For operations such as validation or merging, it is unnecessary to query
  // the alleles.  For operations such as performing intersections it is
  // necessary to populate a data structure with the new allele positions, whereas
  // sometimes (e.g. filtering) it is only necessary to determine the variant type
  // (e.g. SNP/MNP/indel).  The storeReducedAlts determines whether or not to
  // populate this structure.
  assessAlts       = alleles;
  storeReducedAlts = store;
}

// Build up a structure containing variants.
bool variant::buildVariantStructure(vcf& v) {
  unsigned int count = 0;

// When variants in the correct reference sequence are found, build the variant
// structure.
  count = 0;
  string tempReferenceSequence = v.variantRecord.referenceSequence;
  while (v.success && count < recordsInMemory && v.variantRecord.referenceSequence == tempReferenceSequence) {
    addVariantToStructure(v.position, v.variantRecord);
    v.success = v.getRecord(tempReferenceSequence); 
    count++;
  }

  return v.success;
}

// Add a variant from the vcf file into the variant structure.
void variant::addVariantToStructure(int position, variantDescription& variant) {
  int altID;
  bool firstAtLocus   = true;
  vector<string> alts = split(variant.altString, ",");

  // Initialise variant type and reducedVariants structure.
  variantType type;
  reducedVariants rVar;

  // Only update the originalVariantsMap if this is the first time a
  // variant at this position is seen.
  if (originalVariantsMap.count(position) == 0) {
    originalVariantsMap[position].referenceSequence = variant.referenceSequence;
    originalVariantsMap[position].position          = position;
    originalVariantsMap[position].rsid              = variant.rsid;
    originalVariantsMap[position].ref               = variant.ref;
    originalVariantsMap[position].alts              = alts;
    originalVariantsMap[position].numberAlts        = alts.size();
    originalVariantsMap[position].quality           = variant.quality;
    originalVariantsMap[position].filters           = variant.filters;
    originalVariantsMap[position].info              = variant.info;
    originalVariantsMap[position].hasGenotypes      = variant.hasGenotypes;
    originalVariantsMap[position].genotypeFormat    = variant.genotypeFormatString;
    originalVariantsMap[position].genotypes         = variant.genotypeString;
  
    // None of the variants should be marked as filtered, so the filtered vector
    // should be set to the size of the alts vector, with false in every field.
    for (int i = 0; i < originalVariantsMap[position].alts.size(); i++) {
      originalVariantsMap[position].filtered.push_back(false);
    }
  } else {
    firstAtLocus = false;
  }

  // If the individual alternate alleles need to be examined to 
  // determine their types (SNP, MNP, indel), this should be done
  // now.
  if (assessAlts) {

    // Find the length of the first and second alternate allele.
    int altASize = alts[0].size();
    int altBSize = (alts.size() > 1) ? alts[1].size() : 0;
    int altCSize = (alts.size() > 2) ? alts[2].size() : 0;

    // Clear the structure defining the variant type.
    clearType(type);

    // Determine if there are multiple alleles.  If the variant is a triallelic SNP,
    // leave the variant as is, otherwise, put each alternate allele in the
    // structure seperately.
    //
    // Single alternate allele.
    if (alts.size() == 0) {
      determineVariantType(variant.referenceSequence, position, variant.ref, alts[0], type, 0, firstAtLocus);

    // Tri-allelic SNP.
    } else if (alts.size() == 2 && altASize == 1 && altBSize == 1) {
      altID = 0;
      type.isTriallelicSnp = true;

      // Add both alleles separately to the originalVariantsMap reduced alleles structure.
      vector<string>::iterator altIter = alts.begin();
      for (; altIter != alts.end(); altIter++) {
        if (firstAtLocus) {
          originalVariantsMap[position].type.push_back(type);
          originalVariantsMap[position].reducedRef.push_back(originalVariantsMap[position].ref);
          originalVariantsMap[position].reducedAlts.push_back(*altIter);
          originalVariantsMap[position].reducedPosition.push_back(position);
          if (position > originalVariantsMap[position].maxPosition) {originalVariantsMap[position].maxPosition = position;}
        }
        if (storeReducedAlts) {
          rVar.originalPosition = position;
          rVar.ref              = variant.ref;
          rVar.alt              = *altIter;
          rVar.altID            = altID;
          variantMap[position].referenceSequence = variant.referenceSequence;
          variantMap[position].snps.push_back(rVar);
          altID++;
        }
      }

    // Quad-allelic SNP.
    } else if (alts.size() == 3 && altASize == 1 && altBSize == 1 && altCSize == 1) {
      altID = 0;
      type.isQuadallelicSnp = true;

      // Add all alleles separately to the originalVariantsMap reduced alleles structure.
      vector<string>::iterator altIter = alts.begin();
      for (; altIter != alts.end(); altIter++) {
        if (firstAtLocus) {
          originalVariantsMap[position].type.push_back(type);
          originalVariantsMap[position].reducedRef.push_back(originalVariantsMap[position].ref);
          originalVariantsMap[position].reducedAlts.push_back(*altIter);
          originalVariantsMap[position].reducedPosition.push_back(position);
          if (position > originalVariantsMap[position].maxPosition) {originalVariantsMap[position].maxPosition = position;}
        }
        if (storeReducedAlts) {
          rVar.originalPosition = position;
          rVar.ref              = variant.ref;
          rVar.alt              = *altIter;
          rVar.altID            = altID;
          variantMap[position].referenceSequence = variant.referenceSequence;
          variantMap[position].snps.push_back(rVar);
          altID++;
        }
      }

    // MNPs, insertions and deletions.
    } else {
  
      // Keep track of the number of SNPs.  If the locus contains a multiallelic
      // SNP as well as a different variant type, the following loop will identify
      // multiple biallelic SNPs but not acknowledge that it is tri- or quad
      // allelic.  This ensures that it is.
      locusHasSnp     = false;
      locusHasTriSnp  = false;
      locusHasQuadSnp = false;

      // Determine the type of each alternate alllele.
      int count = 0;
      for (vector<string>::iterator aa = alts.begin(); aa != alts.end(); aa++) {
        clearType(type);
        determineVariantType(variant.referenceSequence, position, variant.ref, *aa, type, count, firstAtLocus);
        count++;
      }

      // If the locus contains a multi allelic SNP, modify the entries in the
      // variantType structure to reflect this.  
      //
      // Since the originalVariantsMap currently only contains the variants
      // from the first record at this locus, only do this step the first
      // time this locus is observed.
      if (firstAtLocus) {
        if (locusHasTriSnp) {
          vector<variantType>::iterator typeIter = originalVariantsMap[position].type.begin();
          for (; typeIter != originalVariantsMap[position].type.end(); typeIter++) {
            if (typeIter->isBiallelicSnp) {
              typeIter->isBiallelicSnp  = false;
              typeIter->isTriallelicSnp = true;
            }
          }
        } else if (locusHasQuadSnp) {
          vector<variantType>::iterator typeIter = originalVariantsMap[position].type.begin();
          for (; typeIter != originalVariantsMap[position].type.end(); typeIter++) {
            if (typeIter->isBiallelicSnp) {
              typeIter->isBiallelicSnp   = false;
              typeIter->isQuadallelicSnp = true;
            }
          }
        }
      }
    }
  }
}

// Initialise variant type.
void variant::clearType(variantType& type) {
  type.isBiallelicSnp   = false;    
  type.isTriallelicSnp  = false;    
  type.isQuadallelicSnp = false;    
  type.isMnp            = false;    
  type.isInsertion      = false;    
  type.isDeletion       = false;    
}

// If two vcf files are being compared and all variants from one of
// the files (for the current reference sequence) have been loaded
// into the data structures, finish parsing the second vcf file and
// flush out all variants for this reference sequence.
//
// Only the originalVariantsMap is used for building output records,
// so the variantMap can be kept clear at all times.
void variant::clearReferenceSequence(vcf& v, intFlags flags, string cRef, output& ofile, bool write) {
  while (originalVariantsMap.size() != 0) {

    // Since the vcf file to compare with has been exhausted, all of
    // the remaining variants will not intersect with any others and
    // so can be omitted unless unique records are required.
    if (!flags.annotate) {
      bool filter = flags.findUnique ? false : true;
      vector<bool>::iterator fIter = ovmIter->second.filtered.begin();
      for (; fIter != ovmIter->second.filtered.end(); fIter++) {*fIter = filter;}
    }

    // Build the output record, removing unwanted alleles and modifying the
    // genotypes if necessary and send to the output buffer.
    if (write) {buildOutputRecord(ofile);}

    // Erase the entries from the maps.
    originalVariantsMap.erase(ovmIter);
    if (variantMap.size() != 0) {variantMap.erase(vmIter);}

    if (v.variantRecord.referenceSequence == cRef && v.success) {
      addVariantToStructure(v.position, v.variantRecord);
      v.success = v.getRecord(cRef);
    }
    if (originalVariantsMap.size() != 0) {ovmIter = originalVariantsMap.begin();}
    if (variantMap.size() != 0) {vmIter = variantMap.begin();}
  }
}

// If a vcf file is being intersected with a bed file and all of the
// bed intervals have been processed, flush out the remaining variants
// in the data structure for this reference sequence.
void variant::clearReferenceSequenceBed(vcf& v, intFlags flags, string cRef, output& ofile) {
  while (originalVariantsMap.size() != 0) {
    if (flags.annotate) {}

    // Since the bed file has been exhausted, all of the remaining
    // records will not intersect with any intervals.  As such, remove
    // all variants unless findUnique = true.
    vector<bool>::iterator filtIter = ovmIter->second.filtered.begin();
    bool filter = flags.findUnique ? false : true;
    for (; filtIter != ovmIter->second.filtered.end(); filtIter++) {*filtIter = filter;}

    // Build the output record, removing unwanted alleles and modifying the
    // genotypes if necessary and send to the output buffer.
    buildOutputRecord(ofile);

    // Erase the entries from the maps.
    originalVariantsMap.erase(ovmIter);
    if (variantMap.size() != 0) {variantMap.erase(vmIter);}
    
    // Update the originalVariants structure.
    if (v.variantRecord.referenceSequence == cRef && v.success) {
      addVariantToStructure(v.position, v.variantRecord);
      v.success = v.getRecord(cRef);
    }
    if (originalVariantsMap.size() != 0) {ovmIter = originalVariantsMap.begin();}
    if (variantMap.size() != 0) {vmIter = variantMap.begin();}
  }
}

// Determine the variant class from the ref and alt alleles.
//void variant::determineVariantType(int position, string ref, string alt, variantDescription& variant, bool isDbsnp) {
void variant::determineVariantType(string refSeq, int position, string ref, string alt, variantType& type, int ID, bool first) {
  reducedVariants rVar;
  int start;
  bool doSmithWaterman = false;
  bool doTrimAlleles   = true;

  // TO MODIFY
  // Define the reference fasta file.  ***A COMMAND LINE
  // PATH SHOULD BE ALLOWED FOR THIS***
  string refFa = "/d2/data/references/build_37/human_reference_v37.fa";

  // SNP.
  if (ref.size() == 1 && (ref.size() - alt.size()) == 0) {
    type.isBiallelicSnp = true;
    if (first) {
      originalVariantsMap[position].type.push_back(type);
      originalVariantsMap[position].reducedRef.push_back(ref);
      originalVariantsMap[position].reducedAlts.push_back(alt);
      originalVariantsMap[position].reducedPosition.push_back(position);
      if (position > originalVariantsMap[position].maxPosition) {originalVariantsMap[position].maxPosition = position;}
    }
    if (storeReducedAlts) {
      rVar.originalPosition = position;
      rVar.ref   = ref;
      rVar.alt   = alt;
      rVar.altID = ID;
      variantMap[position].referenceSequence = refSeq;
      variantMap[position].snps.push_back(rVar);
    }
  } else {

    // Multi-base variants have the alt allele aligned to the ref allele
    // to unambiguously determine the variant type and the start
    // position.  This can be done with a Smith-Waterman alignment
    // between the ref and alt, or just trimming the ref and alt until
    // just matching sequence is left.
    string alRef         = ref;
    string alAlt         = alt;
    string originalRef   = ref;

    if (doSmithWaterman) {
      int pos = position;
      start = alignAlternate(refSeq, pos, ref, alt, alRef, alAlt, refFa); // vcf_aux.cpp
      if (pos != start) {
        cerr << "WARNING: Modified variant position from " << refSeq;
        cerr << ":" << pos << " to " << refSeq << ":" << start << endl;
      }
    } else if (doTrimAlleles) {
      start = trimAlleles(refSeq, position, ref, alt, alRef, alAlt); // trim_alleles.cpp
      if (start != position) {
        cerr << "WARNING: Modified variant position from " << refSeq;
        cerr << ":" << position << " to " << refSeq << ":" << start << endl;
      }
    }

    // Populate the structure rVar with the modified variant.
    rVar.originalPosition = position;
    rVar.ref              = alRef;
    rVar.alt              = alAlt;
    rVar.altID            = ID;

    // SNP.
    if (alRef.size() == 1 && (alRef.size() - alAlt.size()) == 0) {
      type.isBiallelicSnp = true;

      // Check if other SNPs have been seen at this locus.
      if (locusHasQuadSnp) {
        cerr << "ERROR: More than three alternate SNP alleles at a single locus." << endl;
        cerr << "Check position " << refSeq << ":" << position << endl;
        cerr << "Program terminated." << endl;
        exit(1);
      } else if (locusHasTriSnp) {
        locusHasTriSnp  = false;
        locusHasQuadSnp = true;
      } else if (locusHasSnp) {
        locusHasSnp    = false;
        locusHasTriSnp = true;
      } else {
        locusHasSnp = true;
      }
      if (first) {
        originalVariantsMap[position].type.push_back(type);
        originalVariantsMap[position].reducedRef.push_back(alRef);
        originalVariantsMap[position].reducedAlts.push_back(alAlt);
        originalVariantsMap[position].reducedPosition.push_back(start);
        if (position > originalVariantsMap[position].maxPosition) {originalVariantsMap[position].maxPosition = position;}
      }
      if (storeReducedAlts) {
        variantMap[position].referenceSequence = refSeq;
        variantMap[start].snps.push_back(rVar);
      }

    // MNP.
    } else if (alRef.size() != 1 && (alRef.size() - alAlt.size()) == 0) {
      type.isMnp = true;
      if (first) {
        originalVariantsMap[position].type.push_back(type);
        originalVariantsMap[position].reducedRef.push_back(alRef);
        originalVariantsMap[position].reducedAlts.push_back(alAlt);
        originalVariantsMap[position].reducedPosition.push_back(start);
        if (position > originalVariantsMap[position].maxPosition) {originalVariantsMap[position].maxPosition = position;}
      }
      if (storeReducedAlts) {
        variantMap[position].referenceSequence = refSeq;
        variantMap[start].mnps.push_back(rVar);
      }
//      if (splitMnps) {
//        string::iterator ovmIter->second.altIter = alt.begin();
//        int variantPosition = position;
//
//        // Append the flag "FROM_MNP" to the info field so it is clear that this
//        // entry was generated from a called MNP.
//        if (variant.info == "") {variant.info = "SNP;FROM_MNP";}
//        else {variant.info += ";SNP;FROM_MNP";}
//
//        for (string::iterator reovmIter->second.fIter = ref.begin(); reovmIter->second.fIter != ref.end(); reovmIter->second.fIter++) {
//
//          // Add this part of the MNP to the map as a SNP.
//          variant.isBiallelicSnp = true;
//          variant.ref = *reovmIter->second.fIter;
//          variant.altString = *ovmIter->second.altIter;
//
//          // Add the variant to the map.
//          variantMap[variantPosition].biSnps.push_back(variant);
//          variantMap[variantPosition].hasBiallelicSnp = true;
//
//          // Update the alt allele and the coordinate of the variant.
//          variantPosition++;
//          ovmIter->second.altIter++;
//        }
//      } else {
//        variant.isMnp = true;
//        variant.ref = alRef;
//        variant.altString = alAlt;
//        variantMap[position].mnps.push_back(variant);
//        variantMap[position].hasMnp = true;
//      }

    // Insertion.
    } else if ( alAlt.size() > alRef.size() ) {
      type.isInsertion = true;
      if (first) {
        originalVariantsMap[position].type.push_back(type);
        originalVariantsMap[position].reducedRef.push_back(alRef);
        originalVariantsMap[position].reducedAlts.push_back(alAlt);
        originalVariantsMap[position].reducedPosition.push_back(start);
        if (position > originalVariantsMap[position].maxPosition) {originalVariantsMap[position].maxPosition = position;}
      }
      if (storeReducedAlts) {
        variantMap[position].referenceSequence = refSeq;
        variantMap[start].indels.push_back(rVar);
      }

    // Deletion.
    } else if ( alRef.size() > alAlt.size() ) {
      type.isDeletion = true;
      if (first) {
        originalVariantsMap[position].type.push_back(type);
        originalVariantsMap[position].reducedRef.push_back(alRef);
        originalVariantsMap[position].reducedAlts.push_back(alAlt);
        originalVariantsMap[position].reducedPosition.push_back(start);
        if (position > originalVariantsMap[position].maxPosition) {originalVariantsMap[position].maxPosition = position;}
      }
      if (storeReducedAlts) {
        variantMap[position].referenceSequence = refSeq;
        variantMap[start].indels.push_back(rVar);
      }

    // None of the known types.
    } else {
      cerr << "Unknown variant class:" << endl;
      cerr << "Coordinates: " << refSeq << ": " << position << endl;
      cerr << "Ref: " << ref << endl << "Alt: " << alt << endl;
      exit(1);
    }
  }
}

// Annotate the variants at this locus with the contents of the bed file.
void variant::annotateRecordBed(bedRecord& b) {
//  // SNPs.
//  if (processSnps) {
//    for (variantIter = vmIter->second.biSnps.begin(); variantIter != vmIter->second.biSnps.end(); variantIter++) {
//      if (variantIter->info != "" && variantIter->info != ".") {variantIter->info += ";" + b.info;}
//      else {variantIter->info = b.info;}
//      buildRecord(vmIter->first, *variantIter); // tools.cpp
//    }
//    for (variantIter = vmIter->second.multiSnps.begin(); variantIter != vmIter->second.multiSnps.end(); variantIter++) {
//      if (variantIter->info != "" && variantIter->info != ".") {variantIter->info += ";" + b.info;}
//      else {variantIter->info = b.info;}
//      buildRecord(vmIter->first, *variantIter); // tools.cpp
//    }
//  }
// 
//  // MNPs.
//  if (processMnps) {
//    for (variantIter = vmIter->second.mnps.begin(); variantIter != vmIter->second.mnps.end(); variantIter++) {
//      if (variantIter->info != "" && variantIter->info != ".") {variantIter->info += ";" + b.info;}
//      else {variantIter->info = b.info;}
//      buildRecord(vmIter->first, *variantIter); // tools.cpp
//    }
//  }
// 
//  // Indels.
//  if (processIndels) {
//    for (variantIter = vmIter->second.indels.begin(); variantIter != vmIter->second.indels.end(); variantIter++) {
//      if (variantIter->info != "" && variantIter->info != ".") {variantIter->info += ";" + b.info;}
//      else {variantIter->info = b.info;}
//      buildRecord(vmIter->first, *variantIter); // tools.cpp
//    }
//  }
}

// Extract the genotypes from each sample.
vector<string> variant::extractGenotypeField(string field) {
  unsigned int i;
  vector<string> values;

//  vector<string> format = split(variantIter->genotypeFormatString, ':');
//  for (i = 0; i < format.size(); i++) {
//    if (format[i] == field) {break;}
//  }
//
//  // Parse the genotype of each sample;
//  vector<string> genotypeString = split(variantIter->genotypeString, "\t");
//  for (vector<string>::iterator gIter = genotypeString.begin(); gIter != genotypeString.end(); gIter++) {
//    vector<string> genotypeFields = split(*gIter, ':');
//    if ( genotypeFields.size() < format.size() ) {values.push_back("0");}
//    else {values.push_back(genotypeFields[i]);}
//  } 

  return values;
}

// From the intersection routine, two variants are found at the same position.
// It is possible that different variants are to be compared with each other
// (e.g. SNPs with indels) and this requires looking ahead in the variant
// structures and some of the variants have variable length.  The following 
// routines deal with the comparisons and then write out the required variants
// (annotated if necessary).
//
// If the variants compared are at the same position.
void variant::compareVariantsSameLocus(variant& var, intFlags flags) {

  // Compare SNPs.
  compareAlleles(vmIter->second.snps, var.vmIter->second.snps, flags, var);

  // Compare MNPs.
  compareAlleles(vmIter->second.mnps, var.vmIter->second.mnps, flags, var);

  // Compare indels.
  compareAlleles(vmIter->second.indels, var.vmIter->second.indels, flags, var);
}

// Compare two arrays of variant alleles of the same type (e.g. all SNPs).
void variant::compareAlleles(vector<reducedVariants>& alleles1, vector<reducedVariants>& alleles2, intFlags flags, variant& var) {
  string rsid;
  vector<reducedVariants>::iterator iter;
  vector<reducedVariants>::iterator compIter;
  bool commonType;
  bool commonAlleles;

  if (alleles1.size() != 0) {
    iter = alleles1.begin();
    for (; iter != alleles1.end(); iter++) {
      commonType    = false;
      commonAlleles = false;
      if (alleles2.size() != 0) {
        compIter  = alleles2.begin();
        commonType = true;
        for (; compIter != alleles2.end(); compIter++) {

          // If the two files share the alleles, keep them only if the common
          // alleles (or union) is required.
          if (iter->alt == compIter->alt && iter->ref == compIter->ref) {
            if (flags.annotate) {rsid = var.originalVariantsMap[compIter->originalPosition].rsid;}
            commonAlleles = true;
            break;
          }
        }
      }

      // Based on whether this variant is common or not, decide whether
      // or not to filter out the allele.
      if (flags.annotate) {
        annotateRecordVcf(var.isDbsnp, iter->originalPosition, rsid, commonType, commonAlleles);
      } else {
        if (commonAlleles) {originalVariantsMap[iter->originalPosition].filtered[iter->altID] = (flags.findUnique) ? true : false;}
        else {originalVariantsMap[iter->originalPosition].filtered[iter->altID] = (flags.findUnique) ? false : true;}
      }
    }
  }
}

// Annotate the variants at this locus with the contents of the vcf file.
void variant::annotateRecordVcf(bool isDbsnp, int position, string& rsid, bool commonType, bool commonAlleles) {

  // If the vcf file used for annotation is a dbSNP vcf file, only compare the
  // relevant variant classes.  Also, parse the info string to check that the
  // variant class as determined by vcfCTools agrees with that listed in the
  // variant class (VC) info field.
  if (isDbsnp) {
    if (commonAlleles) {originalVariantsMap[position].rsid = rsid;}
  } else {

  }

// If the vcf file used for annotation is a dbSNP vcf file, only compare the
// relevant variant classes.  Also, parse the info string to check that the
// variant class as determined by vcfCTools agrees with that listed in the
// variant class (VC) info field.
//  if (isDbsnp) {
//
//    // Annotate SNPs.
//    if (vmIter->second.biSnps.size() != 0) {
//      for (iter = v.biSnps.begin(); iter != v.biSnps.end(); iter++) {
//  
//        // Check that the VC class lists this entry as a SNP.
//        annInfo.processInfoFields(iter->info);
//        annInfo.getInfo(string("VC"), v.referenceSequence, vmIter->first);
//        if (annInfo.values[0] != "SNP") {
//          cerr << "WARNING: dbSNP entry listed as " << annInfo.values[0] << ", expected SNP.";
//          cerr << " Coordinate (" << vmIter->second.referenceSequence << ":" << vmIter->first << ")" << endl;
//        }
//  
//        // Annotate the SNPs with the rsid value.
//        for (variantIter = vmIter->second.biSnps.begin(); variantIter != vmIter->second.biSnps.end(); variantIter++) {
//          variantIter->rsid = iter->rsid;
//          vcfInfo.processInfoFields(variantIter->info);
//          if (vcfInfo.infoTags.count("dbSNP") == 0) {variantIter->info += ";dbSNP";}
//          buildRecord(vmIter->first, *variantIter); // tools.cpp
//        }
//      }
//    }
//
//    // Annotate multiallelic SNPs.
//    if (vmIter->second.multiSnps.size() != 0) {
//      for (iter = v.multiSnps.begin(); iter != v.multiSnps.end(); iter++) {
//  
//        // Check that the VC class lists this entry as a SNP.
//        annInfo.processInfoFields(iter->info);
//        annInfo.getInfo(string("VC"), v.referenceSequence, vmIter->first);
//        if (annInfo.values[0] != "SNP") {
//          cerr << "WARNING: dbSNP entry listed as " << annInfo.values[0] << ", expected SNP.";
//          cerr << " Coordinate (" << vmIter->second.referenceSequence << ":" << vmIter->first << ")" << endl;
//        }
//  
//        // Annotate the SNPs with the rsid value.
//        for (variantIter = vmIter->second.multiSnps.begin(); variantIter != vmIter->second.multiSnps.end(); variantIter++) {
//          variantIter->rsid = iter->rsid;
//          vcfInfo.processInfoFields(variantIter->info);
//          if (vcfInfo.infoTags.count("dbSNP") == 0) {variantIter->info += ";dbSNP";}
//          buildRecord(vmIter->first, *variantIter); // tools.cpp
//        }
//      }
//    }
//
//    // Annotate MNPs.
//    if (vmIter->second.mnps.size() != 0) {
//      for (iter = v.mnps.begin(); iter != v.mnps.end(); iter++) {
//  
//        // Check that the VC class lists this entry as an MNP.
//        annInfo.processInfoFields(iter->info);
//        annInfo.getInfo(string("VC"), v.referenceSequence, vmIter->first);
//        if (annInfo.values[0] != "MULTI-BASE") {
//          cerr << "WARNING: dbSNP entry listed as " << annInfo.values[0] << ", expected MNP (MULTI-BASE).";
//          cerr << " Coordinate (" << vmIter->second.referenceSequence << ":" << vmIter->first << ")" << endl;
//        }
//  
//        // Annotate the SNPs with the rsid value.
//        for (variantIter = vmIter->second.mnps.begin(); variantIter != vmIter->second.mnps.end(); variantIter++) {
//          vcfInfo.processInfoFields(variantIter->info);
//          if (vcfInfo.infoTags.count("dbSNP") == 0) {variantIter->info += ";dbSNP";}
//          variantIter->rsid = iter->rsid;
//          buildRecord(vmIter->first, *variantIter); // tools.cpp
//        }
//      }
//    }
//
//    // Annotate indels.
//    if (vmIter->second.indels.size() != 0) {
//      for (iter = v.indels.begin(); iter != v.indels.end(); iter++) {
//  
//        // Check that the VC class lists this entry as an indel.
//        annInfo.processInfoFields(iter->info);
//        annInfo.getInfo(string("VC"), v.referenceSequence, vmIter->first);
//        if (annInfo.values[0] != "INDEL") {
//          cerr << "WARNING: dbSNP entry listed as " << annInfo.values[0] << ", expected INDEL.";
//          cerr << " Coordinate (" << vmIter->second.referenceSequence << ":" << vmIter->first << ")" << endl;
//        }
//  
//        // Annotate the SNPs with the rsid value.
//        for (variantIter = vmIter->second.indels.begin(); variantIter != vmIter->second.indels.end(); variantIter++) {
//          vcfInfo.processInfoFields(variantIter->info);
//          if (vcfInfo.infoTags.count("dbSNP") == 0) {variantIter->info += ";dbSNP";}
//          variantIter->rsid = iter->rsid;
//          buildRecord(vmIter->first, *variantIter); // tools.cpp
//        }
//      }
//    }
//
//// If a vcf file is provided for performing annotations that is not a dbSNP
//// vcf file, add the info string from this file to the vcf file being
//// annotated.
//  } else {
//
//    // Annotate SNPs.
//    for (iter = v.biSnps.begin(); iter != v.biSnps.end(); iter++) {
//
//      // Add the filter string from the annotation file to the info string of
//      // the vcf file.  First check that this does not already exist to avoid
//      // multiple instances of the same flag.
//      for (variantIter = vmIter->second.biSnps.begin(); variantIter != vmIter->second.biSnps.end(); variantIter++) {
//        vcfInfo.processInfoFields(variantIter->info);
//        if (vcfInfo.infoTags.count(iter->filters) == 0) {variantIter->info += ";" + iter->filters;}
//        buildRecord(vmIter->first, *variantIter);
//      }
//    }
//
//    // Annotate multiallelic SNPs.
//    for (iter = v.multiSnps.begin(); iter != v.multiSnps.end(); iter++) {
//
//      // Add the filter string from the annotation file to the info string of
//      // the vcf file.  First check that this does not already exist to avoid
//      // multiple instances of the same flag.
//      for (variantIter = vmIter->second.multiSnps.begin(); variantIter != vmIter->second.multiSnps.end(); variantIter++) {
//        vcfInfo.processInfoFields(variantIter->info);
//        if (vcfInfo.infoTags.count(iter->filters) == 0) {variantIter->info += ";" + iter->filters;}
//        buildRecord(vmIter->first, *variantIter);
//      }
//    }
//
//    // Annotate MNPs.
//    for (iter = v.mnps.begin(); iter != v.mnps.end(); iter++) {
//
//      // Add the filter string from the annotation file to the info string of
//      // the vcf file.  First check that this does not already exist to avoid
//      // multiple instances of the same flag.
//      for (variantIter = vmIter->second.mnps.begin(); variantIter != vmIter->second.mnps.end(); variantIter++) {
//        vcfInfo.processInfoFields(variantIter->info);
//        if (vcfInfo.infoTags.count(iter->filters) == 0) {variantIter->info += ";" + iter->filters;}
//        buildRecord(vmIter->first, *variantIter);
//      }
//    }
//
//    // Annotate indels.
//    for (iter = v.indels.begin(); iter != v.indels.end(); iter++) {
//
//      // Add the filter string from the annotation file to the info string of
//      // the vcf file.  First check that this does not already exist to avoid
//      // multiple instances of the same flag.
//      for (variantIter = vmIter->second.indels.begin(); variantIter != vmIter->second.indels.end(); variantIter++) {
//        vcfInfo.processInfoFields(variantIter->info);
//        if (vcfInfo.infoTags.count(iter->filters) == 0) {variantIter->info += ";" + iter->filters;}
//        buildRecord(vmIter->first, *variantIter);
//      }
//    }
//  }
}
  
// For variants that are known to be unique to a single vcf file when
// performing a comparison, set all of the variants to filtered.  This
// is only called if common alleles are requested.
void variant::filterUnique() {
  vector<reducedVariants>::iterator iter;

  // SNPs
  iter = vmIter->second.snps.begin();
  for (; iter != vmIter->second.snps.end(); iter++) {originalVariantsMap[iter->originalPosition].filtered[iter->altID] = true;}

  // MNPs
  iter = vmIter->second.mnps.begin();
  for (; iter != vmIter->second.mnps.end(); iter++) {originalVariantsMap[iter->originalPosition].filtered[iter->altID] = true;}

  // Indels
  iter = vmIter->second.indels.begin();
  for (; iter != vmIter->second.indels.end(); iter++) {originalVariantsMap[iter->originalPosition].filtered[iter->altID] = true;}
}

// Build up the output record.  Depending on preceding actions
// this may require breaking up the genotypes and info string.
void variant::buildOutputRecord(output& ofile) {
  bool removeAllele;
  bool modified;
  bool firstAlt      = true;
  bool hasAltAlleles = false;
  int position       = ovmIter->first;
  string altAlleles  = "";
  string refAllele   = ovmIter->second.ref;;
  string filter      = ovmIter->second.filters;
  ostringstream sPosition, sQuality;

  // Loop over the alternate alleles and check if they have been filtered
  // out or not.  Reconstruct the alt field, the filter string, then info
  // string and the genotypes if any of the alleles have been removed.  If
  // only outputting a particular variant type, use the reduced alleles.
  vector<string>::iterator aIter         = ovmIter->second.alts.begin();
  vector<string>::iterator modifiedIter  = ovmIter->second.reducedAlts.begin();
  vector<int>::iterator posIter          = ovmIter->second.reducedPosition.begin();
  vector<string>::iterator refIter       = ovmIter->second.reducedRef.begin();
  vector<bool>::iterator fIter           = ovmIter->second.filtered.begin();
  vector<variantType>::iterator typeIter = ovmIter->second.type.begin();

  // Create an array of integers that represent the genotype identifiers
  // for the output record.  For example, if the original record has an
  // MNP and a SNP AG -> TC,TG, the genotypes are AG = 0, TC = 1 and
  // TG = 2.  If the MNP is filtered out, the genotypes need to be
  // modified such that AG = 0, TC = 0 and TG = 1.  This would ensure that 
  // samples with the MNP allele would just have a homozygous reference
  // genotype (as the MNP allele no longer exists in the output) and
  // samples with the SNP allele are 0/0, 0/1 or 1/1.  Keeping the original
  // genotype of 0/2 would now be nonsensical as there aren't enough
  // alternate alleles for 2 to mean anything.
  vector<int> modifiedAlleles;
  int alleleID = 1;

  // Set the first entry in modifiedAlleles to 0.  This represents the reference
  // allele.  Now the array index can be used to determine the new allele ID.
  // For example if a genotype 0/2 is encountered, this can be replaced with
  // modifiedAlleles[0]/modifiedAlleles[2].
  //
  // If the alleles weren't explicitly interrogated, the types etc. will not be known
  // and the entire record is required.
  modifiedAlleles.push_back(0);
  if (assessAlts) {
    modified = false;
    for (; aIter != ovmIter->second.alts.end(); aIter++) {

      // Check if this allele is to be kept or removed.
      removeAllele = false;
      if (*fIter) {
        removeAllele = true;
      } else {
        if ( (!processSnps && (typeIter->isBiallelicSnp || typeIter->isTriallelicSnp || typeIter->isQuadallelicSnp)) ||
             (!processMnps && typeIter->isMnp) ||
             (!processIndels && (typeIter->isInsertion || typeIter->isDeletion)) )
        {removeAllele = true;}
      }
  
      // Now update the modifiedAlleles vector depending on whether
      // the allele is to be removed or not.
      if (removeAllele) {
        modified = true;
        modifiedAlleles.push_back(-1);
      } else {
        hasAltAlleles = true;
        modifiedAlleles.push_back(alleleID);
        alleleID++;
        if (firstAlt) {firstAlt = false;}
        else {altAlleles += ",";}
        if (processSnps && !processMnps && !processIndels) {
          altAlleles += *modifiedIter;
          refAllele   = *refIter;
          position    = *posIter;
        } else {
          altAlleles += *aIter;
        }
      }
      fIter++;
      typeIter++;
      modifiedIter++;
      refIter++;
      posIter++;
    }
  } else {
    hasAltAlleles = true;
    altAlleles = *aIter;
    aIter++;
    for (; aIter != ovmIter->second.alts.end(); aIter++) {altAlleles += "," + *aIter;}
  }

  // Check if any alleles remain.  If all are filtered out, there is no
  // record to output.
  if (hasAltAlleles) {

    // The position and quality need to be converted to strings via
    // a stringstream in order to included in the string.
    sPosition << position;
    sQuality << ovmIter->second.quality;

    // If alt alleles have been removed, the info field needs to be modified
    // so that the fields with a value per alternate have the correct number
    // of entries.
    variantInfo info(ovmIter->second.info, headerInfoFields);
    if (modified) {info.modifyInfo(modifiedAlleles);}

    // Build up the standard fields from the constituent parts.
    ofile.outputRecord = ovmIter->second.referenceSequence + "	" +
                         sPosition.str() + "	" +
                         ovmIter->second.rsid + "	" +
                         refAllele + "	" +
                         altAlleles + "	" +
                         sQuality.str() + "	" +
                         filter + "	" +
                         info.infoString;

    // Now, if genotypes exist, modify them if necessary and then add
    // to the record.
    if (ovmIter->second.hasGenotypes && !removeGenotypes) {
      genotypeInfo gen(ovmIter->second.genotypeFormat, ovmIter->second.genotypes, headerFormatFields);

      // If there are multiple alleles at this locus and some of them
      // are being removed/filtered out, the sample genotypes need to be
      // modified to be consistent with the number of alternate alleles
      // being output.
      if (modified) {gen.modifyGenotypes(modifiedAlleles);}
      ofile.outputRecord += "	" + gen.genotypeFormat + "	" + gen.genotypeString;
    }

    // Flush the output record to the output buffer.
    ofile.flushToBuffer(ovmIter->first, ovmIter->second.referenceSequence);
  }
}
