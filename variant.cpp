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

  // If this is not the first time that this reference sequence has been seen
  // when building, the structure then the records for this reference sequence
  // are not contiguous in the vcf file.  For a number of tools (especially
  // the intersect tool), this will cause erroneous results to be generated.
  if (referenceSequenceInfo.count(tempReferenceSequence) != 0) {
    referenceSequenceInfo[tempReferenceSequence].contiguous = false;
  }
  while (v.success && count < recordsInMemory && v.variantRecord.referenceSequence == tempReferenceSequence) {

    // Update the information about observed reference sequences.
    if (referenceSequenceInfo.count(tempReferenceSequence) == 0) {
      referenceSequenceInfo[tempReferenceSequence].numberRecords    = 1;
      referenceSequenceInfo[tempReferenceSequence].usedInComparison = false;
      referenceSequenceInfo[tempReferenceSequence].contiguous       = true;
    } else {
      referenceSequenceInfo[tempReferenceSequence].numberRecords++;
    }

    // Add the variant into the structure.
    addVariantToStructure(v.position, v.variantRecord);
    v.success = v.getRecord(); 
    count++;
  }

  return v.success;
}

// Add a variant from the vcf file into the variant structure.
void variant::addVariantToStructure(int position, variantDescription& variant) {
  int altID;
  vector<string> alts = split(variant.altString, ",");
  vector<string>::iterator altIter;

  // Initialise variant type and reducedVariants structure.
  variantType type;
  reducedVariants rVar;

  // Since there can be multiple records for the same locus, define an
  // originalVariants structure and populate it with information about
  // the current record.  When complete, push this into the vector of
  // records at this locus.
  originalVariants ov;

  // First, update all vectors whose size is the number of records at this locus.
  ov.rsid           = variant.rsid;
  ov.ref            = variant.ref;
  ov.altString      = variant.altString;
  ov.quality        = variant.quality;
  ov.filters        = variant.filters;
  ov.info           = variant.info;
  ov.hasGenotypes   = variant.hasGenotypes;
  ov.genotypeFormat = variant.genotypeFormatString;
  ov.genotypes      = variant.genotypeString;
  ov.maxPosition    = position;

  // Now update information based on whether this is the first record at this
  // locus.
  if (originalVariantsMap[position].size() == 0) {
    ov.hasMultipleRecords     = false;
    ov.numberOfRecordsAtLocus = 1;

    // The reference sequence and position will be the same for all records at
    // this locus and so are not stored in vectors.
    ov.referenceSequence = variant.referenceSequence;
    ov.position          = position;
    ov.numberAlts        = alts.size();
  
  } else {
    ov.hasMultipleRecords     = true;
    ov.numberOfRecordsAtLocus = originalVariantsMap[position].size() + 1;
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
      determineVariantType(variant.referenceSequence, position, variant.ref, alts[0], type, 0, ov);

    // Tri-allelic SNP.
    } else if (alts.size() == 2 && altASize == 1 && altBSize == 1) {
      altID = 0;
      type.isTriallelicSnp = true;

      // Add both alleles separately to the originalVariantsMap reduced alleles structure.
      for (altIter = alts.begin(); altIter != alts.end(); altIter++) {
        updateVariantMaps(*altIter, type, variant.ref, *altIter, position, variant.referenceSequence, ov);

        // Now store the reduced description in the secondary structure if required.
        if (storeReducedAlts) {
          rVar.recordNumber     = ov.numberOfRecordsAtLocus;
          rVar.originalPosition = position;
          rVar.ref              = variant.ref;
          rVar.alt              = *altIter;
          rVar.altID            = altID;
          variantMap[position].snps.push_back(rVar);
          altID++;
        }
      }

    // Quad-allelic SNP.
    } else if (alts.size() == 3 && altASize == 1 && altBSize == 1 && altCSize == 1) {
      altID = 0;
      type.isQuadallelicSnp = true;

      // Add all alleles separately to the originalVariantsMap reduced alleles structure.
      for (altIter = alts.begin(); altIter != alts.end(); altIter++) {
        updateVariantMaps(*altIter, type, variant.ref, *altIter, position, variant.referenceSequence, ov);
        if (storeReducedAlts) {
          rVar.recordNumber     = ov.numberOfRecordsAtLocus;
          rVar.originalPosition = position;
          rVar.ref              = variant.ref;
          rVar.alt              = *altIter;
          rVar.altID            = altID;
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
        determineVariantType(variant.referenceSequence, position, variant.ref, *aa, type, count, ov);
        count++;
      }

      // If the locus contains multiple SNP alleles, modify the type vector to reflect
      // that each SNP allele actually belongs to a multiallelic SNP.
      if (locusHasTriSnp) {
        vector<variantType>::iterator typeIter = ov.type.begin();
        for (; typeIter != ov.type.end(); typeIter++) {
          if (typeIter->isBiallelicSnp) {
            typeIter->isBiallelicSnp  = false;
            typeIter->isTriallelicSnp = true;
          }
        }
      } else if (locusHasQuadSnp) {
        vector<variantType>::iterator typeIter = ov.type.begin();
        for (; typeIter != ov.type.end(); typeIter++) {
          if (typeIter->isBiallelicSnp) {
            typeIter->isBiallelicSnp   = false;
            typeIter->isQuadallelicSnp = true;
          }
        }
      }
    }
  }

  // Update the originalVariantsMap.
  originalVariantsMap[position].push_back(ov);

  // The maxPosition value is used by the intersection routine to determine
  // if all variants have been processed.  If there are multiple records at
  // the same locus, the only value of maxPosition that is checked is that
  // associated with the first record at this locus.  Thus, if this isn't
  // the first record at this locus and the value of maxPosition is greater
  // than that of the first record, modify the first record to have this
  // value.
  if (ov.hasMultipleRecords) {
    for (unsigned int i = 0; i < ov.numberOfRecordsAtLocus; i++) {
      if (originalVariantsMap[position][i].maxPosition > originalVariantsMap[position][0].maxPosition) {
        originalVariantsMap[position][0].maxPosition = originalVariantsMap[position][i].maxPosition;
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

// Determine the variant class from the ref and alt alleles.
//void variant::determineVariantType(int position, string ref, string alt, variantDescription& variant, bool isDbsnp) {
void variant::determineVariantType(string refSeq, int position, string ref, string alt, variantType& type, int ID, originalVariants& ov) {
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
    updateVariantMaps(alt, type, ref, alt, position, refSeq, ov);

    if (storeReducedAlts) {
      rVar.recordNumber     = ov.numberOfRecordsAtLocus;
      rVar.originalPosition = position;
      rVar.ref   = ref;
      rVar.alt   = alt;
      rVar.altID = ID;
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
    rVar.recordNumber     = ov.numberOfRecordsAtLocus;
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

      // Update the structures.
      updateVariantMaps(alt, type, alRef, alAlt, start, refSeq, ov);
      if (storeReducedAlts) {variantMap[start].snps.push_back(rVar);}

    // MNP.
    } else if (alRef.size() != 1 && (alRef.size() - alAlt.size()) == 0) {
      type.isMnp = true;
      updateVariantMaps(alt, type, alRef, alAlt, start, refSeq, ov);
      if (storeReducedAlts) {variantMap[start].mnps.push_back(rVar);}

    // Insertion.
    } else if ( alAlt.size() > alRef.size() ) {
      type.isInsertion = true;
      updateVariantMaps(alt, type, alRef, alAlt, start, refSeq, ov);
      if (storeReducedAlts) {variantMap[start].indels.push_back(rVar);}

    // Deletion.
    } else if ( alRef.size() > alAlt.size() ) {
      type.isDeletion = true;
      updateVariantMaps(alt, type, alRef, alAlt, start, refSeq, ov);
      if (storeReducedAlts) {variantMap[start].indels.push_back(rVar);}

    // None of the known types.
    } else {
      cerr << "Unknown variant class:" << endl;
      cerr << "Coordinates: " << refSeq << ": " << position << endl;
      cerr << "Ref: " << ref << endl << "Alt: " << alt << endl;
      exit(1);
    }
  }
}

// Update the variant maps with the information about individual
// alternate alleles.
void variant::updateVariantMaps(string alt, variantType type, string alRef, string alAlt, int position, string refSeq, originalVariants& ov) {
  ov.alts.push_back(alt);
  ov.filtered.push_back(false);
  ov.type.push_back(type);
  ov.reducedRef.push_back(alRef);
  ov.reducedAlts.push_back(alAlt);
  ov.reducedPosition.push_back(position);

  if (position > ov.maxPosition) {ov.maxPosition = position;}
  if (storeReducedAlts) {variantMap[position].referenceSequence = refSeq;}
}

// Loop over all variants contained in the originalVariantsMap and send them to
// the output.
void variant::clearOriginalVariants(intFlags& flags, output& ofile, bool write) {
  while (originalVariantsMap.size() != 0)  {

    // If this is the second vcf file and an annotation task is being performed,
    // do not send these records to the output.
    if (write || flags.findUnion) {buildOutputRecord(ofile);}
    originalVariantsMap.erase(ovmIter);
    if (originalVariantsMap.size() != 0) {ovmIter = originalVariantsMap.begin();}
  }
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
    // so can be omitted unless unique/union of records is required.
    if (!flags.annotate) {
      bool filter = (!flags.findCommon) ? false : true;
      for (ovIter = ovmIter->second.begin(); ovIter != ovmIter->second.end(); ovIter++) {
        vector<bool>::iterator fIter = ovIter->filtered.begin();
        for (; fIter != ovIter->filtered.end(); fIter++) {*fIter = filter;}
      }
    }

    // Build the output record, removing unwanted alleles and modifying the
    // genotypes if necessary and send to the output buffer.
    if (write || flags.findUnion) {buildOutputRecord(ofile);}

    // Erase the entries from the maps.
    originalVariantsMap.erase(ovmIter);
    if (variantMap.size() != 0) {variantMap.erase(vmIter);}

    if (v.variantRecord.referenceSequence == cRef && v.success) {
      addVariantToStructure(v.position, v.variantRecord);
      v.success = v.getRecord();
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
    bool filter = flags.findUnique ? false : true;
    for (ovIter = ovmIter->second.begin(); ovIter != ovmIter->second.end(); ovIter++) {
      vector<bool>::iterator filtIter = ovIter->filtered.begin();
      for (; filtIter != ovIter->filtered.end(); filtIter++) {*filtIter = filter;}
    }

    // Build the output record, removing unwanted alleles and modifying the
    // genotypes if necessary and send to the output buffer.
    buildOutputRecord(ofile);

    // Erase the entries from the maps.
    originalVariantsMap.erase(ovmIter);
    if (variantMap.size() != 0) {variantMap.erase(vmIter);}
    
    // Update the originalVariants structure.
    if (v.variantRecord.referenceSequence == cRef && v.success) {
      addVariantToStructure(v.position, v.variantRecord);
      v.success = v.getRecord();
    }
    if (originalVariantsMap.size() != 0) {ovmIter = originalVariantsMap.begin();}
    if (variantMap.size() != 0) {vmIter = variantMap.begin();}
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
  bool write;
  string rsid;
  vector<bool> commonA (alleles1.size(), false);
  vector<bool> commonB (alleles2.size(), false);
  vector<bool>::iterator aIter;
  vector<bool>::iterator bIter;
  vector<reducedVariants>::iterator iter;
  vector<reducedVariants>::iterator compIter;

  if (alleles1.size() != 0) {
    iter  = alleles1.begin();
    aIter = commonA.begin();
    for (; iter != alleles1.end(); iter++) {
      if (alleles2.size() != 0) {
        compIter = alleles2.begin();
        bIter    = commonB.begin();
        for (; compIter != alleles2.end(); compIter++) {

          // If the two files share the alleles, keep them only if the common
          // alleles (or union) is required.
          if (iter->alt == compIter->alt && iter->ref == compIter->ref) {
            if (flags.annotate && var.isDbsnp) {rsid = var.originalVariantsMap[compIter->originalPosition][compIter->recordNumber - 1].rsid;}
            *aIter = true;
            *bIter = true;
            break;
          }
          bIter++;
        }
      }
      aIter++;
    }
  }

  // Update the filtered vectors in both originalVariantsMaps (each file) as long
  // as this isn't an annotation task.
  //
  // Start with the first file.
  iter  = alleles1.begin();
  aIter = commonA.begin();
  for (; iter != alleles1.end(); iter++) {
    if (*aIter) {
      if (flags.annotate) {
        annotateRecordVcf(var.isDbsnp, iter->originalPosition, iter->recordNumber - 1, rsid, true, *aIter);
      } else {
        write = (flags.findUnique) ? true : !flags.writeFromFirst;
        originalVariantsMap[iter->originalPosition][iter->recordNumber - 1].filtered[iter->altID] = write;
      }
    } else {

      // Determine if this allele is to be filtered.
      if (flags.findCommon && !flags.annotate) {
        write = true;
      } else if (flags.findUnique) {
        write = !flags.writeFromFirst;
      } else {
        write = false;
      }
      originalVariantsMap[iter->originalPosition][iter->recordNumber - 1].filtered[iter->altID] = write;
    }
    aIter++;
  }

  // Then the second file.
  if (!flags.writeFromFirst || flags.findUnion) {
    iter  = alleles2.begin();
    aIter = commonB.begin();
    for (; iter != alleles2.end(); iter++) {
      if (*aIter) {
        write = (flags.findUnique) ? true : flags.writeFromFirst;
        var.originalVariantsMap[iter->originalPosition][iter->recordNumber - 1].filtered[iter->altID] = write;
      } else {

        // Determine if this allele is to be filtered.
        if (flags.findCommon) {
          write = true;
        } else if (flags.findUnique) {
          write = flags.writeFromFirst;
        } else if (flags.findUnion) {
          write = false;
        }
        var.originalVariantsMap[iter->originalPosition][iter->recordNumber - 1].filtered[iter->altID] = write;
      }
      aIter++;
    }
  }
}

// Annotate the variants at this locus with the contents of the vcf file.
void variant::annotateRecordVcf(bool isDbsnp, int position, unsigned int record, string& rsid, bool commonType, bool commonAlleles) {

  // If the vcf file used for annotation is a dbSNP vcf file, only compare the
  // relevant variant classes.  Also, parse the info string to check that the
  // variant class as determined by vcfCTools agrees with that listed in the
  // variant class (VC) info field.
  if (isDbsnp) {
    if (commonAlleles) {originalVariantsMap[position][record].rsid = rsid;}
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

  // SNPs.
  iter = vmIter->second.snps.begin();
  for (; iter != vmIter->second.snps.end(); iter++) {
    originalVariantsMap[iter->originalPosition][iter->recordNumber - 1].filtered[iter->altID] = true;
  }

  // MNPs.
  iter = vmIter->second.mnps.begin();
  for (; iter != vmIter->second.mnps.end(); iter++) {
    originalVariantsMap[iter->originalPosition][iter->recordNumber - 1].filtered[iter->altID] = true;
  }

  // Indels.
  iter = vmIter->second.indels.begin();
  for (; iter != vmIter->second.indels.end(); iter++) {
    originalVariantsMap[iter->originalPosition][iter->recordNumber - 1].filtered[iter->altID] = true;
  }
}

// Build up the output record.  Depending on preceding actions this may require
// breaking up the genotypes and info string.  Also, if the locus had multiple
// records in the input vcf, output the same multiple records.
void variant::buildOutputRecord(output& ofile) {
  bool hasAltAlleles = false;
  bool removedAllele = false;
  int alleleID       = 1;
  int position;  
  string altAlleles  = "";
  string filter;
  string refAllele;
  ostringstream sPosition, sQuality;
  vector<int> modifiedAlleles;
  modifiedAlleles.push_back(0);

  // Loop over the individual records at this locus.  Merging these into a single
  // record is not implemented.
  for (ovIter = ovmIter->second.begin(); ovIter != ovmIter->second.end(); ovIter++) {
    vector<string>::iterator aIter = ovIter->alts.begin();

    // Define the reference allele and the position.
    filter    = ovIter->filters;
    position  = ovmIter->first;
    refAllele = ovIter->ref;

    // If the individual alternate alleles were not interrogated, the variant types
    // will be unknown etc.  In this case, reconstruct the record as it appeared in
    // the original vcf file.  Otherwise, check each alternate allele to see if it
    // needs to be removed and consequent info and genotype modifications are
    // required.
    if (!assessAlts) {
      hasAltAlleles = true;
      altAlleles = ovIter->altString + ",";
    } else { 

      // Define iterators.
      vector<bool>::iterator fIter           = ovIter->filtered.begin();
      vector<string>::iterator modifiedIter  = ovIter->reducedAlts.begin();
      vector<int>::iterator posIter          = ovIter->reducedPosition.begin();
      vector<string>::iterator refIter       = ovIter->reducedRef.begin();
      vector<variantType>::iterator typeIter = ovIter->type.begin();
  
      // Loop over all the alternate alleles and determine if they need to be removed
      // or not.  Build up the new alt allele string and generate a new array of allele
      // IDs for use in the genotypes.  For example if there are three alt alleles, the
      // values 0, 1, 2 and 3 will appear in the genotypes.  If the allele represented
      // by 2 is removed, all genotypes containing 2 should be replaced with a '.' and
      // all genotypes with a 3 should be replaced with a 2 as there are now only 2
      // alternate alleles.
      for (; aIter != ovIter->alts.end(); aIter++) {

        // Based on the variant type, check if the allele should be kept.
        if (typeIter->isBiallelicSnp || typeIter->isTriallelicSnp || typeIter->isQuadallelicSnp) {if (!processSnps) {*fIter = true;}}
        else if (typeIter->isMnp) {if (!processMnps) {*fIter = true;}}
        else if (typeIter->isInsertion || typeIter->isDeletion) {if (!processIndels) {*fIter = true;}}

        if (!*fIter) {
          hasAltAlleles = true;
  
          // If only SNPs are being output, the alleles in the output record
          // should be the reduced alleles.
          if (processSnps && !processMnps && !processIndels) {
            altAlleles += *modifiedIter + ",";
            position    = *posIter;
            refAllele   = *refIter;
          } else {
            altAlleles += *aIter + ",";
          }
          modifiedAlleles.push_back(alleleID);
          alleleID++;
        } else {
          removedAllele = true;
          modifiedAlleles.push_back(-1);
        }
        fIter++;
        modifiedIter++;
        posIter++;  
        refIter++;
        typeIter++;
      }
    }

    // Check if any alleles remain.  If all are filtered out, there is no
    // record to output.
    if (hasAltAlleles) {

      // The position and quality need to be converted to strings via
      // a stringstream in order to included in the string.
      sPosition << position;
      sQuality << ovIter->quality;

      // If alt alleles have been removed, the info field needs to be modified
      // so that the fields with a value per alternate have the correct number
      // of entries.
      variantInfo info(ovIter->info, headerInfoFields);
      if (removedAllele) {info.modifyInfo(modifiedAlleles);}

      // In constructing the alt allele string, a comma is added after every
      // alt, so the last character (the trailing comma) needs to be removed.
      altAlleles = altAlleles.substr(0, altAlleles.length() - 1);

      // Build up the standard fields from the constituent parts.
      ofile.outputRecord = ovIter->referenceSequence + "	" +
                           sPosition.str() + "	" +
                           ovIter->rsid + "	" +
                           refAllele + "	" +
                           altAlleles + "	" +
                           sQuality.str() + "	" +
                           filter + "	" +
                           info.infoString;

      // Now, if genotypes exist, modify them if necessary and then add
      // to the record.
      if (ovIter->hasGenotypes && !removeGenotypes) {
        genotypeInfo gen(ovIter->genotypeFormat, ovIter->genotypes, headerFormatFields);

        // If there are multiple alleles at this locus and some of them
        // are being removed/filtered out, the sample genotypes need to be
        // modified to be consistent with the number of alternate alleles
        // being output.
        if (removedAllele) {gen.modifyGenotypes(modifiedAlleles);}
        ofile.outputRecord += "	" + gen.genotypeFormat + "	" + gen.genotypeString;
      }

      // Flush the output record to the output buffer.
      ofile.flushToBuffer(ovmIter->first, ovIter->referenceSequence);
    }
  }
}
