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
  isDbsnp               = false;
  processAll            = false;
  processComplex        = false;
  processIndels         = false;
  processMnps           = false;
  processRearrangements = false;
  processSvs            = false;
  processSnps           = false;
  recordsInMemory       = 1000;
  removeGenotypes       = false;
  splitMnps             = false;
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
void variant::determineVariantsToProcess(bool snps, bool mnps, bool indels, bool complexV, bool svs, bool rearrangements, bool isSplitMnps, bool alleles, bool store) {
  if (snps) {processSnps = true;}
  if (mnps) {processMnps = true;}
  if (indels) {processIndels = true;}
  if (complexV) {processComplex = true;}
  if (svs) {processSvs = true;}
  if (rearrangements) {processRearrangements = true;}
  if (!snps && !mnps && !indels && !complexV && !svs && !rearrangements) {
    processAll            = true;
    processComplex        = true;
    processMnps           = true;
    processIndels         = true;
    processRearrangements = true;
    processSvs            = true;
    processSnps           = true;
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
  ov.referenceSequence = variant.referenceSequence;
  ov.position          = position;
  ov.numberAlts        = alts.size();
  ov.rsid              = variant.rsid;
  ov.ref               = variant.ref;
  ov.altString         = variant.altString;
  ov.quality           = variant.quality;
  ov.filters           = variant.filters;
  ov.info              = variant.info;
  ov.hasGenotypes      = variant.hasGenotypes;
  ov.genotypeFormat    = variant.genotypeFormatString;
  ov.genotypes         = variant.genotypeString;
  ov.maxPosition       = position;

  // Now update information based on whether this is the first record at this
  // locus.
  if (originalVariantsMap[position].size() == 0) {
    ov.hasMultipleRecords     = false;
    ov.numberOfRecordsAtLocus = 1;
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

    // MNPs, insertions, deletions and structural variants.
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

  //if (position == 137939) {exit(0);}

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
  type.isComplex        = false;
  type.isSv             = false;
  type.isRearrangement  = false;
}

// Determine the variant class from the ref and alt alleles.
//void variant::determineVariantType(int position, string ref, string alt, variantDescription& variant, bool isDbsnp) {
void variant::determineVariantType(string refSeq, int position, string ref, string alt, variantType& type, int ID, originalVariants& ov) {
  reducedVariants rVar;
  size_t containsAngleBracket   = alt.find('<');
  size_t containsSquareBracketL = alt.find('[');
  size_t containsSquareBracketR = alt.find(']');

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

  // Structural variants.
  } else if (containsAngleBracket != string::npos) {
    type.isSv = true;
    updateVariantMaps(alt, type, ref, alt, position, refSeq, ov);

    if (storeReducedAlts) {
      rVar.recordNumber     = ov.numberOfRecordsAtLocus;
      rVar.originalPosition = position;
      rVar.ref   = ref;
      rVar.alt   = alt;
      rVar.altID = ID;
      variantMap[position].svs.push_back(rVar);
    }

  // Complex rearrangments.
  } else if (containsSquareBracketL != string::npos || containsSquareBracketR != string::npos) {
    type.isRearrangement = true;
    updateVariantMaps(alt, type, ref, alt, position, refSeq, ov);

    if (storeReducedAlts) {
      rVar.recordNumber     = ov.numberOfRecordsAtLocus;
      rVar.originalPosition = position;
      rVar.ref   = ref;
      rVar.alt   = alt;
      rVar.altID = ID;
      variantMap[position].rearrangements.push_back(rVar);
    }

  // MNPs and indels.
  } else {

    // Multi-base variants have the alt allele aligned to the ref allele
    // to unambiguously determine the variant type and the start
    // position.  This can be done with a Smith-Waterman alignment
    // between the ref and alt, or just trimming the ref and alt until
    // just matching sequence is left.
    modifyAlleles mod(refSeq, position, ref, alt);

    // TO MODIFY
    // Define the reference fasta file.  ***A COMMAND LINE
    // PATH SHOULD BE ALLOWED FOR THIS***
    mod.fasta = "/d2/data/references/build_37/human_reference_v37.fa";

    mod.trim();
    if (mod.modifiedPosition != mod.originalPosition) {
      cerr << "WARNING: Modified variant position from " << refSeq;
      cerr << ":" << mod.originalPosition << " to " << refSeq << ":" << mod.modifiedPosition << endl;
    }

    // Populate the structure rVar with the modified variant.
    rVar.recordNumber     = ov.numberOfRecordsAtLocus;
    rVar.originalPosition = mod.modifiedPosition;
    rVar.ref              = mod.modifiedRef;
    rVar.alt              = mod.modifiedAlt;
    rVar.altID            = ID;

    // SNP.
    if (mod.modifiedRef.size() == 1 && (mod.modifiedRef.size() - mod.modifiedAlt.size()) == 0) {
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
      updateVariantMaps(alt, type, mod.modifiedRef, mod.modifiedAlt, mod.modifiedPosition, refSeq, ov);
      if (storeReducedAlts) {variantMap[position].snps.push_back(rVar);}

    // MNP.
    } else if (mod.modifiedRef.size() != 1 && (mod.modifiedRef.size() - mod.modifiedAlt.size()) == 0) {
      type.isMnp = true;
      updateVariantMaps(alt, type, mod.modifiedRef, mod.modifiedAlt, mod.modifiedPosition, refSeq, ov);
      if (storeReducedAlts) {variantMap[position].mnps.push_back(rVar);}

    // Indels are checked to ensure that they are left aligned.  A variant
    // is considered to be an insertion if the reference allele is a
    // single base, the alternate allele is greater that one base and both
    // alleles share the first base.
    //
    // Insertion.
    } else if (mod.modifiedRef.size() == 1 && mod.modifiedRef[0] == mod.modifiedAlt[0]) {
      type.isInsertion = true;
      mod.type = type;
      mod.stepAlleles();
      updateVariantMaps(alt, type, mod.modifiedRef, mod.modifiedAlt, mod.modifiedPosition, refSeq, ov);
      if (mod.originalPosition != mod.modifiedPosition) {
        cerr << "WARNING: Modified insertion locus from " << refSeq;
        cerr << ":" << mod.originalPosition << " to " << refSeq << ":" << mod.modifiedPosition << endl;
      }
      if (storeReducedAlts) {variantMap[mod.modifiedPosition].insertions.push_back(rVar);}

    // Deletion.
    } else if (mod.modifiedAlt.size() == 1 && mod.modifiedRef[0] == mod.modifiedAlt[0]) {
      type.isDeletion = true;
      mod.type = type;
      mod.stepAlleles();
      if (mod.originalPosition != mod.modifiedPosition) {
        cerr << "WARNING: Modified deletion locus from " << refSeq;
        cerr << ":" << mod.originalPosition << " to " << refSeq << ":" << mod.modifiedPosition << endl;
      }

      updateVariantMaps(alt, type, mod.modifiedRef, mod.modifiedAlt, mod.modifiedPosition, refSeq, ov);
      if (storeReducedAlts) {variantMap[mod.modifiedPosition].deletions.push_back(rVar);}

    // Remaining variants are in the complex class.
    } else {
      type.isComplex = true;
      mod.type = type;
      //mod.extendAlleles();
      //mod.alignAlleles();
      updateVariantMaps(alt, type, mod.modifiedRef, mod.modifiedAlt, mod.modifiedPosition, refSeq, ov);
      if (storeReducedAlts) {variantMap[mod.modifiedPosition].complexVariants.push_back(rVar);}
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
void variant::clearOriginalVariants(vcfHeader& header, intFlags& flags, output& ofile, bool write) {
  while (originalVariantsMap.size() != 0)  {

    // If this is the second vcf file and an annotation task is being performed,
    // do not send these records to the output.
    if (write || flags.findUnion) {buildOutputRecord(ofile, header);}
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
void variant::clearReferenceSequence(vcfHeader& header, vcf& v, intFlags flags, string cRef, output& ofile, bool write) {
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
    if (write || flags.findUnion) {buildOutputRecord(ofile, header);}

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

// If a vcf file is being intersected with a bed file and all of the
// bed intervals have been processed, flush out the remaining variants
// in the data structure for this reference sequence.
void variant::clearReferenceSequenceBed(vcfHeader& header, vcf& v, intFlags flags, string cRef, output& ofile) {
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
    buildOutputRecord(ofile, header);

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
  vector<reducedVariants>::iterator iter;
  
  // Compare variant types individually.
  compareAlleles(vmIter->second.snps, var.vmIter->second.snps, flags, var);
  compareAlleles(vmIter->second.mnps, var.vmIter->second.mnps, flags, var);
  compareAlleles(vmIter->second.insertions, var.vmIter->second.insertions, flags, var);
  compareAlleles(vmIter->second.deletions, var.vmIter->second.deletions, flags, var);
  compareAlleles(vmIter->second.complexVariants, var.vmIter->second.complexVariants, flags, var);

  // Structural variation and rearrangement events are currently removed from the
  // output file when performing intersections (regardless of the actual operation).
  iter = vmIter->second.svs.begin();
  for (; iter != vmIter->second.svs.end(); iter++) {
    originalVariantsMap[iter->originalPosition][iter->recordNumber - 1].filtered[iter->altID] = true;
  }
  iter = vmIter->second.rearrangements.begin();
  for (; iter != vmIter->second.rearrangements.end(); iter++) {
    originalVariantsMap[iter->originalPosition][iter->recordNumber - 1].filtered[iter->altID] = true;
  }
}

// Compare two arrays of variant alleles of the same type (e.g. all SNPs).
void variant::compareAlleles(vector<reducedVariants>& alleles1, vector<reducedVariants>& alleles2, intFlags flags, variant& var) {
  bool write;
  string infoAdd;
  string rsid;
  vector<bool> commonA (alleles1.size(), false);
  vector<bool> commonB (alleles2.size(), false);
  vector<bool>::iterator aIter;
  vector<bool>::iterator bIter;
  vector<reducedVariants>::iterator iter;
  vector<reducedVariants>::iterator compIter;

  // If it is only necessary for the alleles to share the same starting location,
  // there is no need to loop over the alleles.  All alleles at this location
  // are for the same variant class and so shouldn't be filtered out.  If only
  // one of the files has variants of this type at this location, filter out
  // all the alleles.
  if (flags.sitesOnly) {
    if (alleles1.size() != 0 && alleles2.size() != 0) {
      aIter = commonA.begin();
      bIter = commonB.begin();
      for (; aIter != commonA.end(); aIter++) {*aIter = true;}
      for (; bIter != commonB.end(); bIter++) {*bIter = true;}
    }
  } else {
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
              if (flags.annotate) {
                if (var.isDbsnp) {
                  rsid = var.originalVariantsMap[compIter->originalPosition][compIter->recordNumber - 1].rsid;
                  infoAdd = "dbSNP";
                }
                else {
                  infoAdd = var.originalVariantsMap[compIter->originalPosition][compIter->recordNumber - 1].filters;
                }
              }
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
  }

  // Update the filtered vectors in both originalVariantsMaps (each file) as long
  // as this isn't an annotation task.
  //
  // Start with the first file.
  iter  = alleles1.begin();
  aIter = commonA.begin();
  bIter = commonB.begin();
  for (; iter != alleles1.end(); iter++) {
    if (*aIter && *bIter) {
      if (flags.annotate) {
        annotateRecordVcf(var, iter->originalPosition, iter->recordNumber - 1, rsid, infoAdd);
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
    bIter++;
  }

  // Then the second file.
  if (!flags.writeFromFirst || flags.findUnion) {
    iter  = alleles2.begin();
    aIter = commonA.begin();
    bIter = commonB.begin();
    for (; iter != alleles2.end(); iter++) {
      if (*aIter && *bIter) {
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
      bIter++;
    }
  }
}

// Annotate the variants at this locus with the contents of the vcf file.
void variant::annotateRecordVcf(variant& var, int position, unsigned int record, string& rsid, string& infoAdd) {
  string oString;

  // If the vcf file used for annotation is a dbSNP vcf file, update the rsid field.
  if (var.isDbsnp) {originalVariantsMap[position][record].rsid = rsid;}

  // If the annotation file is not a dbSNP file, add the contents of the filter field
  // from the annotation vcf to the info string of the first vcf file.
  else {
    //cerr << var.originalVariantsMap[position][record].filters << endl;
  }

  // Add the contents of the variable infoAdd to the info string.
  oString = originalVariantsMap[position][record].info;
  originalVariantsMap[position][record].info = (oString == ".") ? infoAdd : oString + ";" + infoAdd;
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

  // Insertions.
  iter = vmIter->second.insertions.begin();
  for (; iter != vmIter->second.insertions.end(); iter++) {
    originalVariantsMap[iter->originalPosition][iter->recordNumber - 1].filtered[iter->altID] = true;
  }

  // Deletions.
  iter = vmIter->second.deletions.begin();
  for (; iter != vmIter->second.deletions.end(); iter++) {
    originalVariantsMap[iter->originalPosition][iter->recordNumber - 1].filtered[iter->altID] = true;
  }

  // Complex events.
  iter = vmIter->second.complexVariants.begin();
  for (; iter != vmIter->second.complexVariants.end(); iter++) {
    originalVariantsMap[iter->originalPosition][iter->recordNumber - 1].filtered[iter->altID] = true;
  }

  // Structural variation events.
  iter = vmIter->second.svs.begin();
  for (; iter != vmIter->second.svs.end(); iter++) {
    originalVariantsMap[iter->originalPosition][iter->recordNumber - 1].filtered[iter->altID] = true;
  }

  // Complex rearrangement events.
  iter = vmIter->second.rearrangements.begin();
  for (; iter != vmIter->second.rearrangements.end(); iter++) {
    originalVariantsMap[iter->originalPosition][iter->recordNumber - 1].filtered[iter->altID] = true;
  }
}

// Build up the output record.  Depending on preceding actions this may require
// breaking up the genotypes and info string.  Also, if the locus had multiple
// records in the input vcf, output the same multiple records.
void variant::buildOutputRecord(output& ofile, vcfHeader& header) {
  bool hasAltAlleles;
  bool removedAllele;
  int alleleID;
  int position;  
  string altAlleles;
  string filter;
  string refAllele;
  vector<int> modifiedAlleles;

  // Loop over the individual records at this locus.  Merging these into a single
  // record is not implemented.
  for (ovIter = ovmIter->second.begin(); ovIter != ovmIter->second.end(); ovIter++) {

    // Reset certain variables for each iteration.
    alleleID   = 1;
    altAlleles = "";
    modifiedAlleles.clear();
    modifiedAlleles.push_back(0);
    hasAltAlleles = false;
    removedAllele = false;

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
        else if (typeIter->isComplex) {if (!processComplex) {*fIter = true;}}
        else if (typeIter->isSv) {if (!processSvs) {*fIter = true;}}
        else if (typeIter->isRearrangement) {if (!processRearrangements) {*fIter = true;}}

        if (!*fIter) {
          hasAltAlleles = true;
  
          // If only SNPs are being output, the alleles in the output record
          // should be the reduced alleles.
          if (processSnps && !processMnps && !processIndels && !processSvs && !processRearrangements) {
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
      ostringstream sPosition, sQuality;
      sPosition << position;
      sQuality << ovIter->quality;

      // If alt alleles have been removed, the info field needs to be modified
      // so that the fields with a value per alternate have the correct number
      // of entries.
      variantInfo info(ovIter->info);
      if (removedAllele) {info.modifyInfo(modifiedAlleles, header);}

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
        genotypeInfo gen(ovIter->genotypeFormat, ovIter->genotypes);

        // If there are multiple alleles at this locus and some of them
        // are being removed/filtered out, the sample genotypes need to be
        // modified to be consistent with the number of alternate alleles
        // being output.
        if (removedAllele) {gen.modifyGenotypes(header, modifiedAlleles);}
        ofile.outputRecord += "	" + gen.genotypeFormat + "	" + gen.genotypeString;
      }

      // Flush the output record to the output buffer.
      ofile.flushToBuffer(ovmIter->first, ovIter->referenceSequence);
    }
  }
}
