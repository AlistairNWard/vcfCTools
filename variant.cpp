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
  recordsInMemory = 1000;
};

// Destructor.
variant::~variant(void) {};

// Determine which variants to process.
void variant::determineVariantsToProcess(bool snps, bool mnps, bool indels) {
  if (snps) {processSnps = true;}
  if (mnps) {processMnps = true;}
  if (indels) {processIndels = true;}
  if (!snps && !mnps && !indels) {
    processSnps = true;
    processMnps = true;
    processIndels = true;
  }
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

// Add a variant from the vcf file into thr variant structure.
void variant::addVariantToStructure(int position, variantDescription& variant) {
  variant.isBiallelicSnp = false;
  variant.isTriallelicSnp = false;
  variant.isQuadallelicSnp = false;
  variant.isMnp = false;
  variant.isInsertion = false;
  variant.isDeletion = false;

// If this is the first variant at this position, create a new entry in the
// information structure.
  if (variantMap.count(position) == 0) {
    variantMap[position].referenceSequence = variant.referenceSequence;
  }

// Determine if there are multiple alleles.  If the variant is a triallelic SNP,
// leave the variant as is, otherwise, put each alternate allele in the
// structure seperately.
  size_t found = variant.altString.find(",");

  // Single alternate allele.
  if (found == string::npos) {
    determineVariantClass(position, variant.ref, variant.altString, variant);
  } else {
    vector<string> alt = split(variant.altString, ",");

    // Tri-allelic SNP.
    if (alt.size() == 2 && alt[0].size() == 1 && alt[1].size() == 1) {
      variant.isTriallelicSnp = true;
      variantMap[position].multiSnps.push_back(variant);

    // Quad-allelic SNP.
    } else if (alt.size() == 3 && alt[0].size() == 1 && alt[1].size() == 1 && alt[2].size() == 1) {
      variantMap[position].multiSnps.push_back(variant);
      variant.isQuadallelicSnp = true;

    // MNPs, insertions and deletions.
    } else {
      for (vector<string>::iterator iter = alt.begin(); iter != alt.end(); iter++) {
        // Clear the genotype string.  Otherwise, this will be kept in memory
        // for each of the alt alleles.
        variant.genotypeString = "";
        determineVariantClass(position, variant.ref, *iter, variant);
      }
    }
  }
}

// There are times when all variants for the reference sequence in the
// variant map, when all of the variants are to be processed and
// possibly written to file.
void variant::clearReferenceSequence(vcf& v, string& cRef, bool write, ostream* output) {
  while (variantMap.size() != 0) {
    if (write) {writeVariants(output);}
    variantMap.erase(vmIter);
    if (v.variantRecord.referenceSequence == cRef && v.success) {
      addVariantToStructure(v.position, v.variantRecord);
      v.success = v.getRecord(cRef);
    }
    if (variantMap.size() != 0) {vmIter = variantMap.begin();}
  }
}

// Determine the variant class from the ref and alt alleles.
void variant::determineVariantClass(int position, string& ref, string& alt, variantDescription& variant) {
  // SNP.
  if (ref.size() == 1 && (ref.size() - alt.size()) == 0) {
    variant.isBiallelicSnp = true;
    variantMap[position].biSnps.push_back(variant);

  } else {

    // Multi-base variants have the alt allele aligned to the ref allele
    // to unambiguously determine the variant type and the start
    // position.
    bool alignAlleles = false;
    if (alignAlleles) {
      int pos = position;
      string alRef, alAlt;
      string refFa = "/d2/data/references/build_37/human_reference_v37.fa";
      int start = alignAlternate(variant.referenceSequence, pos, ref, alt, alRef, alAlt, refFa); // vcf_aux.cpp
      if (pos != start) {cerr << "WARNING: Modified variant position from  " << pos << " to " << start << endl;}
    }

    // MNP.
    if (ref.size() != 1 && (ref.size() - alt.size()) == 0) {
      variant.isMnp = true;
      variantMap[position].mnps.push_back(variant);

    // Insertion.
    } else if ( alt.size() - ref.size() ) {
      variant.isInsertion = true;
      variantMap[position].indels.push_back(variant);

    // Deletion.
    } else if ( ref.size() > alt.size() ) {
      variant.isDeletion = true;
      variantMap[position].indels.push_back(variant);

    // None of the known types.
    } else {
      cerr << "Unknown variant class:" << endl;
      cerr << "Coordinates: " << variant.referenceSequence << ": " << position << endl;
      cerr << "Ref: " << ref << endl << "Alt: " << alt << endl;
      exit(1);
    }
  }
}

// Write out variants to the output.  Only process the requested variant types.
void variant::writeVariants(ostream* output) {
  // SNPs.
  if (processSnps) {
    for (variantIter = vmIter->second.biSnps.begin(); variantIter != vmIter->second.biSnps.end(); variantIter++) {
      *output << variantIter->record << endl;
    }
    for (variantIter = vmIter->second.multiSnps.begin(); variantIter != vmIter->second.multiSnps.end(); variantIter++) {
      *output << variantIter->record << endl;
    }
  }
 
  // MNPs.
  if (processMnps) {
    for (variantIter = vmIter->second.mnps.begin(); variantIter != vmIter->second.mnps.end(); variantIter++) {
      *output << variantIter->record << endl;
    }
  }
 
  // Indels.
  if (processIndels) {
    for (variantIter = vmIter->second.indels.begin(); variantIter != vmIter->second.indels.end(); variantIter++) {
      *output << variantIter->record << endl;
    }
  }
}

