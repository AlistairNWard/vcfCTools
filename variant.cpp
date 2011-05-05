// ******************************************************
// vcfCTools (c) 2011 Alistair Ward
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ------------------------------------------------------
// Last modified: 18 February 2011
// ------------------------------------------------------
// vcfClass describes the vcf class and all operations.
// ******************************************************

#include "info.h"
#include "tools.h"
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

  variantMap[position].hasBiallelicSnp = false;
  variantMap[position].hasMultiallelicSnp = false;
  variantMap[position].hasMnp = false;
  variantMap[position].hasIndel = false;

// If this is the first variant at this position, create a new entry in the
// information structure.
  variantMap[position].referenceSequence = variant.referenceSequence;

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
      variantMap[position].hasMultiallelicSnp = true;

    // Quad-allelic SNP.
    } else if (alt.size() == 3 && alt[0].size() == 1 && alt[1].size() == 1 && alt[2].size() == 1) {
      variantMap[position].multiSnps.push_back(variant);
      variantMap[position].hasMultiallelicSnp = true;
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
void variant::clearReferenceSequence(vcf& v, string cRef, bool write, ostream* output) {
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
    variant.altString = alt;
    variantMap[position].biSnps.push_back(variant);
    variantMap[position].hasBiallelicSnp = true;

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
      variant.altString = alt;
      variantMap[position].mnps.push_back(variant);
      variantMap[position].hasMnp = true;

    // Insertion.
    } else if ( alt.size() > ref.size() ) {
      variant.isInsertion = true;
      variant.altString = alt;
      variantMap[position].indels.push_back(variant);
      variantMap[position].hasIndel = true;

    // Deletion.
    } else if ( ref.size() > alt.size() ) {
      variant.isDeletion = true;
      variant.altString = alt;
      variantMap[position].indels.push_back(variant);
      variantMap[position].hasIndel = true;

    // None of the known types.
    } else {
      cerr << "Unknown variant class:" << endl;
      cerr << "Coordinates: " << variant.referenceSequence << ": " << position << endl;
      cerr << "Ref: " << ref << endl << "Alt: " << alt << endl;
      exit(1);
    }
  }
}

// Annotate the variants at this locus with the contents of the vcf file.
void variant::annotateRecordVcf(variantsAtLocus& v, bool isDbsnp) {
  variantInfo annInfo;
  variantInfo vcfInfo;
  vector<variantDescription>::iterator iter;

// If the vcf file used for annotation is a dbSNP vcf file, only compare the
// relevant variant classes.  Also, parse the info string to check that the
// variant class as determined by vcfCTools agrees with that listed in the
// variant class (VC) info field.
  if (isDbsnp) {

    // Annotate SNPs.
    if (vmIter->second.biSnps.size() != 0) {
      for (iter = v.biSnps.begin(); iter != v.biSnps.end(); iter++) {
  
        // Check that the VC class lists this entry as a SNP.
        annInfo.processInfoFields(iter->info);
        annInfo.getInfo(string("VC"), v.referenceSequence, vmIter->first);
        if (annInfo.values[0] != "SNP") {
          cerr << "WARNING: dbSNP entry listed as " << annInfo.values[0] << ", expected SNP.";
          cerr << " Coordinate (" << vmIter->second.referenceSequence << ":" << vmIter->first << ")" << endl;
        }
  
        // Annotate the SNPs with the rsid value.
        for (variantIter = vmIter->second.biSnps.begin(); variantIter != vmIter->second.biSnps.end(); variantIter++) {
          variantIter->rsid = iter->rsid;
          vcfInfo.processInfoFields(variantIter->info);
          if (vcfInfo.infoTags.count("dbSNP") == 0) {variantIter->info += ";dbSNP";}
          buildRecord(vmIter->first, *variantIter); // tools.cpp
        }
      }
    }

    // Annotate multiallelic SNPs.
    if (vmIter->second.multiSnps.size() != 0) {
      for (iter = v.multiSnps.begin(); iter != v.multiSnps.end(); iter++) {
  
        // Check that the VC class lists this entry as a SNP.
        annInfo.processInfoFields(iter->info);
        annInfo.getInfo(string("VC"), v.referenceSequence, vmIter->first);
        if (annInfo.values[0] != "SNP") {
          cerr << "WARNING: dbSNP entry listed as " << annInfo.values[0] << ", expected SNP.";
          cerr << " Coordinate (" << vmIter->second.referenceSequence << ":" << vmIter->first << ")" << endl;
        }
  
        // Annotate the SNPs with the rsid value.
        for (variantIter = vmIter->second.multiSnps.begin(); variantIter != vmIter->second.multiSnps.end(); variantIter++) {
          variantIter->rsid = iter->rsid;
          vcfInfo.processInfoFields(variantIter->info);
          if (vcfInfo.infoTags.count("dbSNP") == 0) {variantIter->info += ";dbSNP";}
          buildRecord(vmIter->first, *variantIter); // tools.cpp
        }
      }
    }

    // Annotate MNPs.
    if (vmIter->second.mnps.size() != 0) {
      for (iter = v.mnps.begin(); iter != v.mnps.end(); iter++) {
  
        // Check that the VC class lists this entry as an MNP.
        annInfo.processInfoFields(iter->info);
        annInfo.getInfo(string("VC"), v.referenceSequence, vmIter->first);
        if (annInfo.values[0] != "MULTI-BASE") {
          cerr << "WARNING: dbSNP entry listed as " << annInfo.values[0] << ", expected MNP (MULTI-BASE).";
          cerr << " Coordinate (" << vmIter->second.referenceSequence << ":" << vmIter->first << ")" << endl;
        }
  
        // Annotate the SNPs with the rsid value.
        for (variantIter = vmIter->second.mnps.begin(); variantIter != vmIter->second.mnps.end(); variantIter++) {
          vcfInfo.processInfoFields(variantIter->info);
          if (vcfInfo.infoTags.count("dbSNP") == 0) {variantIter->info += ";dbSNP";}
          variantIter->rsid = iter->rsid;
          buildRecord(vmIter->first, *variantIter); // tools.cpp
        }
      }
    }

    // Annotate indels.
    if (vmIter->second.indels.size() != 0) {
      for (iter = v.indels.begin(); iter != v.indels.end(); iter++) {
  
        // Check that the VC class lists this entry as an indel.
        annInfo.processInfoFields(iter->info);
        annInfo.getInfo(string("VC"), v.referenceSequence, vmIter->first);
        if (annInfo.values[0] != "INDEL") {
          cerr << "WARNING: dbSNP entry listed as " << annInfo.values[0] << ", expected INDEL.";
          cerr << " Coordinate (" << vmIter->second.referenceSequence << ":" << vmIter->first << ")" << endl;
        }
  
        // Annotate the SNPs with the rsid value.
        for (variantIter = vmIter->second.indels.begin(); variantIter != vmIter->second.indels.end(); variantIter++) {
          vcfInfo.processInfoFields(variantIter->info);
          if (vcfInfo.infoTags.count("dbSNP") == 0) {variantIter->info += ";dbSNP";}
          variantIter->rsid = iter->rsid;
          buildRecord(vmIter->first, *variantIter); // tools.cpp
        }
      }
    }

// If a vcf file is provided for performing annotations that is not a dbSNP
// vcf file, add the info string from this file to the vcf file being
// annotated.
  } else {

    // Annotate SNPs.
    for (iter = v.biSnps.begin(); iter != v.biSnps.end(); iter++) {

      // Add the filter string from the annotation file to the info string of
      // the vcf file.  First check that this does not already exist to avoid
      // multiple instances of the same flag.
      for (variantIter = vmIter->second.biSnps.begin(); variantIter != vmIter->second.biSnps.end(); variantIter++) {
        vcfInfo.processInfoFields(variantIter->info);
        if (vcfInfo.infoTags.count(iter->filters) == 0) {variantIter->info += ";" + iter->filters;}
        buildRecord(vmIter->first, *variantIter);
      }
    }

    // Annotate multiallelic SNPs.
    for (iter = v.multiSnps.begin(); iter != v.multiSnps.end(); iter++) {

      // Add the filter string from the annotation file to the info string of
      // the vcf file.  First check that this does not already exist to avoid
      // multiple instances of the same flag.
      for (variantIter = vmIter->second.multiSnps.begin(); variantIter != vmIter->second.multiSnps.end(); variantIter++) {
        vcfInfo.processInfoFields(variantIter->info);
        if (vcfInfo.infoTags.count(iter->filters) == 0) {variantIter->info += ";" + iter->filters;}
        buildRecord(vmIter->first, *variantIter);
      }
    }

    // Annotate MNPs.
    for (iter = v.mnps.begin(); iter != v.mnps.end(); iter++) {

      // Add the filter string from the annotation file to the info string of
      // the vcf file.  First check that this does not already exist to avoid
      // multiple instances of the same flag.
      for (variantIter = vmIter->second.mnps.begin(); variantIter != vmIter->second.mnps.end(); variantIter++) {
        vcfInfo.processInfoFields(variantIter->info);
        if (vcfInfo.infoTags.count(iter->filters) == 0) {variantIter->info += ";" + iter->filters;}
        buildRecord(vmIter->first, *variantIter);
      }
    }

    // Annotate indels.
    for (iter = v.indels.begin(); iter != v.indels.end(); iter++) {

      // Add the filter string from the annotation file to the info string of
      // the vcf file.  First check that this does not already exist to avoid
      // multiple instances of the same flag.
      for (variantIter = vmIter->second.indels.begin(); variantIter != vmIter->second.indels.end(); variantIter++) {
        vcfInfo.processInfoFields(variantIter->info);
        if (vcfInfo.infoTags.count(iter->filters) == 0) {variantIter->info += ";" + iter->filters;}
        buildRecord(vmIter->first, *variantIter);
      }
    }
  }
}

// Annotate the variants at this locus with the contents of the bed file.
void variant::annotateRecordBed(bedRecord& b) {
  // SNPs.
  if (processSnps) {
    for (variantIter = vmIter->second.biSnps.begin(); variantIter != vmIter->second.biSnps.end(); variantIter++) {
      if (variantIter->info != "" && variantIter->info != ".") {variantIter->info += ";" + b.info;}
      else {variantIter->info = b.info;}
      buildRecord(vmIter->first, *variantIter); // tools.cpp
    }
    for (variantIter = vmIter->second.multiSnps.begin(); variantIter != vmIter->second.multiSnps.end(); variantIter++) {
      if (variantIter->info != "" && variantIter->info != ".") {variantIter->info += ";" + b.info;}
      else {variantIter->info = b.info;}
      buildRecord(vmIter->first, *variantIter); // tools.cpp
    }
  }
 
  // MNPs.
  if (processMnps) {
    for (variantIter = vmIter->second.mnps.begin(); variantIter != vmIter->second.mnps.end(); variantIter++) {
      if (variantIter->info != "" && variantIter->info != ".") {variantIter->info += ";" + b.info;}
      else {variantIter->info = b.info;}
      buildRecord(vmIter->first, *variantIter); // tools.cpp
    }
  }
 
  // Indels.
  if (processIndels) {
    for (variantIter = vmIter->second.indels.begin(); variantIter != vmIter->second.indels.end(); variantIter++) {
      if (variantIter->info != "" && variantIter->info != ".") {variantIter->info += ";" + b.info;}
      else {variantIter->info = b.info;}
      buildRecord(vmIter->first, *variantIter); // tools.cpp
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

