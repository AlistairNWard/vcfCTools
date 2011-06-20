// ******************************************************
// vcfCTools (c) 2011 Alistair Ward
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ------------------------------------------------------
// Last modified: 18 February 2011
// ------------------------------------------------------
// Define the variant class.
// ******************************************************

#ifndef VARIANT_H
#define VARIANT_H

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <map>
#include <vector>

#include "bed.h"
#include "info.h"
#include "trim_alleles.h"
#include "vcf.h"

using namespace std;

namespace vcfCTools {

// Create astructure to hold all of the flags required to determine the intersection
// operations to be performed.
struct intFlags {
  bool sitesOnly;
  bool annotate;
  bool findCommon;
  bool findUnion;
  bool findUnique;
};

// At each locus, many different variants can exist.  The
// variantsAtLocus structure holds all of the variants
// present at this locus.
struct variantsAtLocus {
  string referenceSequence;
  vector<variantDescription> biSnps;
  vector<variantDescription> multiSnps;
  vector<variantDescription> mnps;
  vector<variantDescription> indels;

  // Boolean flags.
  bool hasBiallelicSnp;
  bool hasMultiallelicSnp;
  bool hasMnp;
  bool hasIndel;
};

class variant {
  public:
    variant(void);
    ~variant(void);
    void determineVariantsToProcess(bool, bool, bool, bool);
    bool buildVariantStructure(vcf&);
    void addVariantToStructure(int, variantDescription&, bool);
    void clearReferenceSequence(vcf&, variant&, intFlags, string, bool, ostream*);
    void determineVariantClass(int, string, string, variantDescription&, bool);
    void annotateRecordVcf(variantsAtLocus&, bool);
    void annotateRecordBed(bedRecord&);
    void compareVariantsSameLocus(variant&, intFlags, string, ostream*);
    void compareVariantsDifferentLocus(variant&, intFlags, bool, ostream*);
    vector<string> extractGenotypeField(string);
    void writeVariants(int, variantsAtLocus&, ostream*);

  public:
    unsigned int recordsInMemory;
    string referenceSequence;
    map<int, variantsAtLocus> variantMap;
    map<int, variantsAtLocus>::iterator vmIter;
    vector<variantDescription>::iterator variantIter;
    map<string, headerInfoStruct> headerInfoFields;

    // Boolean flags.
    bool processSnps;
    bool processMnps;
    bool processIndels;
    bool splitMnps;
};

} // namespace vcfCTools

#endif // VARIANT_H
