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
#include "vcf.h"

using namespace std;

namespace vcfCTools {

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
    void determineVariantsToProcess(bool, bool, bool);
    bool buildVariantStructure(vcf&);
    void addVariantToStructure(int, variantDescription&);
    void clearReferenceSequence(vcf&, string, bool, ostream*);
    void determineVariantClass(int, string, string, variantDescription&);
    void annotateRecordVcf(variantsAtLocus&, bool);
    void annotateRecordBed(bedRecord&);
    void writeVariants(ostream*);

  public:
    unsigned int recordsInMemory;
    string referenceSequence;
    map<int, variantsAtLocus> variantMap;
    map<int, variantsAtLocus>::iterator vmIter;
    vector<variantDescription>::iterator variantIter;

    // Boolean flags.
    bool processSnps;
    bool processMnps;
    bool processIndels;
};

} // namespace vcfCTools

#endif // VARIANT_H
