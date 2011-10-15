// ******************************************************
// vcfCTools (c) 2011 Alistair Ward
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ------------------------------------------------------
// Last modified: 18 February 2011
// ------------------------------------------------------
// Tools for handling symbolic alternate allele
// descriptions for structural variants.
// ******************************************************

#include "symbolic_alternates.h"

using namespace std;
using namespace vcfCTools;

// Check that the alternate alleles conform to the spec.
void validateAlternateAlleles(vcfHeader&, variant& var) {

  // Loop over all records at this locus.
  vector<string>::iterator aIter         = var.ovIter->alts.begin();
  vector<variantType>::iterator typeIter = var.ovIter->type.begin();
  for (; aIter != var.ovIter->alts.end(); aIter++) {
    if (typeIter->isSv) {
      cerr << *aIter << " " << typeIter->isSv << endl;
    } else if (typeIter->isRearrangement) {
      cerr << *aIter << " " << typeIter->isSv << endl;
    }
    typeIter++;
  }
}
