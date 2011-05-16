// ******************************************************
// vcfCTools (c) 2011 Alistair Ward
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ------------------------------------------------------
// Last modified: 25 February 2011
// ------------------------------------------------------
// Trim the ref and alt alleles to find the correct
// variant type.
// ******************************************************

#include "trim_alleles.h"

unsigned int trimAlleles(string& referenceSequence, unsigned int position, string& ref, string& alt, string& alRef, string& alAlt) {
  unsigned int start = position;
  string fRef;
  string fAlt;
  alRef = ref;
  alAlt = alt;
  bool isMnp = false;

  // Determine if the variant is an MNP or an indel.
  if (ref.size() == alt.size()) {isMnp = true;}

  // Start at the beginning and work forwards.
  string::iterator refIter = alRef.begin();
  string::iterator altIter = alAlt.begin();
  fRef = *refIter;
  fAlt = *altIter;
  while ( ( refIter != alRef.end() ) && ( *refIter == *altIter) ) {
    start++;
    fRef = *refIter;
    fAlt = *altIter;
    refIter = alRef.erase(refIter);
    altIter = alAlt.erase(altIter);
  }

  // Start at the end and work backwards.
  string::reverse_iterator rRefIter = alRef.rbegin();
  string::reverse_iterator rAltIter = alAlt.rbegin();
  while ( ( rRefIter != alRef.rend() ) && ( *rRefIter == *rAltIter) ) {
    alRef.erase( --(rRefIter.base()) );
    rRefIter = alRef.rbegin();
    alAlt.erase( --(rAltIter.base()) );
    rAltIter = alAlt.rbegin();
  }

  // If the variant is an indel, replace the first base at the
  // beginning of both the ref and the alt allele.
  if (!isMnp) {
    alRef = fRef + alRef;
    alAlt = fAlt + alAlt;
    start--;
  }

  return start;
}
