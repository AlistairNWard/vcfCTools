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
  //bool isMnp = false;

  // Start at the end and work backwards.
  string::reverse_iterator rRefIter = alRef.rbegin();
  string::reverse_iterator rAltIter = alAlt.rbegin();
  while ( ( rRefIter != alRef.rend() ) && ( *rRefIter == *rAltIter) && alRef.length() > 1 && alAlt.length() > 1) {
    alRef.erase( --(rRefIter.base()) );
    rRefIter = alRef.rbegin();
    alAlt.erase( --(rAltIter.base()) );
    rAltIter = alAlt.rbegin();
  }

  // Start at the beginning and work forwards.
  string::iterator refIter = alRef.begin();
  string::iterator altIter = alAlt.begin();
  while ( ( refIter != alRef.end() ) && ( *refIter == *altIter) && alRef.length() > 1 && alAlt.length() > 1 ) {
    start++;
    refIter = alRef.erase(refIter);
    altIter = alAlt.erase(altIter);
  }

  return start;
}
