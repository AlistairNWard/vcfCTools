// ******************************************************
// vcfCTools (c) 2011 Alistair Ward
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ------------------------------------------------------
// Last modified: 25 February 2011
// ------------------------------------------------------
// Tools to deal with parsing a fasta reference.
// ******************************************************

#include "fastaReference.h"

using namespace std;

// Define a fasta reference.
void defineFastaReference(string faRef) {
  FastaReference* fr = new FastaReference(faRef);
}

// Use fastahack to find the reference flanking the ref allele.
void getFlankingReference(string referenceSequence, int position, string ref, string alt, unsigned int flankLength, string& flankFront, string& flankEnd, string refFa) {
  int frontPos, length;

  // Add flanking sequence to the ref and alt allele in order to perform a Smith
  // Waterman alignemnt.
  if (ref.length() > 100) {
    cout << "reference sequence in excess of 100bp.  Not yet handled." << endl;
    exit(1);
  } else {
    if (flankLength == 0) {
      int maxLength = (ref.length() > alt.length()) ? ref.length() : alt.length();
      flankLength = (ref.length() < 10) ? 20 : 2 * maxLength;
    } else {
      flankLength += flankLength;
    }
    frontPos = position - flankLength - 1;
    length = ref.length() + 2 * flankLength;
    if (frontPos <= 0) {frontPos = 1;}
  }

  FastaReference* fr = new FastaReference(refFa);
  string flank = fr->getSubSequence(referenceSequence, frontPos, length);

  flankFront = flank.substr(0, flankLength);
  flankEnd   = flank.substr(flankLength + ref.length(), flankLength);
}
