// ******************************************************
// vcfCTools (c) 2011 Alistair Ward
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ------------------------------------------------------
// Last modified: 25 February 2011
// ------------------------------------------------------
// Auxillary tools for vcf.
// ******************************************************

//#include "fasta_reference.h"
#include "vcf_aux.h"

using namespace std;

// Align the alt to the ref allele using a Smith-Waterman algorithm
// and find the start position of the reference
unsigned int alignAlternate(string referenceSequence, int position, string& ref, string& alt, string& alRef, string& alAlt, string refFa) {
  unsigned int start;
  unsigned int flankLength = 0;
  bool swSuccess = false;
  unsigned int repetitions = 0;
  string::iterator refIter, altIter;
  string::reverse_iterator refRevIter;
  string flankFront = "", flankEnd = "", swRef = "", swAlt = "";

// Check that the flanking sequences of the ref and the alt still match.  If not
// the Smith-Waterman algorithm generated an incorrect alignment.  Allow an increase
// in the size of the flanking sequences, but only this to be increased a specified
// number of times before terminating with an error.
  while (!swSuccess) {
    getFlankingReference(referenceSequence, position, ref, alt, flankLength, flankFront, flankEnd, refFa);
    swRef = flankFront + ref + flankEnd;
    swAlt = flankFront + alt + flankEnd;
    smithWaterman(referenceSequence, position, swRef, swAlt, alRef, alAlt);
    flankLength = flankFront.length();

    //if ( alRef.substr(0, flankLength) != alAlt.substr(0,flankLength) ||
    //     (alRef.substr(alRef.length() - flankLength, flankLength)) != (alAlt.substr(alAlt.length() - flankLength, flankLength)) ||
    //     alRef.substr(flankLength, 1) == "-" || alAlt.substr(flankLength, 1) == "-") {
    if ( alRef.substr(0, flankLength) != flankFront ||
         alAlt.substr(0,flankLength)  != flankFront ||
         (alRef.substr(alRef.length() - flankLength, flankLength)) != flankEnd || 
         (alAlt.substr(alAlt.length() - flankLength, flankLength)) != flankEnd ||
         alRef.substr(flankLength, 1) == "-" || alAlt.substr(flankLength, 1) == "-") {
      repetitions++;
      if (repetitions > 2) {
        cerr << "Difference in flanking sequences after performing SW alignment." << endl;
        cerr << "Either the alignment is incorrect or the vcf record can be" << endl;
        cerr << "further left aligned.  Please check position:" << endl;
        cerr << referenceSequence << ":" << position << endl;
        cerr << endl;
        cerr << "Ref and alt alleles with flanking reference sequence:" << endl;
        cerr << endl;
        string flank = string(flankLength, '_');
        flank.insert(flankLength, alRef.length() - 2 * flankLength, ' ');
        flank.insert(flankLength + alRef.length() - 2 * flankLength, flankLength, '_');
        cerr << "Flank: " << flank << endl;
        cerr << "Ref:   " << alRef << endl;
        cerr << "Alt:   " << alAlt << endl;
        exit(1);
      }
    } else {
      swSuccess = true;
      alRef = alRef.substr(flankLength, alRef.length() - 2 * flankLength);
      alAlt = alAlt.substr(flankLength, alAlt.length() - 2 * flankLength);

      //Find the left-most unambigous start location for the reference allele..
      start = position;
      altIter = alAlt.begin();
      for (refIter = alRef.begin(); refIter != alRef.end(); refIter++) {
        if (*refIter != *altIter) {
          break;
        }
        altIter++;
        start++;
      }
      if ( (start - position) != 0) {
        start--;
        ref.erase(0, start - position);
        alt.erase(0, start - position);
      }
      //alRef.erase(0, start - position);
      //alAlt.erase(0, start - position);
    }
  }

  return start;
}

// Use fastahack to find the reference flanking the ref allele.
void getFlankingReference(string referenceSequence, int position, string ref, string alt, unsigned int flankLength, string& flankFront, string& flankEnd, string refFa) {
  int frontPos, length;

  // Add flanking sequence to the ref and alt allele in order to perform a Smith
  // Waterman alignemnt.
  if (ref.length() > 100) {
    cerr << "reference sequence in excess of 100bp.  Not yet handled." << endl;
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

// Align two sequences using a Smith-Waterman alignment.
void smithWaterman(string referenceSequence, int position, string ref, string alt, string& alRef, string& alAlt) {

// Initialise Smith-Waterman parameters.
  unsigned int referencePos;
  float matchScore = 10.0f;
  float mismatchScore = -9.0f;
  float gapOpenPenalty = 15.0f;
  float gapExtendPenalty = 6.66f;

  //char* reference = new char[ref.size() + 1];
  //copy(ref.begin(), ref.end(), reference);
  //reference[ref.length()] = '\0';
  const char* reference = ref.c_str();
  const char* query = alt.c_str();
  //char* query = new char[alt.size() + 1];
  //copy(alt.begin(), alt.end(), query);
  //query[alt.length()] = '\0';

  const unsigned int referenceLen = strlen(reference);
  const unsigned int queryLen     = strlen(query);
  alRef = "", alAlt = "";
  CSmithWatermanGotoh sw(matchScore, mismatchScore, gapOpenPenalty, gapExtendPenalty);
  sw.Align(referencePos, alRef, alAlt, reference, referenceLen, query, queryLen);
}
