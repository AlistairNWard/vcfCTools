// ******************************************************
// vcfCTools (c) 2011 Alistair Ward
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ------------------------------------------------------
// Last modified: 18 February 2011
// ------------------------------------------------------
// Tools for modifying the alleles.  This includes
// trimming the reference and alternate alleles as well
// as using the Smith-Waterman algorithm for aligning
// alleles to each other.
// ******************************************************

#define GAP '-'
#define MATCH '='
#define MISMATCH 'X'
#define INSERTION 'I'
#define DELETION 'D'
#include "modify_alleles.h"

using namespace std;
using namespace vcfCTools;

// Constructor
modifyAlleles::modifyAlleles(string& refSeq, int pos, string ref, string alt) {
  aligned              = false;
  leftAligned          = false;
  maxAllowedAlignments = 5;
  originalPosition     = pos;
  originalRef          = ref;
  originalAlt          = alt;
  referenceSequence    = refSeq;
};

// Destructor.
modifyAlleles::~modifyAlleles(void) {};

// Trim the reference and alternate alleles to produce an unambiguous
// description of the alleles.  For example, a SNP could be represented
// as CAT -> CAG as a result of other variants appearing in the same
// record.  If this is the case, the unambiguous description results when
// the alleles are trimmed to reveal T -> G.
void modifyAlleles::trim() {
  modifiedPosition = originalPosition;
  workingRef  = originalRef;
  workingAlt  = originalAlt;

  // Start at the end and work backwards.
  string::reverse_iterator rRefIter = workingRef.rbegin();
  string::reverse_iterator rAltIter = workingAlt.rbegin();
  while ( ( rRefIter != workingRef.rend() ) && ( *rRefIter == *rAltIter) && workingRef.length() > 1 && workingAlt.length() > 1) {
    workingRef.erase( --(rRefIter.base()) );
    rRefIter = workingRef.rbegin();
    workingAlt.erase( --(rAltIter.base()) );
    rAltIter = workingAlt.rbegin();
  }

  // Start at the beginning and work forwards.
  string::iterator refIter = workingRef.begin();
  string::iterator altIter = workingAlt.begin();
  while ( ( refIter != workingRef.end() ) && ( *refIter == *altIter) && workingRef.length() > 1 && workingAlt.length() > 1 ) {
    modifiedPosition++;
    refIter = workingRef.erase(refIter);
    altIter = workingAlt.erase(altIter);
  }

  // Update the modified alleles.
  modifiedRef = workingRef;
  modifiedAlt = workingAlt;
}

// Extend the reference and alternate alleles using both the reference sequence
// and a string of Z's.
void modifyAlleles::extendAlleles() {

  // Set the boolean flag aligned to false.  Only after the alignment process
  // has finished and a left aligned variant has been found is this set to 
  // true.  It is possible that multiple Smith-Waterman calls will be made as
  // it may be necessary to extend the flanking sequences to fully left-align
  // an indel.
  aligned          = false;
  numberAlignments = 0;

  // Set the variant type variable equal to that carried to the
  // routine from the variant class.
  unsigned int variantLength;

  // Determine the length of the flanks to add to the alleles.

  // SNPs.
  if (type.isBiallelicSnp || type.isTriallelicSnp || type.isQuadallelicSnp) {

  // MNPs.
  } else if (type.isMnp) {

  // Insertions/Deletions.
  } else if (type.isInsertion) {
     variantLength = workingAlt.size() - 1;
     flankLength = (variantLength < 5) ? variantLength * 5 : variantLength * 2;

  // Deletions.
  } else if (type.isDeletion) {
     variantLength = workingRef.size() - 1;
     flankLength = (variantLength < 5) ? variantLength * 5 : variantLength * 2;

  // Complex variants.
  } else if (type.isComplex) {
  }
  getFlankingReference();

  // Build the alleles to align.
  workingRef = "ZZZZZZZZZZ" + flankFront + modifiedRef + flankEnd + "ZZZZZZZZZZ";
  workingAlt = "ZZZZZZZZZZ" + flankFront + modifiedAlt + flankEnd + "ZZZZZZZZZZ";
}

// Align the alleles to each other.
void modifyAlleles::alignAlleles() {

  while (!aligned) {

    // Align the alleles to each other.
    smithWaterman();

    // Generate the CIGAR string separating matches and mismatches.
    generateCigar();

    // The variant type is known, so check that the CIGAR contains the expected
    // values.
    checkAlignment();

    // Determine if the allele is fully left aligned and realign if not.  Find
    // the genomic coordinate of the new alleles.
    processAlignment();
  }
}

// Use fastahack to find the reference flanking the ref allele.
void modifyAlleles::getFlankingReference() {

  // Add flanking sequence to the ref and alt allele in order to perform a Smith
  // Waterman alignemnt.
  int frontPos = modifiedPosition - flankLength - 1;
  int length   = modifiedRef.length() + 2 * flankLength;
  if (frontPos <= 0) {frontPos = 1;}

  FastaReference* fr = new FastaReference(fasta);
  string flank = fr->getSubSequence(referenceSequence, frontPos, length);

  flankFront = flank.substr(0, flankLength);
  flankEnd   = flank.substr(flankLength + modifiedRef.length(), flankLength);
}

// Set parameters and use the Smith-Waterman algorithm to align the reference
// to the alternate allele.
void modifyAlleles::smithWaterman() {

// Initialise Smith-Waterman parameters.
  unsigned int referencePos;
  float matchScore       = 10.0f;
  float mismatchScore    = -9.0f;
  float gapOpenPenalty   = 15.0f;
  float gapExtendPenalty = 6.66f;

  const char* reference           = workingRef.c_str();
  const char* query               = workingAlt.c_str();
  const unsigned int referenceLen = strlen(reference);
  const unsigned int queryLen     = strlen(query);

  // Call the Smith-Waterman routines.
  CSmithWatermanGotoh sw(matchScore, mismatchScore, gapOpenPenalty, gapExtendPenalty);
  sw.Align(referencePos, workingRef, workingAlt, reference, referenceLen, query, queryLen);
}

// After the alignment has taken place, generate the CIGAR string.
// Distinguish between matches and mismatches ('=' and 'X') instead
// of the standard 'M' for all.
void modifyAlleles::generateCigar() {
  bool gapRegion      = false;
  bool matchRegion    = false;
  bool mismatchRegion = false;
  numberMatch         = 0;
  numberMismatch      = 0;
  numberDeletion      = 0;
  numberInsertion     = 0;
  matchGroups         = 0;
  mismatchGroups      = 0;
  insertionGroups     = 0;
  deletionGroups      = 0;
  ostringstream oCigar;
  unsigned int m = 0, x = 0, i = 0, d = 0;

  // Strip the leading and lagging Z's.
  workingRef = workingRef.substr(10, workingRef.size() - 20);
  workingAlt = workingAlt.substr(10, workingAlt.size() - 20);

  string::iterator refIter = workingRef.begin();
  string::iterator altIter = workingAlt.begin();
  for (; refIter != workingRef.end(); refIter++) {

    // If both ref and alt alleles are bases and not gaps.
    if (*refIter != GAP && *altIter != GAP) {

      // Write out to the CIGAR string if the sequence has
      // changed from gapped to match/mismatch.
      if (gapRegion) {
        if (d != 0) {
          oCigar << d << "D";
          deletionGroups++;
        } else {
          oCigar << i << "I";
          insertionGroups++;
        }
      }
      gapRegion = false;

      // If the bases match, add one to the m value.
      if (*refIter == *altIter) {
        if (!matchRegion && x != 0) {
          oCigar << x << "X";
          mismatchGroups++;
        }
        m++;
        numberMatch++;
        x = 0;
        matchRegion    = true;
        mismatchRegion = false;

      // If the bases mismatch.
      } else {
        if (!mismatchRegion && m != 0) {
          oCigar << m << "=";
          matchGroups++;
        }
        x++;
        numberMismatch++;
        m = 0;
        matchRegion    = false;
        mismatchRegion = true;
      }

      // Reset the deletion and insertion lengths to zero since the bases matched/mismatched
      // at this position.
      d = 0;
      i = 0;

    // If the ref and alt both have a gap at this position.
    } else {
      if (!gapRegion) {

        // If not already in a gap region, flush the matches/mismatches to
        // the CIGAR.
        if (m != 0) {
          oCigar << m << '=';
          matchGroups++;
        } else if (x != 0) {
          oCigar << x << 'X';
          mismatchGroups++;
        }
      }
      gapRegion      = true;
      matchRegion    = false;
      mismatchRegion = false;
      m = 0;
      x = 0;
      if (*refIter == GAP ) {
        if (d != 0) {
          oCigar << d << 'D';
          deletionGroups++;
        }
        i++;
        numberInsertion++;
        d = 0;
      } else {
        if (i != 0) {
          oCigar << i << 'I';
          insertionGroups++;
        }
        d++;
        numberDeletion++;
        i = 0;
      }
    }
    altIter++;
  }
  if (m != 0) {
    oCigar << m << '=';
    matchGroups++;
  } else if (x != 0) {
    oCigar << d << 'X';
    mismatchGroups++;
  } else if (d != 0) {
    oCigar << d << 'D';
    deletionGroups++;
  } else if (i != 0) {
    oCigar << i << 'I';
    insertionGroups++;
  }

  cigar = oCigar.str();
}

// Check that the alignment has produced the expected results.  For example,
// if the variant is a deletion, there should only be matches and deletions
// in the CIGAR string.  There should also only be a single contiguous
// deletion.
void modifyAlleles::checkAlignment() {

  // SNPs.
  if (type.isBiallelicSnp || type.isTriallelicSnp || type.isQuadallelicSnp) {
    if ( (numberInsertion + numberDeletion) != 0 && numberMismatch == 1) {
      cerr << "ERROR: Alignment of SNP failed.  CIGAR: " << cigar << endl;
      exit(1);
    }

  // MNPs.
  } else if (type.isMnp) {
    if ( (numberInsertion + numberDeletion) != 0) {
      cerr << "ERROR: Alignment of MNP failed.  CIGAR: " << cigar << endl;
      exit(1);
    }

  // Insertions.
  } else if (type.isInsertion) {
    if ( (numberMismatch + numberDeletion) != 0 && insertionGroups == 1) {
      cerr << "ERROR: Alignment of insertion failed.  CIGAR: " << cigar << endl;
      exit(1);
    }

  // Deletions.
  } else if (type.isDeletion) {
    if ( (numberMismatch + numberInsertion) != 0 && deletionGroups == 1) {
      cerr << "ERROR: Alignment of deletion failed.  CIGAR: " << cigar << endl;
      exit(1);
    }

  // Complex events.
  } else if (type.isComplex) {
  }
}

// After aligning the alleles to each other remove the flanks and
// interrogate the resulting alignment.
void modifyAlleles::processAlignment() {
  unsigned int start;
  string startingNumber = "";

  // Find the first entry in the CIGAR string.
  string::iterator cigarIter = cigar.begin();
  for (; cigarIter != cigar.end(); cigarIter++) {

    // If the CIGAR string begins with matches, determine the left aligned
    // ref and alt alleles and the starting position.
    if (*cigarIter == MATCH) {
      start = atoi(startingNumber.c_str());

      // If start equals 1, then this means that the CIGAR string is of the
      // form 1=9D, for a deletion example.  Since there was flanking
      // sequence added to the alleles, this means that the alleles have
      // been left aligned all the way to the beginning of the new alleles.
      // It is thus possible that if the flanking sequence was longer, the
      // alleles could be further left aligned.  Thus do not mark the
      // alignment as complete and increase the size of the flanks.
      if (start == 1) {
        updateWorkingAlleles();
      } else {
        aligned = true;

        // Determine the starting coordinate of the modified alleles.
        workingPosition = modifiedPosition - flankLength + start - 1;
        if (workingPosition != modifiedPosition) {
          leftAligned = true;
          modifiedPosition = workingPosition;
        }
  
        // Trim off all but one of the matching bases at the start of the
        // alleles.
        workingRef = workingRef.substr(start - 1, workingRef.length() - start + 1);
        workingAlt = workingAlt.substr(start - 1, workingAlt.length() - start + 1);
  
        // Now trim of all the bases after the end of the alleles.  This
        // depends on the type of variant observed.  This marks the end
        // of the alignment and so put the results into the modifiedRef
        // and modifiedAlt alleles.
        if (type.isBiallelicSnp || type.isTriallelicSnp || type.isQuadallelicSnp) {
        } else if (type.isMnp) {
        } else if (type.isInsertion) {
          modifiedRef = workingRef.substr(0, 1);
          modifiedAlt = workingAlt.substr(0, numberInsertion + 1);
        } else if (type.isDeletion) {
          modifiedRef = workingRef.substr(0, numberDeletion + 1);
          modifiedAlt = workingAlt.substr(0, 1);
        } else if (type.isComplex) {
        }
      }

      break;
    } else if (*cigarIter == MISMATCH) {
      break;
    } else if (*cigarIter == INSERTION) {
      break;

    // If the variant is a deletion (any other class of variant would have failed
    // the checkAlignment if it contained a deletion) and the first entry in the
    // CIGAR is a deletion, then the alleles can be further left aligned.  Update
    // the working alleles and continue the alignment.
    } else if (*cigarIter == DELETION) {
      updateWorkingAlleles();
      break;

    } else {
      startingNumber += *cigarIter;
    }
  }
}

// If the alleles can be further left aligned, update the working
// alleles and send back to the alignment process.
void modifyAlleles::updateWorkingAlleles() {

  // Double the size of the flanking reference sequence.
  flankLength = 2 * flankLength;
  getFlankingReference();

  // Build the alleles to align.
  workingRef = "ZZZZZZZZZZ" + flankFront + modifiedRef + flankEnd + "ZZZZZZZZZZ";
  workingAlt = "ZZZZZZZZZZ" + flankFront + modifiedAlt + flankEnd + "ZZZZZZZZZZ";
}
