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
  string anchor;

  anchor.assign(50, 'Z');
  flankLength   = 50;
  getFlankingReference();

  // Build the alleles to align.
  workingRef = anchor + flankFront + modifiedRef + flankEnd + anchor;
  workingAlt = anchor + flankFront + modifiedAlt + flankEnd + anchor;
}

// Align the alleles to each other.
void modifyAlleles::alignAlleles() {

  // Align the alleles to each other.
  smithWaterman();

  // Generate the CIGAR string separating matches and mismatches.
  generateCigar();
  cerr << modifiedRef << endl << modifiedAlt << endl << endl;
  cerr << workingRef << endl << workingAlt << endl << endl;
  cerr << cigar << endl;
  exit(0);

  // The variant type is known, so check that the CIGAR contains the expected
  // values.
  //checkAlignment();

  // Determine if the allele is fully left aligned and realign if not.  Find
  // the genomic coordinate of the new alleles.
  //processAlignment();
  //leftAlign();
}

// Use fastahack to find the reference flanking the ref allele.
void modifyAlleles::getFlankingReference() {

  // Add flanking sequence to the ref and alt allele in order to perform a Smith
  // Waterman alignemnt.
  int frontPos = modifiedPosition - flankLength - 1;
  int length   = modifiedRef.length() + 2 * flankLength;
  if (frontPos <= 0) {frontPos = 1;}

  FastaReference* fr = new FastaReference(fasta);
  flank = fr->getSubSequence(referenceSequence, frontPos, length);
  delete fr;

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
  workingRef = workingRef.substr(15, workingRef.size() - 30);
  workingAlt = workingAlt.substr(15, workingAlt.size() - 30);

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
    if ( (numberInsertion + numberDeletion) != 0 || numberMismatch != 1) {
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
    if ( (numberMismatch + numberDeletion) != 0 || insertionGroups != 1) {
      cerr << "ERROR: Alignment of insertion failed.  CIGAR: " << cigar << endl;
      exit(1);
    }

  // Deletions.
  } else if (type.isDeletion) {
    if ( (numberMismatch + numberInsertion) != 0 || deletionGroups != 1) {
      cerr << "ERROR: Alignment of deletion failed at " << referenceSequence << ":";
      cerr << originalPosition << ".  CIGAR: " << cigar << endl;
      exit(1);
    }

  // Complex events.
  } else if (type.isComplex) {
  }
}

// Left align the allele.
void modifyAlleles::leftAlign() {
  int altpos = 0;
  int refpos = 0;
  int len;
  string slen;
  vector<pair<int, char> > cigarData;

  cerr << "START: " << endl << workingRef << endl << workingAlt << endl << cigar << endl;

  for (string::iterator c = cigar.begin(); c != cigar.end(); ++c) {
  cerr << *c << endl;
  switch (*c) {
    case 'I':
      len = atoi(slen.c_str());
      slen.clear();
      cigarData.push_back(make_pair(len, *c));
      cerr << "CASE I: " << workingAlt.substr(altpos, len) << " " <<  refpos + workingPosition << endl;
      //variants.push_back(VariantAllele("", alternateQuery.substr(altpos, len), refpos - paddingLen + position));
      altpos += len;
      break;
    case 'D':
      len = atoi(slen.c_str());
      slen.clear();
      cigarData.push_back(make_pair(len, *c));
      cerr << "CASE D: " << workingRef.substr(refpos, len) << " " << refpos + workingPosition << endl;
      //variants.push_back(VariantAllele(reference.substr(refpos, len), "", refpos - paddingLen + position));
      refpos += len;
      break;
    case '=':
      cerr << "MATCH" << endl;
      {
        len = atoi(slen.c_str());
        cerr << "LEN=" << len << endl;
        slen.clear();
        cigarData.push_back(make_pair(len, *c));
        //string refmatch = reference.substr(refpos, len);
        //string altmatch = alternateQuery.substr(altpos, len);
        string refmatch = workingRef.substr(refpos, len);
        string altmatch = workingAlt.substr(altpos, len);
        cerr << refmatch << " " << altmatch << endl;
        bool inmismatch = false;
        int mismatchStart = 0;
        for (int i = 0; i < refmatch.size(); ++i) {
          cerr << "COMPARE: " << refmatch.at(i) << " " << altmatch.at(i) << endl;
          if (refmatch.at(i) == altmatch.at(i)) {
            cerr << "EQUAL " << inmismatch << endl;
            if (inmismatch) {
              //variants.push_back(VariantAllele(
                //refmatch.substr(mismatchStart, i - mismatchStart),
                //altmatch.substr(mismatchStart, i - mismatchStart),
                //mismatchStart - paddingLen + position));
              cerr << "CASE M: " << refmatch.substr(mismatchStart, i - mismatchStart) << " ";
              cerr << altmatch.substr(mismatchStart, i - mismatchStart) << " ";
              cerr << mismatchStart + workingPosition << endl;
            }
            inmismatch = false;
          } else {
            cerr << "NOT EQUAL" << endl;
            if (!inmismatch) {
              mismatchStart = i;
              inmismatch = true;
            }
          }
          ++refpos;
          ++altpos;
        }
      }
      break;
    case 'S':
      len = atoi(slen.c_str());
      slen.clear();
      cigarData.push_back(make_pair(len, *c));
      refpos += len;
      altpos += len;
      break;
    default:
      len = 0;
      slen += *c;
      break;
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
  workingRef = "ZZZZZZZZZZZZZZZ" + flankFront + modifiedRef + flankEnd + "ZZZZZZZZZZZZZZZ";
  workingAlt = "ZZZZZZZZZZZZZZZ" + flankFront + modifiedAlt + flankEnd + "ZZZZZZZZZZZZZZZ";
}

// Attempt to place the indel further left than reported in the vcf file.
// Use fastahack to get flanking reference sequence.
void modifyAlleles::stepAlleles() {
  bool leftAligned = false;
  bool reset;
  unsigned int count;
  unsigned int refPosition;
  unsigned int indelPosition;
  string anchor;
  string laggingBase;

  // Get the flanking reference sequence.  The variable sequence is
  // populated with the inserted/deleted bases.
  FastaReference* fr = new FastaReference(fasta);
  if (type.isInsertion) {
    flankLength = 30 * modifiedAlt.length();
    sequence    = modifiedAlt.substr(1, modifiedAlt.length());
    flank       = fr->getSubSequence(referenceSequence, originalPosition - 1 - flankLength, flankLength);
    flank      += modifiedAlt;
    laggingBase = fr->getSubSequence(referenceSequence, originalPosition, 1);
  } else if (type.isDeletion) {
    flankLength = 30 * modifiedRef.length();
    sequence    = modifiedRef.substr(1, modifiedRef.length());
    flank       = fr->getSubSequence(referenceSequence, originalPosition - 1 - flankLength, flankLength);
    flank      += modifiedRef;
    laggingBase = fr->getSubSequence(referenceSequence, originalPosition + sequence.length(), 1);
  }

  // Try stepping the inserted/deleted alleles backwards through the
  // reference sequence.  The first base is an anchor base and may be
  // permitted to change.  For example, consider the deletion of an
  // AC in a dinucleotide repeat.  If the alleles are represented as
  // CAC -> C, the starting C in the ref allele may be subject to
  // change.  If the context is ACGACAC and the coordinate referes to
  // deletion of the second AC, this can be replaced by GAC -> G and
  // thus the starting C is not preserved.
  count           = 0;
  refPosition     = flank.length() - sequence.length() - 1;
  workingPosition = modifiedPosition;
  while (count < sequence.length()) {
    indelPosition = 0;
    while (sequence[indelPosition] == flank[refPosition + indelPosition] && indelPosition < sequence.length()) {

      // If the whole sequence matches the reference then the sequence
      // can be moved to this position.
      if (indelPosition == sequence.length() - 1) {
        workingPosition  = workingPosition - count - 1;
        count            = 0;
        anchor           = flank[refPosition - 1];

        // Ensure that replacing the reported allele representation with
        // the left-aligned one results in the same alleles.
        string oldAllele, newAllele;
        if (type.isInsertion) {
          oldAllele = flank.substr(refPosition - 1, flankLength - refPosition + 2) + sequence + laggingBase;
          newAllele = anchor[0] + flank.substr(refPosition - 1 + sequence.length() - 1, flankLength - refPosition + 1) + sequence + laggingBase;
        } else if (type.isDeletion) {
          oldAllele = flank.substr(refPosition - 1, flankLength - refPosition + 2) + laggingBase;
          newAllele = anchor[0] + flank.substr(refPosition - 1 + sequence.length() - 1, flankLength - refPosition + 1) + laggingBase;
        }
        if (oldAllele == newAllele) {
          leftAligned = true;
          reset       = true;
        }
      }
      indelPosition++;
    }
    refPosition--;
    if (reset) {reset = false;}
    else {count++;}
  }

  if (leftAligned) {

    // Change the first base in the reference and alternate alleles to
    // the anchor base.  This is the leftmost base that both alleles
    // share.
    modifiedRef[0]   = anchor[0];
    modifiedAlt[0]   = anchor[0];
    modifiedPosition = workingPosition;
  }

  delete fr;
}
