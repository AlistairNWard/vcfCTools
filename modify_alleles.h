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

#ifndef MODIFY_ALLELES_H
#define MODIFY_ALLELES_H

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <map>
#include <vector>

#include "Fasta.h"
#include "SmithWatermanGotoh.h"
#include "split.h"
#include "structures.h"

using namespace std;

namespace vcfCTools {

class modifyAlleles {
  public:
    modifyAlleles(string&, int, string, string);
    ~modifyAlleles(void);
    void alignAlleles();
    void checkAlignment();
    void extendAlleles();
    void generateCigar();
    void getFlankingReference();
    void leftAlign();
    void smithWaterman();
    void stepAlleles();
    void trim();
    void updateWorkingAlleles();

  public:
    bool leftAligned;
    int originalPosition;
    int modifiedPosition;
    int workingPosition;
    string fasta;
    string modifiedRef;
    string modifiedAlt;
    string originalRef;
    string originalAlt;
    string referenceSequence;
    string workingRef;
    string workingAlt;
    variantType type;

    // Variables for extending the alleles.
    unsigned int flankLength;
    unsigned int numberInsertion;
    unsigned int numberDeletion;
    unsigned int numberMatch;
    unsigned int numberMismatch;
    unsigned int matchGroups;
    unsigned int mismatchGroups;
    unsigned int insertionGroups;
    unsigned int deletionGroups;
    string cigar;
    string flank;
    string flankFront;
    string flankEnd;

    // Variables for finding repeat structure and stepping
    // alleles.
    bool hasRepetition;
    bool keepSplitting;
    unsigned int minRepeatLength;
    string deleted;
    string repeatElement;
    string sequence;
    string strippedBase;
};

} // namespace vcfCTools

#endif // MODIFY_ALLELES_H
