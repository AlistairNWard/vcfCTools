// ******************************************************
// vcfCTools (c) 2011 Alistair Ward
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ------------------------------------------------------
// Last modified: 18 February 2011
// ------------------------------------------------------
// Define the bedStructure class.  This is essentially the
// equivalent of the variant class, but for bed files and
// as such is significantly simpler.
// ******************************************************

#ifndef BEDSTRUCTURE_H
#define BEDSTRUCTURE_H

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <map>
#include <vector>

#include "bed.h"
#include "structures.h"

using namespace std;

namespace vcfCTools {

class bedStructure {
  public:
    bedStructure(void);
    ~bedStructure(void);
    void addIntervalToStructure(bedRecord&);
    bool buildBedStructure(bed&);
    void generateCommonInterval(bedRecord&, bedRecord&, int, int);
    void initialiseBedMap(bed&, intFlags&);
    void resolveOverlaps(bedRecord&);
    void resolveOverlaps(bool);

  public:
    bool lastBedInterval;
    unsigned int lastBedIntervalEnd;
    unsigned int recordsInMemory;
    map<int, bedRecord> bedMap;
    map<int, bedRecord>::iterator bmIter;
    map<int, bedRecord>::iterator bmNext;
};

} // namespace vcfCTools

#endif // BEDSTRUCTURE_H
