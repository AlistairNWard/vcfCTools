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

using namespace std;

namespace vcfCTools {

class bedStructure {
  public:
    bedStructure(void);
    ~bedStructure(void);
    bool buildBedStructure(bed&);
    void addIntervalToStructure(bedRecord&);
    void resolveOverlaps(bedRecord&);
    void resolveOverlaps(bool);
    void generateCommonInterval(bedRecord&, bedRecord&, int, int);

  public:
    unsigned int recordsInMemory;
    map<int, bedRecord> bedMap;
    map<int, bedRecord>::iterator bmIter;
};

} // namespace vcfCTools

#endif // BEDSTRUCTURE_H
