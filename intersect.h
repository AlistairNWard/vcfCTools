// ******************************************************
// vcfCTools (c) 2011 Alistair Ward
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ------------------------------------------------------
// Last modified: 18 February 2011
// ------------------------------------------------------
// Perform file intersections
// ******************************************************

#ifndef INTERSECT_H
#define INTERSECT_H

#include <cstdio>
#include <iostream>
#include <string>
#include <getopt.h>
#include <stdlib.h>

#include "bed.h"
#include "bedStructure.h"
#include "header.h"
#include "output.h"
#include "structures.h"
#include "tools.h"
#include "variant.h"
#include "vcf.h"
#include "vcfCTools_tool.h"

using namespace std;
//using namespace vcfCTools;

namespace vcfCTools {

class intersect {
  public:
    intersect(void);
    ~intersect(void);
    void beyondInterval();
    void checkReferenceSequences(variant&, variant&);
    void intersectVcf(vcfHeader&, vcfHeader&, vcf&, variant&, vcf&, variant&, output&);
    void intersectVcfBed(vcfHeader&, vcf&, variant&, bed&, bedStructure&, output&);
    void iterateBedFile(bed&, bedStructure&);
    void iterateVcfFile(vcfHeader&, vcf&, variant&, output&);
    void nextReferenceSequence(vcf&, variant&, bed&, bedStructure&);
    bool priorToInterval(intFlags&);
    void setBooleanFlags(bool, bool, bool, bool, bool, bool);
    bool withinInterval(intFlags&);

  public:
    bool iterateBed;
    bool iterateVcf;
    string currentReferenceSequence;
    intFlags flags;
    map<string, map<int, unsigned int> > distanceDist;
};

} // namespace vcfCTools

#endif
