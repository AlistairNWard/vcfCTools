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
#include "output.h"
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
    void checkReferenceSequences(variant&, variant&);
    void intersectVcf(vcf&, variant&, vcf&, variant&, output&);
    void intersectVcfBed(vcf&, variant&, bed&, bedStructure&, output&);
    void setBooleanFlags(bool, bool, bool, bool, bool, bool);

  public:
    intFlags flags;
    map<string, map<int, unsigned int> > distanceDist;
};

} // namespace vcfCTools

#endif
