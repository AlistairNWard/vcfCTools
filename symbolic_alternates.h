// ******************************************************
// vcfCTools (c) 2011 Alistair Ward
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ------------------------------------------------------
// Last modified: 18 February 2011
// ------------------------------------------------------
// Tools for handling the symbolic alternate allele
// descriptions for structural variants.
// ******************************************************

#ifndef SYMBOLIC_ALTERNATES_H
#define SYMBOLIC_ALTERNATES_H

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <map>
#include <vector>

#include "header.h"
#include "variant.h"

using namespace std;
using namespace vcfCTools;

void validateAlternateAlleles(vcfHeader&, variant&);

#endif // SYMBOLIC_ALTERNATES_H
