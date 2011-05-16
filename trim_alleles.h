// ******************************************************
// vcfCTools (c) 2011 Alistair Ward
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ------------------------------------------------------
// Last modified: 18 February 2011
// ------------------------------------------------------
// Trim the reference and alternate alleles to determine
// the variant type.
// ******************************************************

#ifndef TRIM_ALLELES_H
#define TRIM_ALLELES_H

#include <cstdio>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <getopt.h>
#include <stdlib.h>
#include <map>
#include <vector>
#include <algorithm>

#include "Fasta.h"
#include "SmithWatermanGotoh.h"
#include "split.h"

using namespace std;

unsigned int trimAlleles(string&, unsigned int, string&, string&, string&, string&);

#endif
