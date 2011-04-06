// ******************************************************
// vcfCTools (c) 2011 Alistair Ward
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ------------------------------------------------------
// Last modified: 18 February 2011
// ------------------------------------------------------
// Auxillary tools for vcf.
// ******************************************************

#ifndef VCF_AUX_H
#define VCF_AUX_H

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

unsigned int alignAlternate(string&, int&, string&, string&, string&, string&, string&);
void getFlankingReference(string&, int&, string&, string&, string&, string&);
void smithWaterman(string&, int&, string&, string&, string&, string&);

#endif
