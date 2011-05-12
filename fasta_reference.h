// ******************************************************
// vcfCTools (c) 2011 Alistair Ward
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ------------------------------------------------------
// Last modified: 18 February 2011
// ------------------------------------------------------
// Tools to extract information from a fasta reference.
// ******************************************************

#ifndef FASTA_REFERENCE_H
#define FASTA_REFERENCE_H

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

void defineFastaReference(string);
void getFlankingReference(string, int, string, string, unsigned int, string&, string&, string);

#endif
