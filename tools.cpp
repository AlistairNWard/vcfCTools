// ******************************************************
// vcfCTools (c) 2011 Alistair Ward
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ------------------------------------------------------
// Last modified: 25 February 2011
// ------------------------------------------------------
// Additional tools.
// ******************************************************

#include "tools.h"

using namespace std;

// Calculate a factorial.
unsigned int fact(unsigned int& x) {
  unsigned int y;
  unsigned int z = 1;

  for (y = 0; y < x; y++) {z = z * (y + 1);}

  return z;
}
