// ******************************************************
// vcfCTools (c) 2011 Alistair Ward
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ------------------------------------------------------
// Last modified: 18 February 2011
// ------------------------------------------------------
// String split method.
// Thanks to Evan Teran, http://stackoverflow.com/
// questions/236129/how-to-split-a-string/236803#236803
// ******************************************************

#ifndef SPLIT_H
#define SPLIT_H

#include <string>
#include <vector>
#include <sstream>
#include <string.h>

// split a string on a single delimiter character (delim)
std::vector<std::string>& split(const std::string &s, char delim, std::vector<std::string> &elems);
std::vector<std::string>  split(const std::string &s, char delim);

// split a string on a single delimiter character, but only find the first n splits.
std::vector<std::string>& split(const std::string &s, char delim, std::vector<std::string> &elems, int &n);
std::vector<std::string>  split(const std::string &s, char delim, int n);

// split a string on any character found in the string of delimiters (delims)
std::vector<std::string>& split(const std::string &s, const std::string& delims, std::vector<std::string> &elems);
std::vector<std::string>  split(const std::string &s, const std::string& delims);

#endif // SPLIT_H
