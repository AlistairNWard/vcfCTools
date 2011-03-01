// ******************************************************
// vcfCTools (c) 2011 Alistair Ward
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ------------------------------------------------------
// Last modified: 18 February 2011
// ------------------------------------------------------
// string split method from
// Evan Teran, http://stackoverflow.com/questions/236129/
// how-to-split-a-string/236803#236803
// ******************************************************

#include "split.h"

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <stdlib.h>

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while(std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    return split(s, delim, elems);
}

// Only split the first n times.
std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems, int &n) {
    std::stringstream ss(s);
    std::string item;
    int counter = 1;
    while(std::getline(ss, item, delim) && counter <= n) {
      elems.push_back(item);
      counter++;
    }

    return elems;
}

std::vector<std::string> split(const std::string &s, char delim, int n) {
    std::vector<std::string> elems;
    return split(s, delim, elems, n);
}

// Split on entries in a delimeter string.
std::vector<std::string> &split(const std::string &s, const std::string& delims, std::vector<std::string> &elems) {
    char* tok;
    char cchars [s.size()+1];
    char* cstr = &cchars[0];
    strcpy(cstr, s.c_str());
    tok = strtok(cstr, delims.c_str());
    while (tok != NULL) {
        elems.push_back(tok);
        tok = strtok(NULL, delims.c_str());
    }
    return elems;
}

std::vector<std::string> split(const std::string &s, const std::string& delims) {
    std::vector<std::string> elems;
    return split(s, delims, elems);
}
