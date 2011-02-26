// ******************************************************
// vcfCTools (c) 2011 Alistair Ward
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ------------------------------------------------------
// Last modified: 18 February 2011
// ------------------------------------------------------
// Base class for tools.
// All must contain Help() and Run()
// ******************************************************

#ifndef VCFCTOOLS_ABSTRACTTOOL_H
#define VCFCTOOLS_ABSTRACTTOOL_H

#include <string>

namespace vcfCTools {

class AbstractTool {
  public:
    AbstractTool( void ) {}
    virtual ~AbstractTool( void ) {}

  public:
    virtual int Help( void ) = 0;
    virtual int Run( int argc, char* argv[] ) = 0;
};

} // namespace vcfCTools

#endif // VCFCTOOLS_ABSTRACTTOOL_H
