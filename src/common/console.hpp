#ifndef __CONSOLE__
#define __CONSOLE__

#include "Options.hpp"
#include "CFID.hpp"

/*!
 * A collection of functions for the command line interface to TARQUIN. They initialise the options
 * and FID structures so the algorithm can run.
 */
namespace tarquin {

    void DisplaySplash();

    void DisplayUsage();

    bool ParseCommandLine(int argc, char* argv[], Options& options, CFID& fid);
}

#endif
