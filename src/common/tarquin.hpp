#ifndef __TARQUINHEADER__
#define __TARQUINHEADER__

#include "CFID.hpp"
#include "CBasis.hpp"
#include "Workspace.hpp"
#include "Options.hpp"
#include "CBoswell.hpp"

namespace tarquin 
{

	/*!
	 * This is the main function for generating the results. 
	 *
	 * \param workspace must be initialised prior to use and when finished will contain the results.
	 * \param log is where all messages go.
	 */
	bool RunTARQUIN(Workspace& work, CBoswell& log);

}

#endif
