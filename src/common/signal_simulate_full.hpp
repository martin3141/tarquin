#ifndef __SIGNALSIMULATE__
#define __SIGNALSIMULATE__

#include "../common/common.hpp"
#include "../common/CSignal.hpp"
#include "CCSVFile.hpp"

namespace tarquin {

	/*! The function to call to simulate a single CSignal object, e.g. a metabolite. The
	 * simulation parameters are chosen from fidMatch and the simulation details are specified by the
	 * CSV file.
	 * \param doubmat is the parsed CSV file of chemicals shifts, etc.
	 * \param signal is the result
	 * \param fidMatch contains the parameters of the simulation.
	 */
    int signal_simulate_full(std::vector<std::vector<double> >& full_doubmat, CSignal& signal, const CFID& fidMatch, const Options& options, std::string name);
	
}

#endif
