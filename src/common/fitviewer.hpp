#ifndef __FITVIEWER__
#define __FITVIEWER__

#include <string>

namespace tarquin 
{

	// prototpe for fitviewer program, should probably get its own file
	int FitView(std::string strFilename, std::string stdPlotSigs, treal ppm_start, treal ppm_end, treal lb, bool pause);

	int FitView(
		const Workspace& workspace, 
		std::string strPlotSigs, 
		treal ppm_start, 
		treal ppm_end, 
		treal lb,
		bool pause_on_screen,
		std::string filename
		);

}

#endif
