#ifndef __WriteDanger__
#define __WriteDanger__

#include <iostream>
#include "../common/common.hpp"
#include "cvm.h"

using namespace std;
//namespace fs = boost::filesystem;
namespace tarquin {

	void WriteDanger(CExperiment expParas, cvm::cvector FID, string strFilename) {

		ofstream fout(strFilename.c_str(), ios::out);
		
		// don't like these formatting options
		// setup the stream so that things are formatted properly
		//fout.setf(ios::scientific |ios::left);
		//fout.precision(8);

		fout << "Dangerplot_version\t" << "1.0" << endl;
		fout << "Number_of_points\t" << expParas.getNumberOfPoints() << endl;
		fout << "Sampling_frequency\t" << expParas.getSamplingFrequency() << endl;
		fout << "Transmitter_frequency\t" << expParas.getTransmitterFrequency() << endl;
		fout << "Phi0\t" << 0 << endl;
		fout << "Phi1\t" << 0 << endl;
		fout << "PPM_reference\t" << expParas.getReferenceOffset()
 << endl;
		fout << "Real_FID\t" << "Imag_FID\t" << endl;
		for(int n = 1; n <= (expParas.getNumberOfPoints()); n++)
			fout << FID(n).real() << "\t" << FID(n).imag() << endl;
		// close file?
	}
}

#endif
