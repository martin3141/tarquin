#include <stdlib.h>
#include "../common/CBasis.hpp"
#include "../common/common.hpp"

#include <stdlib.h>
#include <algorithm>
#include <complex>
#include <sstream>

using namespace tarquin;
using namespace std;

int main(int argc, char* argv[])
{
	// cout << argv[1];
	// create the dummy FID object
	CFID fid;

	fid.SetSamplingFrequency(1000.0);
	//fid.SetSamplingFrequency(7200.0);
	fid.SetTransmitterFrequency(12.7789e7);
	//fid.SetTransmitterFrequency(600e6);
	//fid.SetEchoTime(35e-3);
	fid.SetEchoTime(0.060);
	//fid.SetEchoTime(0);
	//fid.SetEchoTime(0);
	fid.SetPPMRef(4.7);
	fid.SetNumberOfPoints(1024*4);
    fid.SetPulseSequence("press");

	// create the basis object
	CBasis basis;

	// simulate a basis to match this FID
	// basis.SimulateCSV("../../basis/phd/Lac.csv", fid);
	basis.SimulateCSV(argv[1], fid);

	// define matrix of signals (no groups) G 
	cvm::cmatrix G = basis.GetGroupMatrix();

	// get labels for plotting	
	std::vector < std::string > signal_names = basis.GetSignalNames();	
	// do fft
	cvm::cmatrix freq_sig(fft(G));
	// fftshift	
	freq_sig = fftshift(freq_sig);
	// get PPM scale
	cvm::rvector freq_scale = fid.GetPPMScale();
	// and plot
	//plot(freq_scale, freq_sig, signal_names);
	plot(freq_scale, freq_sig);

	return 0;
}
