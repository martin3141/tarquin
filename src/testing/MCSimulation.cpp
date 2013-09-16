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
    // create the dummy FID object
    CFID fid;

    fid.SetSamplingFrequency(1000.0);
    fid.SetTransmitterFrequency(6.3866e7);
    fid.SetEchoTime(30e-3);
    fid.SetPPMRef(4.7);
    fid.SetNumberOfPoints(1024);

    // create the basis object
    CBasis basis;

    // simulate a basis to match this FID
    //basis.Simulate("../../basis/phd/", fid);
	
	// load an xml basis file
	//load_xml(basis, "../../basis/phd.xml");
	//basis.GenerateNonSTLstruc();

    // define matrix of signals (no groups) S 
    const cvm::cmatrix& Q = basis.GetBasisMatrix();
    cvm::cmatrix S = Q;

    // initialise random number generator to the system time
    srand(time(NULL));

    // a represents signal amplitudes, generate from gauss dist. rand numbers
    cvm::cvector a(S.nsize());

    // randgen - number of random spectra to generate
    int randgen = 100;
    cvm::cvector y(S.msize());

	for(int m = 1; m <= randgen ;m++) {
		for(int n = 1; n < S.nsize(); n++) {	

			// mean 10 stdev 5 (amplitude of basis column)
			a(n) = GeneratePosGauss(10,15);

		}

		/*// apply some alpha peturbations
		  for(int c = 1; c < S.nsize(); c++ ) {

		//treal alpha = GeneratePosGauss(50, 2) -48;
		treal alpha = 5;

		treal dt = 1.0 / fid.GetSamplingFrequency();

		for( int n = 1; n < S.msize(); n++ ) 
		S(n,c) = S(n,c) * exp( -alpha * (n-1) * dt );
		}*/

		// apply some frequency perturbations
		/*for(int c = 1; c < S.nsize(); c++ ) {

		  treal dt = 1.0 / fid.GetSamplingFrequency();

		  for( int n = 1; n < S.msize(); n++ ) {

		  tcomplex jomega = tcomplex(0, 2.0 * M_PI * 0.5 * dt * (n-1));

		  S(n,c) = S(n,c) * exp( jomega );
		  }
		  }*/

		// add components together
		y = (S*a);

		// change phase of the signal	
		double phi0 = 0;
		double phi1 = 0.001;	
		// number of points in the signal
		integer N = y.size();
		// make the DFT of the time signal
		cvm::cvector Y(N);
		fft(y, Y);
		Y = fftshift(Y);
		// get the shifted [-fs/2, ..., 0, ... fs/2) frequency range
		cvm::rvector freq_range = fid.GetFreqScale();
		for( integer n = 0; n < N; n++ ) {
			Y[n+1] = Y[n+1] * exp( tcomplex(0, phi0 + phi1*2.0*M_PI*freq_range[n+1]) );
		}
		Y = fftshift(Y);
		// go back to normal range
		ifft(Y, y);


		// add noise
		y += GenerateNoise(y.size(),2);

		// convert m to string
		stringstream ss;
		string FileNum;
		ss << m;
		ss >> FileNum;

		// now save y as .dpt format
		string strFilename = "../../data/simulated/";
		strFilename.append(FileNum);
		strFilename.append(".dpt");

		CFID dptout = fid;
		// change by greg to reflect new structure
		//for( int n = 0; n < S.msize(); n++)
		//	dptout.appendSample(y(n+1));
		dptout.InitialiseFromVector(y);

		dptout.SaveToFile(strFilename);

		// fn creates dangerplot file at strFilename	
		//		WriteDanger(expParams, y, strFilename);

		// redefine strFilename and save signal aplitudes
		strFilename = "../../data/simulated/";
		strFilename.append(FileNum);
		strFilename.append(".txt");
		ofstream fout(strFilename.c_str(), ios::out);
		for(integer n = 1; n < S.nsize(); n++) 
			fout << basis.GetSignalName(n) << "\t" << a(n).real() << endl;
		// close file?

	}
	// plot the last one	
	// do fft
	//cvm::cmatrix freq_sig(fft(y));
	// fftshift	
    //freq_sig = fftshift(freq_sig);
    // get PPM scale
    //cvm::rvector freq_scale = fid.GetPPMScale();
    // and plot
    //plot(freq_scale, freq_sig);


    return 0;
}
