#include "common.hpp"
#include "CFID.hpp"
#include "Workspace.hpp"
#include "cvm_util.hpp"

namespace tarquin {

    std::string g_strGnuPlot = "gnuplot";

    // this should probably shift by the reference amount as well
    void PlotFIDs(const CFID& fidA, const CFID& fidB)
    {
        const cvm::cvector& yA = fidA.GetVectorFID(0); //TODO
        const cvm::cvector& yB = fidB.GetVectorFID(0); //TODO
        assert( yA.size() == yB.size() );

        integer N = yA.size();

        cvm::cvector YA(N);
        cvm::cvector YB(N);

        fft(yA, YA);
        fft(yB, YB);

        // obtain PPM scale from fidB
        coord vox(1, 1, 1); 
        cvm::rvector freq_scale = fidB.GetPPMScale(vox);

        // make labels	
        std::vector<std::string> labels;
        labels.push_back("preprocessed");
        labels.push_back("processed");

        cvm::cmatrix B(N, 2);
        B(1) = YA;
        B(2) = YB;
        B = fftshift(B);
        plot(freq_scale, B, labels);
    }

    void lb(cvm::cvector& y, const CFID& fid, double lb)
    {
        if ( lb > 0 )
        {
            // apply some line broadening
            int N = fid.GetNumberOfPoints();
            // the sampling interval (time step)
            double dt = 1.0 / fid.GetSamplingFrequency(); 

            // generate time signal
            cvm::rvector t(N);
            for(int n = 0; n < N; n++)
                t(n+1) = n*dt;

            // don't do the last 5 data points
            for(int n = 0; n < N-5; n++)
                y(n+1) = y(n+1)*exp(-lb*M_PI*t(n+1));
        }
    }
    
    void ZeroPad(cvm::cvector& y, int factor)
	{
		// zero fill if less than 4096 points
		//if( y.size() < 4096 )
		//	zf = round(4096/y.size());

		if( factor > 1 )
            y.resize(y.size()*factor);

		// copy last pts points of FID to end of zfilled fid so we don't get
		// any nasty Gibbs artifacts
		int pts = 5;
		if( factor > 1 )
		{
			for( int n = 1; n < pts + 1; ++n )
			{
				y(y.size()-pts+n)    = y(y.size()/factor-pts+n);
				y(y.size()/factor-pts+n) = 0;
			}
		}
	}



    

}
