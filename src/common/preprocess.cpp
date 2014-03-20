#include "preprocess.hpp"
#include "Workspace.hpp"
#include "CBoswell.hpp"
#include "common.hpp"
#include "hsvd.hpp"
#include "Options.hpp"
#include "td_conv_ws.hpp"
//#include "../levmar/levmar-2.1.3/lm.h"
#include "../levmar/levmar-2.5/levmar.h"
#include <limits>
#include "signal_simulate_full.hpp"


namespace tarquin
{
void phasing_objective(treal* pp, treal* pYP, integer N, void* pParams);
void residual_objective_ref(treal* pp, treal* pyhat, integer nParams, integer activeNdbl, void* pParams);
    /*!
     * Attempt to automatically phase the passed in FID. 
     */
    void AutoPhaseSimple(const coord& proc_coord, CFID& fid, CBoswell& log);
    
    /*!
     * Apply phase paramters to the FID. 
     */
    void ApplyPhase(const coord& proc_coord, CFID& fid, CBoswell& log);
    void ApplyPhase(const coord& proc_coord, CFID& fid, CFID& av_fid, CBoswell& log);

    /*!
     * Attempt to automatically reference the passed in FID.
     */
    bool AutoReference(const coord& proc_coord, CFID& fid, CBoswell& log);

    /*!
     * Attempt to automatically reference the passed in FID.
     */
    bool AutoReferenceCorr(const coord& proc_coord, Options& options, CFID& fid, bool dyn_mode, bool add_h2o, bool skip_beta_guess, CBoswell& bos_log);

    /*!
     * Do eddy current correction using water reference signal.
     */
    void WaterEddy(const coord& proc_coord, CFID& fidproc, const CFID& fidwater, CBoswell& log);
    
    void WaterEddySimple(const coord& proc_coord, CFID& fidproc, const CFID& fidwater, treal water_freq, CBoswell& log);

	/*!
     * Find the amplitude of a fitted single Voigt peak, this is the normalisation value.
     */
	treal ComputeWaterNormalisation(const coord& proc_coord, const CFID& fid, CBoswell& log);

}

tarquin::Preprocessor::Preprocessor(Workspace& workspace, CBoswell& log) :
	m_workspace(workspace),
	m_log(log) 
{ 
}


void tarquin::Preprocessor::operator() ()
{
	Task progress(m_log, "Preprocessing");

	const CFID& fidraw            = m_workspace.GetFIDRaw();
	Options& options              = m_workspace.GetOptions();
	CFID& fidproc                 = m_workspace.GetFIDProc();
	CFID& fidwater                = m_workspace.GetFIDWater();

	// initialise preprocessed FID with raw FID's data (all members are copied)
	fidproc = fidraw;
	
	const coord_vec& fit_list = options.GetFitList();

    bool combine_preproc = options.GetCombinePreproc();
    
    bool lipid_filter = options.GetLipidFilter();

    // calculate an average spectrum for preprocessing
    CFID fidav = fidproc;
    if ( combine_preproc )
    {
        cvm::cvector yav( fidproc.GetNumberOfPoints() );
        for( coord_vec::const_iterator i = fit_list.begin(); i != fit_list.end(); ++i )	
            yav += fidproc.GetVectorFID(*i);
    
        fidav.SetRows(1);
        fidav.SetCols(1);
        fidav.SetSlices(1);
        cvec_stdvec& avfid = fidav.GetVectorFID();
        avfid.clear();
        avfid.push_back(yav);
        const coord start_spec(1, 1, 1); 
        
        // are we doing water SVD based removal?
		if( options.GetWaterWindow() > 0 ) 
		{
			cvm::rvector frequencies;
			cvm::rvector dampings;
			cvm::cmatrix basis;
			cvm::cvector ahat;

			int pts; 
			if ( fidraw.GetNumberOfPoints() > 1024 )
				pts = 1024;
			else
				pts = fidraw.GetNumberOfPoints();

            if ( pts == fidraw.GetNumberOfPoints() )
                pts = pts * 0.9;

			// decompose into linear combination of complex exponentials
			//DecomposeHSVD(start_spec, fidav, frequencies, dampings, basis, ahat, pts, 50, m_log);
			DecomposeHSVD(yav, fidraw.GetSamplingFrequency(), frequencies, dampings, basis, ahat, pts, 50, m_log);

			assert( ahat.size() == frequencies.size() );
			assert( basis.msize() == yav.size() );

			// set all frequencies outside of window to zero 
			// (we are building a model of the water to subtract)
			treal width = options.GetWaterWindow();
            
            treal ref = fidproc.GetPPMRef(0);
            //std::cout << "REF = " << ref << std::endl;
            treal trans_freq = fidproc.GetTransmitterFrequency();
            treal lip_freq = ( ref - options.GetLipFilterFreq() ) * trans_freq * 1e-6;
            //treal lip_freq = ( ref - 0.0 ) * trans_freq * 1e-6;

			for( integer i = 0; i < frequencies.size(); i++ )
            {
                if ( lipid_filter )
                {
                    if ( ( frequencies[i+1] > width ) && ( frequencies[i+1] < lip_freq ))
                        ahat[i+1] = 0;
                }
                else
                {
                    if( frequencies[i+1] > width ) 
                        ahat[i+1] = 0;
                }

				// kill increasing signals
                //if ( !options.GetFullEcho() )
                if( dampings[i+1] > 0 )
                    ahat[i+1] = 0;
			}

			// make the model of the water
			cvm::cvector yw = basis*ahat;

			// subtract the water
			yav = yav - yw;
		}

		// are we doing convolution based water removal?
		if( options.GetConvWindowWidth()>0 ) 
		{
			Task progress(m_log, "water removal");

			// estimated water signal
			cvm::cvector yw(yav.size());
			int K = options.GetConvWindowWidth();

			// make the model of the water yw
			td_conv_ws(yav, yw, K, 16);

			// subtract the water
			yav = yav - yw;
		}

        avfid.clear();
        avfid.push_back(yav);
        
       
        /*
        if ( options.Getlb() != 0 )
        {
            // apply some line broadening
			int N = fidproc.GetNumberOfPoints();
			// the sampling interval (time step)
			treal dt = 1.0 / fidproc.GetSamplingFrequency(); 

			// generate time signal
			cvm::rvector t(N);
			for(integer n = 0; n < N; n++)
                t(n+1) = n*dt;

            treal lb = options.Getlb();
            for(integer n = 0; n < N; n++)
                yav(n+1) = yav(n+1)*exp(-lb*M_PI*t(n+1));
        }
        */
        
        
        if( options.GetAutoPhase() ) 
		{
			//AutoPhaseSimple(start_spec, fidav, m_log);
			//AutoPhase(start_spec, fidav, options, m_log);
			AutoPhaseNew(start_spec, fidav, options, m_log);
		}

        /*
        cvm::cvector yav_phased = fidav.GetVectorFID(start_spec);
        cvm::cvector Yav(yav.size());
        fft(yav_phased, Yav);
        Yav = fftshift(Yav);
        plot(Yav);
        */
    }

	int voxel_num = 0;
	// iterate over the fit list	
	for( coord_vec::const_iterator i = fit_list.begin(); i != fit_list.end(); ++i )	
	{
		voxel_num++;
		m_log.LogMessage(LOG_INFO, "Preprocessing fid %u of %u", voxel_num, fit_list.size());


		cvm::cvector& y = fidproc.GetVectorFID(*i);

        // decide on normalisation value
		// it is the biggest point in the magnitude spectrum
		if( options.GetFilenameWater() != "" ) 
		{
            // compute the value
			treal norm_val = ComputeWaterNormalisation(*i, fidwater, m_log);
            // get a reference
			std::vector<double>& norm_val_vec = m_workspace.GetNormalisationValue();
            // append value to vector
            norm_val_vec.push_back(norm_val);
            
            // phase water signal
			//AutoPhaseNew(*i, fidwater, options, m_log, true);

            // this one seems to do a better job for water
			AutoPhase(*i, fidwater, options, m_log, true);

            // find the FWHM
            cvm::cvector y = fidwater.GetVectorFID(*i);
            integer zf = 8;
            // apply zero filling
            integer N = y.size() * zf;
            y.resize(N);
            treal fs = fidwater.GetSamplingFrequency();
            treal ft = fidwater.GetTransmitterFrequency();
            // DFT is at this frequency resolution
            treal df = fs / N;

            cvm::cvector Y(N);
            fft(y, Y);
            Y = fftshift(Y);
            cvm::rvector Yre = Y.real();

            // find the maximum data point in spectrum
            treal max_pt = Yre.indofmax();
            treal max_val = Yre(max_pt);
            treal hh = max_val/2.0;

            // find the right hh point
            treal right_hh_pt = -1;
            for ( int n = max_pt; n < Yre.size(); n++ )
            {
                if ( Yre(n) >= hh && Yre(n+1) < hh )
                {
                    treal m = ( Yre(n) - Yre(n+1) ) / ( n - (n+1) );
                    treal c = Yre(n) - m * n;
                    right_hh_pt = (hh - c) / m;
                    break;
                }
            }

            // find the left hh point
            treal left_hh_pt = -1;
            for ( int n = max_pt; n > 1; n-- )
            {
                if ( Yre(n) >= hh && Yre(n-1) < hh )
                {
                    treal m = ( Yre(n) - Yre(n-1) ) / ( n - (n-1) );
                    treal c = Yre(n) - m * n;
                    left_hh_pt = (hh - c) / m;
                    break;
                }
            }

            treal width_hz = ( right_hh_pt - left_hh_pt ) * fs / N;

	        cvm::rvector freq_range = fidwater.GetFreqScale(zf);
            treal water_freq = freq_range(max_pt);

            //std::cout << std::endl << "Water width (Hz) = " << width_hz << std::endl;
            //plot(Y);

            m_workspace.AppendWaterWidth(width_hz);
            m_workspace.AppendWaterFreq(water_freq);

		}

		// we will be doing no normalisation
		else 
		{
            // get a reference
			std::vector<double>& norm_val_vec = m_workspace.GetNormalisationValue();
            // append value to vector
            norm_val_vec.push_back(1.0);
            m_workspace.AppendWaterWidth(-1);
            m_workspace.AppendWaterFreq(0);
		}

        
        // are we using the water for eddy current correction?
		if( options.GetWaterEddy() && options.GetFilenameWater() != "")
		{
            treal water_freq = m_workspace.GetWaterFreq(voxel_num-1);
            //std::cout << std::endl << "Water freq (Hz) = " << water_freq << std::endl;

			WaterEddySimple(*i, fidproc, fidwater, water_freq, m_log);
		} 
        
        // pre shift spectra to center the water peak?
        treal old_ref = fidproc.GetPPMRef(*i);
        treal new_ref;
        
        if ( options.GetPreWsShift() )
        {
            AutoReferenceCorr(*i, options, fidproc, false, true, true, m_log);
            new_ref = fidproc.GetPPMRef(*i);
            fidproc.ShiftRef( old_ref, *i);
        }
        
        /*
        std::cout << old_ref << std::endl;
        std::cout << new_ref << std::endl;
        */
        
        // are we doing convolution based water removal?
		if( options.GetConvWindowWidth()>0 ) 
		{
			Task progress(m_log, "water removal");

			// estimated water signal
			cvm::cvector yw(y.size());
			int K = options.GetConvWindowWidth();

			// make the model of the water yw
			td_conv_ws(y, yw, K, 16);

			// subtract the water
			y = y - yw;
		}

        
        // are we doing a pre SVD based removal?
		if( options.GetPreHSVD() ) 
		{
			cvm::rvector frequencies;
			cvm::rvector dampings;
			cvm::cmatrix basis;
			cvm::cvector ahat;

			int pts; 
			if ( fidraw.GetNumberOfPoints() > 1024 )
				pts = 1024;
			else
				pts = fidraw.GetNumberOfPoints();

            if ( pts == fidraw.GetNumberOfPoints() )
                pts = pts * 0.9;

			// decompose into linear combination of complex exponentials
			DecomposeHSVD(y, fidraw.GetSamplingFrequency(), frequencies, dampings, basis, ahat, pts, 5, m_log);

			assert( ahat.size() == frequencies.size() );
			assert( basis.msize() == y.size() );

			// set all frequencies outside of window to zero 
			// (we are building a model of the water to subtract)
			treal width = options.GetWaterWindow();
			//treal lip_filt = options.GetWaterWindow();

            treal ref = fidproc.GetPPMRef(*i);
            treal trans_freq = fidproc.GetTransmitterFrequency();
            treal lip_freq = ( ref - options.GetLipFilterFreq() ) * trans_freq * 1e-6;

			for( integer i = 0; i < frequencies.size(); i++ ) 
            {
                
                if( frequencies[i+1] > width ) 
                    ahat[i+1] = 0;

				// kill increasing signals
                //if ( !options.GetFullEcho() )
                if( dampings[i+1] > 0 )
                    ahat[i+1] = 0;
			}

			// make the model of the water
			cvm::cvector yw = basis*ahat;

			// subtract the water
			y = y - yw;
		}


		// are we doing water SVD based removal?
		if( options.GetWaterWindow() != 0 ) 
		{
			cvm::rvector frequencies;
			cvm::rvector dampings;
			cvm::cmatrix basis;
			cvm::cvector ahat;

			int pts; 
			if ( fidraw.GetNumberOfPoints() > 1024 )
				pts = 1024;
			else
				pts = fidraw.GetNumberOfPoints();
            
            if ( pts == fidraw.GetNumberOfPoints() )
                pts = pts * 0.9;

			// decompose into linear combination of complex exponentials
			DecomposeHSVD(y, fidraw.GetSamplingFrequency(), frequencies, dampings, basis, ahat, pts, 50, m_log);

			assert( ahat.size() == frequencies.size() );
			assert( basis.msize() == y.size() );

			// set all frequencies outside of window to zero 
			// (we are building a model of the water to subtract)
			treal width = options.GetWaterWindow();
			//treal lip_filt = options.GetWaterWindow();

            treal ref = fidproc.GetPPMRef(*i);
            treal trans_freq = fidproc.GetTransmitterFrequency();
            treal lip_freq = ( ref - options.GetLipFilterFreq() ) * trans_freq * 1e-6;

			for( integer i = 0; i < frequencies.size(); i++ ) 
            {
                
                if ( lipid_filter )
                {
                    if ( ( frequencies[i+1] > width ) && ( frequencies[i+1] < lip_freq ))
                        ahat[i+1] = 0;
                }
                else
                {
                    if( frequencies[i+1] > width ) 
                        ahat[i+1] = 0;
                }

				// kill increasing signals
                //if ( !options.GetFullEcho() )
                if( dampings[i+1] > 0 )
                    ahat[i+1] = 0;
			}

			// make the model of the water
			cvm::cvector yw = basis*ahat;

            y = y - yw;

            /*
            if ( options.GetPreHSVD() )
            {
                cvm::cvector y_new = y - yw;
                cvm::cvector Y_new(y.size());
                fft(y_new, Y_new);
                cvm::cvector Y(y.size());
                fft(y, Y);
                std::cout << std::endl << Y_new.norm() << std::endl;
                std::cout << std::endl << Y.norm() << std::endl;
                if ( Y_new.norm() < Y.norm() )
                    y = y - yw;
            }
            else
            {
                y = y - yw;
            }
            */

		}

		        
        // now the water supression has been done
        // try and figure out nEnd if necessary
        if ( options.GetRangeEnd() == -1 )
        {
            cvm::rvector decay(y.size());
            for ( int n = 1; n < y.size() + 1; n++ )
                decay(n) = abs(y(n));

            cvm::rvector decay_smooth(y.size());
            td_conv_ws(decay, decay_smooth, 50, 10);
            plot(decay);
            plot(decay_smooth);
            
            treal max_decay = decay_smooth(decay_smooth.indofmax());
            treal min_decay = decay_smooth(decay_smooth.indofmin());

            int end_pt = -1;
            for ( int n = 1; n < decay_smooth.size(); n++ )
            {
                if ( ( decay_smooth(n) - min_decay > (max_decay-min_decay)*0.01 ) && ( decay_smooth(n+1) -min_decay <= (max_decay-min_decay)*0.01 ) )
                {
                    end_pt = n;
                    break;
                }
            }
            
            if ( end_pt == -1 )
                end_pt = y.size() - 10;

		    m_log.LogMessage(LOG_INFO, "Auto end pt set to : %i", end_pt);
		    options.SetRangeEnd(end_pt);
        }

        // are we broadening the spectrum?
        /*
        if ( options.Getlb() != 0 )
        {
            // apply some line broadening
			int N = fidproc.GetNumberOfPoints();
			// the sampling interval (time step)
			treal dt = 1.0 / fidproc.GetSamplingFrequency(); 

			// generate time signal
			cvm::rvector t(N);
			for(integer n = 0; n < N; n++)
                t(n+1) = n*dt;

            treal lb = options.Getlb();
            for(integer n = 0; n < N; n++)
                y(n+1) = y(n+1)*exp(-lb*M_PI*t(n+1));
        }*/

        // shift back to old value 
        if ( options.GetPreWsShift() )
        {
            fidproc.ShiftRef( new_ref, voxel_num-1);
            fidproc.SetPPMRef(*i, old_ref);
        }

      	// are we doing automatic phasing?
		if( options.GetAutoPhase() && !combine_preproc ) 
		{
			//AutoPhaseSimple(*i, fidproc, m_log);

            // this algo works pretty well for most
            // data changed to default on 
            // 27th Jun 12
			AutoPhaseNew(*i, fidproc, options, m_log);

            /*
			if ( fidraw.GetNumberOfPoints() > 8192 )
				AutoPhaseSimple(*i, fidproc, m_log);
			else
				AutoPhase(*i, fidproc, options, m_log);
                */
		}
		else if( options.GetAutoPhase() && combine_preproc ) 
        {
			ApplyPhase(*i, fidproc, fidav, m_log);
        }
		// no, manual phasing
		else 
		{
			ApplyPhase(*i, fidproc, m_log);
		}

		// are we doing automatic referencing?
		// this can be done before phasing if corr is used
		// this way zero filling can be used without danger of
		// final point problems
	    /*		
		if( options.GetAutoReference() ) 
		{
			AutoReferenceCorr(*i, options, fidproc, m_log);
			//AutoReference(fidproc, log);
		}*/
        
        // For now always do this to make sure the init_beta values are correctly set
        // the ref will not be changed if !options.GetAutoReference()
        AutoReferenceCorr(*i, options, fidproc, false, false, false, m_log);

		
		/*
		   cvm::cvector yy = y;
		   cvm::cvector YY(yy.size());
		   fft(yy, YY);
		   cvm::rvector YY_mag(yy.size());
		   for ( int n=1; n < YY.size()+1; n++ )
		   YY_mag(n) = pow(pow(YY(n).real(),2) + pow(YY(n).imag(),2),0.5);

		   YY_mag = fftshift(YY_mag);
		// find the average of the first 50 freq pts
		treal bg_mag = 0;
		for ( int n = 1 ; n < 50 + 1; n++)
		{
		bg_mag += YY_mag(n)/50;
		}

		for ( int n=1; n < yy.size()+1; n++ )
		YY_mag(n) -= bg_mag;

		cvm::rvector YY_mag_smo(y.size());
		td_conv_ws_noext( YY_mag, YY_mag_smo, 50, 10);	

		cvm::rmatrix mat(YY_mag.size(),2);
		mat(1) = YY_mag;
		mat(2) = YY_mag_smo;
		plot(mat);

		// find biggest dif between YY_mag and YY_mag_smo
		treal max_n_diff = 0;
		treal max_diff = 0;
		for ( int n=1; n < YY.size()+1; n++ )
		if ( YY_mag_smo(n) - YY_mag(n) > max_diff )
		{
		max_diff = YY_mag_smo(n) - YY_mag(n);
		max_n_diff = n;
		}

		treal scale =  YY_mag(max_n_diff) / YY_mag_smo(max_n_diff);

		for ( int n=1; n < YY.size()+1; n++ )
		YY_mag_smo(n) = YY_mag_smo(n) * scale;

		mat(1) = YY_mag;
		mat(2) = YY_mag_smo;
		plot(mat);

		std::cout << max_n_diff << std::endl;
		 */


		/*
		// get the magnitude spectrum
		cvm::cvector yy = y;
		cvm::rvector yy_mag(yy.size());
		for ( int n=1; n < yy.size()+1; n++ )
		yy_mag(n) = pow(pow(yy(n).real(),2) + pow(yy(n).imag(),2),0.5);

		cvm::rvector yy_mag_smo(y.size());
		td_conv_ws_noext( yy_mag, yy_mag_smo, 10, 10);	

		//plot(yy_mag_smo);
		//plot(yy_mag);


		cvm::cvector YY(yy.size());
		cvm::rvector YY_mag(yy.size());
		cvm::rvector YY_mag_smo(yy.size());

		treal mag_sum = 0;
		int steps = 200;
		cvm::rvector step_vec(steps);

		for ( int m=0; m < steps; m++ ) 
		{
		yy(1) = yy(1) * 0.5;
		fft(yy, YY);
		for ( int n=1; n < yy.size()+1; n++ )
		YY_mag(n) = pow(pow(YY(n).real(),2) + pow(YY(n).imag(),2),0.5);

		YY_mag = fftshift(YY_mag);
		//plot(YY_mag);
		//td_conv_ws_noext( YY_mag, YY_mag_smo, 200, 10);	
		//cvm::rmatrix mat(YY_mag.size(),2);
		//mat(1) = YY_mag;
		//mat(2) = YY_mag_smo;
		//plot(mat);


		// find the average of the first 50 freq pts
		//treal bg = 0;
		//for ( int n = 1 ; n < 50 + 1; n++)
		//{
		//    bg += YY_mag(n)/50;
		//}

		//for ( int n=1; n < yy.size()+1; n++ )
		//    YY_mag(n) -= bg;


		//plot(YY_mag);

		// calc sum
		mag_sum = 0;
		for ( int n=1; n < YY_mag.size()+1; n++ )
		mag_sum = mag_sum + YY_mag(n);
		//std::cout << mag_sum << std::endl;
		step_vec(m+1) = mag_sum;
		// remove first element of yy
		yy.erase(yy.begin());
		yy.resize(yy.size());
		YY_mag.resize(yy.size());
		//plot(YY_mag);
		}
		cvm::rvector step_vec_smo(steps);
		cvm::rvector der_step_vec(steps-1);
		td_conv_ws_noext( step_vec, step_vec_smo, 4, 10);	
		plot(step_vec_smo);
		plot(step_vec);

		treal dt = 1.0 / fidproc.GetSamplingFrequency(); 
		for ( int m=1; m < steps; m++ ) 
		der_step_vec(m) = (-step_vec_smo(m)+step_vec_smo(m+1))/dt;

		//plot(der_step_vec);
		*/

			/*
			// apply some line broadening
			int N = fidproc.GetNumberOfPoints();
			// the sampling interval (time step)
			treal dt = 1.0 / fidproc.GetSamplingFrequency(); 

			// generate time signal
			cvm::rvector t(N);
			for(integer n = 0; n < N; n++)
			t(n+1) = n*dt;

			treal lb = options.Getlb();
			treal beta = pow((lb/2),2)/-0.69314718;
			std::cout << beta << std::endl;

			for(integer n = 0; n < N; n++)
			y(n+1) = y(n+1)*exp(beta*t(n+1)*t(n+1));
			 */

			
		/*int zf = 1;
		// zero fill if less than 4096 points
		if ( y.size() < 4096 )
			zf = static_cast<int>(4096/y.size());
            */

		//zf = 1;

        // Guess the SNR of the signal to find a good value for nStart
        // measure of fit quality
        cvm::cvector yz = y; 

        // zero fill
		yz.resize(yz.size()*options.GetZF());

		// copy last pts points of real fid to end of zfilled fid
		int pts = 5;
		if ( options.GetZF() > 1 )
		{
			for ( int n = 1 ; n < pts + 1; n++)
			{
				yz(yz.size()-pts+n) = yz(yz.size()/options.GetZF()-pts+n);
				yz(yz.size()/options.GetZF()-pts+n) = 0;
			}
		}

		cvm::cvector Y(y.size());
		Y.resize(yz.size());

		fft(yz, Y);

		Y = fftshift(Y);

        //plot(Y);

		// find the average of the first 10 freq pts
		treal bg = 0;
		for ( int n = 5 ; n < 15 + 1; n++)
		{
			bg += abs(Y(n))/10;
		}

		//std::cout << std::endl << bg << std::endl << std::flush;
		//plot(Y);	

		cvm::rvector freq_scale = fidproc.GetPPMScale(*i, options.GetZF());

		// find points corresponding to ppm start and ppm end
		int left = 1, right = freq_scale.size()-1;
		for ( int n = 1; n < (freq_scale.size()); n++ ) 
		{
			if ( ( freq_scale(n) > options.GetPPMend() ) && ( freq_scale(n+1) <= options.GetPPMend() ) )
				left = n;
			if ( ( freq_scale(n) > options.GetPPMstart() ) && ( freq_scale(n+1) <= options.GetPPMstart() ) )
				right = n;
		}

        if ( left == 1 || right == (freq_scale.size() - 1) )
            m_log.LogMessage(LOG_INFO, "Warning, one of the PPM limits is outside the bandwidth.");

		double Ymax = 0;
		// find max of y for SNR calculation
		for ( int n = left; n < right+1; n++ ) 
		{
			//if ( Ymax < Y(n).real() )
		    //	Ymax = Y(n).real();

			if ( Ymax < abs(Y(n)) )
				Ymax = abs(Y(n));
		}

        Ymax = Ymax - bg;
		//cvm::cvector BASELINE;
		//td_conv_ws( Y, BASELINE, options.GetBL()*options.GetZF(), 10);	
        //plot(Y);
        //Y = Y-BASELINE;

        cvm::rvector plot_vec(Y.size());
        cvm::rvector plot_vec_re(Y.size());
        for ( int n = 1; n < Y.size(); n++ )
        {
            plot_vec(n) = abs(Y(n)) - bg;
            plot_vec_re(n) = Y(n).real() - bg;
        }

        //plot(plot_vec);
        //plot(plot_vec_re);

		double noise_min = std::numeric_limits<double>::infinity();
		double noise_temp = 0;
		int block_size = 200;
		for ( int n = 1; n < (Y.size()+1) - block_size; n = n + block_size ) 
		{
			noise_temp = stdev(Y.real(),n,n+block_size-1);
			if ( noise_temp < noise_min )
				noise_min = noise_temp;
		}

		//m_log.DebugMessage(DEBUG_LEVEL_1, "stddev of noise = %.4f", noise_min);

        //std::cout << "Ymax = " << Ymax << std::endl;
		double SNR = Ymax/(2*noise_min);

		//fid.SetSNR(SNR);

        m_log.LogMessage(LOG_INFO, "Ymax       = %f", Ymax);
		m_log.LogMessage(LOG_INFO, "sdev noise = %f", noise_min);
		m_log.LogMessage(LOG_INFO, "SNR guess  = %.4f", SNR);

		treal fs = fidproc.GetSamplingFrequency();
		treal ft = fidproc.GetTransmitterFrequency();

		//int nStartAdpt = round(log(SNR) * fs / 500);
		//int nStartAdpt = round(log(SNR*0.3) * fs / (ft*1e-6*2.8));
		//int nStartAdpt = round(log(SNR*0.2) * fs / (250)); //works ok
		//int nStartAdpt = round(log(SNR*0.3) * fs / (310));
		//int nStartAdpt = round(log(SNR*0.4) * fs / (350));
		int nStartAdpt = round(log(SNR*0.5) * fs / (385));
        if ( nStartAdpt < options.GetMinRange() )
            nStartAdpt = options.GetMinRange();
        
        //std::cout << options.GetRangeStart() << std::endl;

	    if( 0 == options.GetRangeStart() )
        {
		    m_log.LogMessage(LOG_INFO, "Adaptive nStart = %i", nStartAdpt);
		    options.SetRangeStart(nStartAdpt);
        }

	}
}


// Do eddy current correction, as described in introduction of:
// "Artifacts Introduced by Zero Order Phase Correction in Proton
// NMR Spectroscopy and a Method of Elimination by Phase Filtering"
// by J. M. Wild, 1998, J. MR, V137, pp 430-436
void tarquin::WaterEddySimple(const coord& proc_coord, CFID& fidproc, const CFID& fidwater, treal water_freq, CBoswell& log)
{
    cvm::cvector& y = fidproc.GetVectorFID(proc_coord);

    Task progress(log, "eddy current correction");
    cvm::cvector yW = fidwater.GetVectorFID(proc_coord);

    log.DebugMessage(DEBUG_LEVEL_1, "Size of y  in WaterEddy is: %d", y.size());
    log.DebugMessage(DEBUG_LEVEL_1, "Size of yW in WaterEddy is: %d", yW.size());

    assert( y.size() == yW.size() );

    // get phase of water signal
    cvm::rvector yphi( yW.size() );

    for( integer n = 1; n <= yW.size(); n++ ) 
        yphi[n] = atan2(imag(yW[n]), real(yW[n]));

    progress(".");

    // now subtract water phase from water suppressed signal phase but miss off last 20% (Siemens W data can have annoying ringing here that causes problems)
    for( integer n = 1; n <= 0.8*yW.size(); n++ ) 
        y[n] = y[n] * exp( -tcomplex(0, yphi[n]) );

    treal fs = fidwater.GetSamplingFrequency();
	treal dt = 1.0 / fs;

    // apply this frequency shift
	for( integer n = 0; n < y.size(); n++ ) 
		y[n+1] = exp(tcomplex(0, 2.0 * M_PI * water_freq * n * dt)) * y[n+1];

    progress(".");
}

// Do eddy current correction, as described in introduction of:
// "Artifacts Introduced by Zero Order Phase Correction in Proton
// NMR Spectroscopy and a Method of Elimination by Phase Filtering"
// by J. M. Wild, 1998, J. MR, V137, pp 430-436
// also, do custom alignment step
void tarquin::WaterEddy(const coord& proc_coord, CFID& fidproc, const CFID& fidwater, CBoswell& log)
{
	cvm::cvector& y = fidproc.GetVectorFID(proc_coord);
	
	// grab a copy for use in the alignment stage
	cvm::cvector yorig = y;

	{
		Task progress(log, "eddy current correction");
		const cvm::cvector& yW = fidwater.GetVectorFID(proc_coord);

		log.DebugMessage(DEBUG_LEVEL_1, "Size of y  in WaterEddy is: %d", y.size());
		log.DebugMessage(DEBUG_LEVEL_1, "Size of yW in WaterEddy is: %d", yW.size());

		assert( y.size() == yW.size() );
		//assert( options.GetRangeEnd() != 0 );

		// get phase of water signal
		cvm::rvector yphi( yW.size() );

		for( integer n = 1; n <= yW.size(); n++ ) 
			yphi[n] = atan2(imag(yW[n]), real(yW[n]));

		progress(".");

		integer nEnd = yW.size(); 

		// after the last point we are interested in, we keep phase the same
		// because otherwise it is just noise
		for( integer n = nEnd; n <= yW.size(); n++ ) 
			yphi[n] = yphi[nEnd];

		progress(".");

		// now subtract water phase from water suppressed signal phase
		for( integer n = 1; n <= yW.size(); n++ ) 
			y[n] = y[n] * exp( -tcomplex(0, yphi[n]) );

		progress(".");
	}

	Task progress(log, "realignment");

	// range over which we compute inner product
	treal range_ppm = 0.2;

	// number of points of evaluation
	integer num_steps = 100;

	treal fs = fidproc.GetSamplingFrequency();
	treal fr = fidproc.GetTransmitterFrequency();
	treal dt = 1.0 / fs;
	treal range_hz  = (range_ppm * fr) / 1e6;
	treal step_size = range_hz / num_steps;

	cvm::cvector yprime = y;

	// the stepsize at which the maximum occurr
	treal currentMax = -std::numeric_limits<treal>::infinity();
	treal fshiftMax  = 0;

	for( integer step = 0; step < num_steps; step++ ) 
	{
		treal fshift = -(range_hz/2) + step_size*step;

		// apply this frequency shift
		for( integer n = 0; n < y.size(); n++ ) 
			yprime[n+1] = exp(tcomplex(0, 2.0 * M_PI * fshift * n * dt)) * y[n+1];

		tcomplex result_cmp = 0;

		// compute inner product (FIXME: work out which cvmlib operator to use)
		for( integer n = 0; n < y.size(); n++ ) 
			result_cmp += yprime[n+1] * conj(yorig[n+1]);

		treal result = abs(result_cmp);

		// is this the biggest yet?
		if( result > currentMax ) 
		{
			currentMax = result;
			fshiftMax = fshift;
		}

		progress(".");
	}

	// apply this frequency shift
	for( integer n = 0; n < y.size(); n++ ) 
	{
		y[n+1] = exp(tcomplex(0, 2.0 * M_PI * fshiftMax * n * dt)) * y[n+1];
	}
}

// parameters structure for the optimiser used in phasing
struct SPhaseParams 
{
	SPhaseParams(cvm::cvector& Y_in, cvm::rvector& freq_range_in) : 
		Y(Y_in), freq_range(freq_range_in)
	{

	}

	cvm::cvector& Y;
	cvm::rvector& freq_range;
};

void tarquin::ApplyPhase(const coord& proc_coord, CFID& fid, CFID& av_fid, CBoswell& log)
{
	cvm::cvector& y = fid.GetVectorFID(proc_coord);
	//treal fs = fid.GetSamplingFrequency();

	cvm::cvector ynew = y;

	assert( y.size() > 0 );
	//assert( fs > 0 );

	Task progress(log, "phase parameter application");

	// number of points in the signal
	integer N = y.size();

	// make the DFT of the time signal
	cvm::cvector Y(N);
	fft(y, Y);

	Y = fftshift(Y);

	// get the shifted [-fs/2, ..., 0, ... fs/2) frequency range
	cvm::rvector freq_range = fid.GetFreqScale();

	// the parameters we are trying to find
	cvm::rvector phi(2);
    const coord start_spec(1, 1, 1); 
	phi(1) = av_fid.GetPhi0(start_spec);
	phi(2) = av_fid.GetPhi1(start_spec);
    
    //std::cout << phi(1) << std::endl;
    //std::cout << phi(2) << std::endl;

	// now apply the parameters we have just found to the outgoing signal
	for( integer n = 0; n < N; n++ ) 
		Y[n+1] = Y[n+1] * exp( tcomplex(0, phi(1) + phi(2)*2.0*M_PI*freq_range[n+1]) );

	Y = fftshift(Y);

	// go back to normal range
	ifft(Y, y);

	//fid.SetPhi0(phi[1]);
	//fid.SetPhi1(phi[2]);
}

void tarquin::ApplyPhase(const coord& proc_coord, CFID& fid, CBoswell& log)
{
	cvm::cvector& y = fid.GetVectorFID(proc_coord);
	//treal fs = fid.GetSamplingFrequency();

	cvm::cvector ynew = y;

	assert( y.size() > 0 );
	//assert( fs > 0 );

	Task progress(log, "phase parameter application");

	// number of points in the signal
	integer N = y.size();

	// make the DFT of the time signal
	cvm::cvector Y(N);
	fft(y, Y);

	Y = fftshift(Y);

	// get the shifted [-fs/2, ..., 0, ... fs/2) frequency range
	cvm::rvector freq_range = fid.GetFreqScale();

	// the parameters we are trying to find
	cvm::rvector phi(2);
	phi(1) = fid.GetPhi0(proc_coord);
	phi(2) = fid.GetPhi1(proc_coord);

	// now apply the parameters we have just found to the outgoing signal
	for( integer n = 0; n < N; n++ ) 
		Y[n+1] = Y[n+1] * exp( tcomplex(0, phi(1) + phi(2)*2.0*M_PI*freq_range[n+1]) );

	Y = fftshift(Y);

	// go back to normal range
	ifft(Y, y);

	//fid.SetPhi0(phi[1]);
	//fid.SetPhi1(phi[2]);
}

// Phase the spectrum by minimising the difference between the magnitude and real spectrum.
// This probably needs to do have different methods of doing things, depending on the echo time,
// because of upside-down peaks, etc.
void tarquin::AutoPhase(const coord& proc_coord, CFID& fid, const Options& opts, CBoswell& log, bool water_file)
{
	cvm::cvector& y = fid.GetVectorFID(proc_coord);
	treal fs = fid.GetSamplingFrequency();

	cvm::cvector ynew = y;

	assert( y.size() > 0 );
	assert( fs > 0 );

	// This doesn't always work very well.
	// I think a crude initial sweep of various values of phi0 to determine a good
	// starting position would help it.

    if ( water_file )
        Task progress(log, "autophasing mag method water");
    else
        Task progress(log, "autophasing mag method");

	// number of points in the signal
	integer N = y.size();

	// make the DFT of the time signal
	cvm::cvector Y(N);
	fft(y, Y);

	Y = fftshift(Y);

	// get the shifted [-fs/2, ..., 0, ... fs/2) frequency range
	cvm::rvector freq_range = fid.GetFreqScale();

	// the magnitude spectrum
	cvm::rvector YMAG(N);
	cvm::rvector YR = Y.real();
	cvm::rvector YI = Y.imag();

	for( integer n = 0; n < N; n++ )
		YMAG[n+1] = sqrt(YR[n+1]*YR[n+1] + YI[n+1]*YI[n+1]);

	// the parameters we are trying to find
	cvm::rvector phi(2);
	phi(1) = 0.0;
	phi(2) = 0.0;

	// lower limits for phi0 and phi1
	cvm::rvector vLowerBounds(2);
	vLowerBounds(1) = -M_PI;
	vLowerBounds(2) = -5 / fs;

	// upper limits for phi0 and phi1
	cvm::rvector vUpperBounds(2);
	vUpperBounds(1) = +M_PI;
	vUpperBounds(2) = +5 / fs;

	// parameters for the optimiser
	SPhaseParams params(Y, freq_range);

    cvm::rvector freq_scale = fid.GetPPMScale(proc_coord, 1);
    // find points corresponding to ppm start and ppm end
    int left = 1, right = freq_scale.size() - 1;
    for ( int n = 1; n < (freq_scale.size()); n++ ) 
    {
        if ( ( freq_scale(n) > opts.GetPPMend() ) && ( freq_scale(n+1) <= opts.GetPPMend() ) )
            left = n;
        if ( ( freq_scale(n) > opts.GetPPMstart() ) && ( freq_scale(n+1) <= opts.GetPPMstart() ) )
            right = n;
    }

    if ( left == 1 || right == (freq_scale.size() - 1) )
        log.LogMessage(LOG_INFO, "Warning, one of the PPM limits is outside the bandwidth.");


    if ( water_file )
    {
        left = 1;
        right = freq_scale.size();
    }

	// get a good starting value for phi0
	// sweep phi0 and find the lowest residual
	cvm::rvector diff(YMAG.size());
	cvm::rvector yp(YMAG.size());
	treal min_norm = std::numeric_limits<treal>::infinity();
	treal min_phi0 = 0.0;

	for( treal phi0 = -M_PI; phi0 < M_PI; phi0 += 2*M_PI/100.0 ) 
	{
		phi(1) = phi0;

		phasing_objective(phi, yp, N, &params);

		diff = YMAG - yp;
        treal val = diff.norm();
        
        /*treal val = 0;
	    for ( int n = left; n < right+1; n++ ) 
        {
            val += diff(n)*diff(n);
        }*/

		if( diff.norm() < min_norm )  
		{
			min_norm = val;
			min_phi0 = phi0;
		}
	}

	phi(1) = min_phi0;

	/*
	   double opts[LM_OPTS_SZ];
	   opts[0] = NUMERICAL_TOL_INIT_MU;
	   opts[1] = NUMERICAL_TOL_INF_JE;
	   opts[2] = NUMERICAL_TOL_L2_DP;
	   opts[3] = NUMERICAL_TOL_L2_RES;
	   opts[4] = NUMERICAL_TOL_DIFF_DELTA;

	// info about result
	double info[LM_INFO_SZ];

	// call the optimiser
	dlevmar_bc_dif(phasing_objective, phi, YMAG, 2, N, vLowerBounds, vUpperBounds, 100, 
	opts, info, NULL, NULL, (void*)&params);

*/

	/*	std::cout << "\nl2 norm of error at initial p  = " << info[0];
		std::cout << "\nl2 norm of error at final p    = " << info[1];
		std::cout << "\nl2 norm of J.'*e at final p    = " << info[2];
		std::cout << "\nl2 norm of D*p at final p      = " << info[3];
		std::cout << "\nnumber of iterations           = " << info[5];
		std::cout << "\nnumber of function evaluations = " << info[7];
		std::cout << "\nnumber of Jacobian evaluations = " << info[8];

		if( 1 == info[6] )
		std::cout << "\nstopped by small gradient";
		else if( 2 == info[6] )
		std::cout << "\nstopped by small D*p";
		else if( 3 == info[6] )
		std::cout << "\nstopped by iteration limit";
		else if( 4 == info[6] )
		std::cout << "\nstopped by singular matrix - restart with bigger mu";
		else if( 5 == info[6] )
		std::cout << "\nno further reduction possible - restart with bigger mu";
		else if( 6 == info[6] )
		std::cout << "\nstopping by small l2 norm of error"; */

	// now apply the parameters we have just found to the outgoing signal

	phi(2) = 0;
    
    for( integer n = 0; n < N; n++ ) 
	{
		Y[n+1] = Y[n+1] * exp( tcomplex(0, phi(1) + phi(2)*2.0*M_PI*freq_range[n+1]) );
	}

    // apply crude baseline correction
    int pts = static_cast<int>(N*0.05);
	float real_av = 0.0f;
	float cplx_av = 0.0f;
	for ( int n = 1; n <= pts; n++) 
	{
		real_av = static_cast<float>(real_av + Y(n).real()/pts);
		cplx_av = static_cast<float>(cplx_av + Y(n).imag()/pts);
	}

    // subtract mean and scale by number of scans 
    //plot(Y);
	//std::complex<double> av(real_av, cplx_av);
	//for ( int n = 1; n <= N; n++) 
	//	Y(n) = ( Y(n) - av );
    //plot(Y);


	/*for( integer n = 0; n < N; n++ ) 
	{
		Y[n+1] = Y[n+1] * exp( tcomplex(0, phi(1) + phi(2)*2.0*M_PI*freq_range[n+1]) );
	}*/

    

    treal maxY = -std::numeric_limits<treal>::infinity();
    treal minY = std::numeric_limits<treal>::infinity();

	//for( integer n = 0; n < N; n++ ) 
	for ( int n = left; n < right+1; n++ ) 
	{
		//Y[n+1] = Y[n+1] * exp( tcomplex(0, phi(1) + phi(2)*2.0*M_PI*freq_range[n+1]) );

        if ( Y[n].real() > maxY )
            maxY = Y[n].real();

        if ( Y[n].real() < minY )
            minY = Y[n].real();

	}
    
    //std::cout << std::endl << "min y " << minY << std::endl;
    //std::cout << "max y " << maxY << std::endl;
    //std::cout << "real_av " << real_av << std::endl;
    //plot(Y);
    
    /* NOT REALLY NEEDED for this phasing algo
    // change phase by 180 degrees if required
    if ( fabs(minY-real_av) > fabs(maxY-real_av) && ( opts.GetPulSeq() != MEGA_PRESS ) )
    {
	    Task progress(log, "flipping data");
        treal delta;
        if ( phi(1) >= 0 )
            delta = -M_PI;
        else 
            delta = +M_PI;

        phi(1) = phi(1) + delta;
        for( integer n = 0; n < N; n++ ) 
		    Y[n+1] = Y[n+1] * exp( tcomplex(0, delta + phi(2)*2.0*M_PI*freq_range[n+1]) );

    }
    
    // change phase by 180 degrees if required
    if ( fabs(minY-real_av) < fabs(maxY-real_av) && ( opts.GetPulSeq() == MEGA_PRESS ) )
    {
	    Task progress(log, "flipping data");
        treal delta;
        if ( phi(1) >= 0 )
            delta = -M_PI;
        else 
            delta = +M_PI;

        phi(1) = phi(1) + delta;
        for( integer n = 0; n < N; n++ ) 
		    Y[n+1] = Y[n+1] * exp( tcomplex(0, delta + phi(2)*2.0*M_PI*freq_range[n+1]) );

    }
    */

	Y = fftshift(Y);

	// go back to normal range
	ifft(Y, y);

	fid.SetPhi0(proc_coord, phi[1]);
	fid.SetPhi1(proc_coord, phi[2]);
}


void tarquin::AutoPhaseNew(const coord& proc_coord, CFID& fid,const Options& opts, CBoswell& log, bool water_file)
{
	cvm::cvector& y = fid.GetVectorFID(proc_coord);

    integer zf = 2;
    // apply zero filling
	integer N = y.size();

	cvm::cvector yzf = y;
    yzf.resize(N*zf);

	assert( y.size() > 0 );

    if ( water_file )
        Task progress(log, "autophasing water (new)");
    else
        Task progress(log, "autophasing (new)");

	// make the DFT of the time signal
	cvm::cvector Yzf(N*zf);
	fft(yzf, Yzf);
	Yzf = fftshift(Yzf);

	// get the shifted [-fs/2, ..., 0, ... fs/2) frequency range
	cvm::rvector freq_range_zf = fid.GetFreqScale(zf);

	// the parameters we are trying to find
	cvm::rvector phi(2);
	phi(1) = 0.0;
	phi(2) = 0.0;

    cvm::rvector freq_scale_zf = fid.GetPPMScale(proc_coord, zf);
    // find points corresponding to ppm start and ppm end
    int left_zf = 1, right_zf = freq_scale_zf.size() - 1;
    for ( int n = 1; n < (freq_scale_zf.size()); n++ ) 
    {
        if ( ( freq_scale_zf(n) > opts.GetPPMend() ) && ( freq_scale_zf(n+1) <= opts.GetPPMend() ) )
            left_zf = n;
        if ( ( freq_scale_zf(n) > opts.GetPPMstart() ) && ( freq_scale_zf(n+1) <= opts.GetPPMstart() ) )
            right_zf = n;
    }
    
    if ( left_zf == 1 || right_zf == (freq_scale_zf.size() - 1) )
        log.LogMessage(LOG_INFO, "Warning, one of the PPM limits is outside the bandwidth.");



    if ( water_file )
    {
        left_zf = 1;
        right_zf = freq_scale_zf.size();
    }

	// get a good starting value for phi0
	// sweep phi0 and find the lowest residual
	cvm::rvector yp(Yzf.size());
	cvm::rvector ypf(Yzf.size());
	treal max_norm = -std::numeric_limits<treal>::infinity();
	treal max_phi0 = 0.0;

	SPhaseParams params(Yzf, freq_range_zf);

	for( treal phi0 = -M_PI; phi0 < M_PI; phi0 += 2*M_PI/100.0 ) 
	{
		phi(1) = phi0;

		phasing_objective(phi, yp, yp.size(), &params);

        // filter
        td_conv_ws_fft(yp, ypf, 200, 16);
        ypf = yp - ypf;
        
        //plot(ypf);
        //ypf = yp;
        //std::cout << phi0 << std::endl;
        //std::cout << ypf.norminf() << std::endl;
        
        double Ymax = -std::numeric_limits<treal>::infinity();
        double Ysum = 0;
		// find max of y for SNR calculation
		for ( int n = left_zf; n < right_zf+1; n++ ) 
		{
            Ysum += ypf(n);
			if ( Ymax < ypf(n) )
				Ymax = ypf(n);
		}

        /*std::cout << "Ysum : " << Ysum << std::endl;
        std::cout << "Ymax : " << Ymax << std::endl;
        std::cout << "Yinf : " << ypf.norminf() << std::endl;
        */
        
        treal norm = 0;

        norm = Ymax;
        //norm = ypf.norminf();

		if( norm > max_norm )  
		{
			max_norm = norm;
			max_phi0 = phi0;
		}
	}

	phi(1) = max_phi0;
    

	// now apply the parameters we have just found to the outgoing signal

    log.DebugMessage(DEBUG_LEVEL_1, "Max Phi0 : %.2f", max_phi0);
	
    // make the DFT of the time signal
	cvm::cvector Y(N);
	fft(y, Y);
	Y = fftshift(Y);

	// get the shifted [-fs/2, ..., 0, ... fs/2) frequency range
	cvm::rvector freq_range = fid.GetFreqScale();

	phi(2) = 0;

    treal maxY = -std::numeric_limits<treal>::infinity();
    treal minY = std::numeric_limits<treal>::infinity();

	for( integer n = 0; n < N; n++ ) 
	{
		Y[n+1] = Y[n+1] * exp( tcomplex(0, phi(1) + phi(2)*2.0*M_PI*freq_range[n+1]) );
    }

    cvm::rvector freq_scale = fid.GetPPMScale(proc_coord);
    // find points corresponding to ppm start and ppm end
    int left = 1, right = freq_scale.size() - 1;
    for ( int n = 1; n < (freq_scale.size()); n++ ) 
    {
        if ( ( freq_scale(n) > opts.GetPPMend() ) && ( freq_scale(n+1) <= opts.GetPPMend() ) )
            left = n;
        if ( ( freq_scale(n) > opts.GetPPMstart() ) && ( freq_scale(n+1) <= opts.GetPPMstart() ) )
            right = n;
    }

    if ( left == 1 || right == (freq_scale.size() - 1) )
        log.LogMessage(LOG_INFO, "Warning, one of the PPM limits is outside the bandwidth.");

    
    if ( water_file )
    {
        left = 1;
        right = freq_scale.size();
    }


	for ( int n = left; n < right+1; n++ ) 
    {
        if ( Y[n].real() > maxY )
            maxY = Y[n].real();

        if ( Y[n].real() < minY )
            minY = Y[n].real();
    }

    /*std::cout << maxY << std::endl;
    std::cout << minY << std::endl;
    
    std::cout << fabs(maxY) << std::endl;
    std::cout << fabs(minY) << std::endl;*/

    // change phase by 180 degrees if required
    if ( fabs(minY) > fabs(maxY) && ( opts.GetPulSeq() != MEGA_PRESS ) )
    {
	    Task progress(log, "flipping data");
        treal delta;
        if ( phi(1) >= 0 )
            delta = -M_PI;
        else 
            delta = +M_PI;

        phi(1) = phi(1) + delta;
        for( integer n = 0; n < N; n++ ) 
		    Y[n+1] = Y[n+1] * exp( tcomplex(0, delta + phi(2)*2.0*M_PI*freq_range[n+1]) );

    }
    
    // change phase by 180 degrees if required
    if ( fabs(minY) < fabs(maxY) && ( opts.GetPulSeq() == MEGA_PRESS ) )
    {
	    Task progress(log, "flipping data");
        treal delta;
        if ( phi(1) >= 0 )
            delta = -M_PI;
        else 
            delta = +M_PI;

        phi(1) = phi(1) + delta;
        for( integer n = 0; n < N; n++ ) 
		    Y[n+1] = Y[n+1] * exp( tcomplex(0, delta + phi(2)*2.0*M_PI*freq_range[n+1]) );

    }

    //plot(Y);

	Y = fftshift(Y);
    
	// go back to normal range
	ifft(Y, y);

	fid.SetPhi0(proc_coord, phi[1]);
	fid.SetPhi1(proc_coord, phi[2]);
}

void tarquin::AutoPhaseSimple(const coord& proc_coord, CFID& fid, CBoswell& log)
{
	cvm::cvector& y = fid.GetVectorFID(proc_coord);
	treal fs = fid.GetSamplingFrequency();

	cvm::cvector ynew = y;

	assert( y.size() > 0 );
	assert( fs > 0 );

	// This doesn't always work very well.
	// I think a crude initial sweep of various values of phi0 to determine a good
	// starting position would help it.

	Task progress(log, "autophasing (simple)");

	// number of points in the signal
	integer N = y.size();

	// make the DFT of the time signal
	cvm::cvector Y(N);
	fft(y, Y);

	Y = fftshift(Y);

	// get the shifted [-fs/2, ..., 0, ... fs/2) frequency range
	cvm::rvector freq_range = fid.GetFreqScale();

	// the magnitude spectrum
	cvm::rvector YMAG(N);
	cvm::rvector YR = Y.real();
	cvm::rvector YI = Y.imag();

	for( integer n = 0; n < N; n++ )
		YMAG[n+1] = sqrt(YR[n+1]*YR[n+1] + YI[n+1]*YI[n+1]);

	// the parameters we are trying to find
	cvm::rvector phi(2);
	phi(1) = 0.0;
	phi(2) = 0.0;

	// lower limits for phi0 and phi1
	cvm::rvector vLowerBounds(2);
	vLowerBounds(1) = -M_PI;
	vLowerBounds(2) = -5 / fs;

	// upper limits for phi0 and phi1
	cvm::rvector vUpperBounds(2);
	vUpperBounds(1) = +M_PI;
	vUpperBounds(2) = +5 / fs;

	// parameters for the optimiser
	SPhaseParams params(Y, freq_range);

	// get a good starting value for phi0
	// sweep phi0 and find the lowest residual
	cvm::rvector yp(YMAG.size());
	cvm::rvector ypf(YMAG.size());
	treal max_norm = -std::numeric_limits<treal>::infinity();
	treal max_phi0 = 0.0;

	for( treal phi0 = -M_PI; phi0 < M_PI; phi0 += 2*M_PI/100.0 ) 
	{
		phi(1) = phi0;

		phasing_objective(phi, yp, N, &params);

        // filter
        //td_conv_ws(yp, ypf, 200, 16);
        td_conv_ws_fft(yp, ypf, 200, 16);
        ypf = yp - ypf;
        /*for ( int n = 1; n < N + 1; n++ )
        {
            if ( ypf(n) < 0 )
                ypf(n) = 0;
        }*/
        
        //plot(ypf);


        //ypf = yp;
        
        //std::cout << phi0 << std::endl;
        //std::cout << ypf.norminf() << std::endl;
        treal norm = 0;
	    //for( integer n = 1; n < N+1; n++ ) 
        //    norm += ypf(n);

        norm = ypf.norminf();

		if( norm > max_norm )  
		{
			max_norm = norm;
			max_phi0 = phi0;
		}

		/*if( ypf.norminf() > max_norm )  
		{
			max_norm = ypf.norminf();
			max_phi0 = phi0;
		}*/
	}

	phi(1) = max_phi0;
	//std::cout << max_phi0 << std::endl;

	// now apply the parameters we have just found to the outgoing signal

	phi(2) = 0;

    treal maxY = -std::numeric_limits<treal>::infinity();
    treal minY = std::numeric_limits<treal>::infinity();

	for( integer n = 0; n < N; n++ ) 
	{
		Y[n+1] = Y[n+1] * exp( tcomplex(0, phi(1) + phi(2)*2.0*M_PI*freq_range[n+1]) );

        if ( Y[n+1].real() > maxY )
            maxY = Y[n+1].real();

        if ( Y[n+1].real() < minY )
            minY = Y[n+1].real();
	}

    /*std::cout << maxY << std::endl;
    std::cout << minY << std::endl;
    
    std::cout << fabs(maxY) << std::endl;
    std::cout << fabs(minY) << std::endl;*/

    // change phase by 180 degrees if required
    if ( fabs(minY) > fabs(maxY) )
    {
        treal delta;
        if ( phi(1) >= 0 )
            delta = -M_PI;
        else 
            delta = +M_PI;

        phi(1) = phi(1) + delta;
        for( integer n = 0; n < N; n++ ) 
		    Y[n+1] = Y[n+1] * exp( tcomplex(0, delta + phi(2)*2.0*M_PI*freq_range[n+1]) );

    }

    //plot(Y);

	Y = fftshift(Y);
    
	// go back to normal range
	ifft(Y, y);

	fid.SetPhi0(proc_coord, phi[1]);
	fid.SetPhi1(proc_coord, phi[2]);
}

// optimiser will attempt to minimise norm(YP - YMAG, 2)
void tarquin::phasing_objective(treal* pp, treal* pYP, integer N, void* pParams)
{
	SPhaseParams& params = *static_cast<SPhaseParams*>(pParams);

	// the current candidate set of variables
	cvm::rvector phi(pp, 2);

	// the output of this function, our FD estimate of the phased signal
	cvm::rvector YP(pYP, N);

	// the signal we are going to modify to produce YP
	cvm::cvector& Y = params.Y;

	// the frequency range over which the parameters are applied
	cvm::rvector& freq_range = params.freq_range;

	// apply the candidate parameters
	for( integer n = 0; n < N; n++ ) 
	{
		tcomplex z = Y[n+1] * exp( tcomplex(0, phi(1) + phi(2)*2.0*M_PI*freq_range[n+1]) );

		YP[n+1] = z.real();
	}
}

namespace tarquin
{
struct PeakRef 
{
	PeakRef()
	{

	}

	PeakRef(treal fppml, treal fppmh, treal rppm) : 
		freq_ppm_lower(fppml), freq_ppm_higher(fppmh), ref_ppm(rppm)
	{

	}

	treal freq_ppm_lower;
	treal freq_ppm_higher;
	treal ref_ppm;
};
}

namespace tarquin
{
struct RefPoint
{
	RefPoint()
	{

	}

	RefPoint(treal center_ppm, treal ref_amp, treal doub_split_hz, int metab_lip) : 
		ppm(center_ppm), amp(ref_amp), split_hz(doub_split_hz), ml(metab_lip)
	{

	}

	treal ppm;
	treal amp;
	treal split_hz;
	int ml;
};
}

//! Return a reference in ppm.
bool tarquin::AutoReferenceCorr(const coord& proc_coord, Options& options, CFID& fid, bool dyn_mode, bool pre_ws_shift, bool skip_beta_guess, CBoswell& bos_log)
{
    // this function calculates the cross correlation
    // of the spectra with a series of reference peaks to
    // find the optimal ref and an estimate of the inital
    // beta variable

    //std::cout << std::endl << options.GetMaxDRef() << std::endl;
    //std::cout << options.GetRefFile() << std::endl;

    std::vector<RefPoint> points;

    if ( !dyn_mode )
    {
        // if the ref file string is empty assume defualt values	
        if ( options.GetRefFile() == "" )
        {
            if ( options.GetRefSignals() == PROTON_NAA_CR_CHO_LIP )
            {
                points.push_back( RefPoint(2.01, 1, 0, 1) ); // NAA
                points.push_back( RefPoint(3.22, 1, 0, 1) ); // TCho
                points.push_back( RefPoint(3.03, 1, 0, 1) ); // Cr
                points.push_back( RefPoint(1.28, 1, 0, 0) ); // Lip
                points.push_back( RefPoint(0.9, 1, 0, 0) ); // Lip
                if ( pre_ws_shift )
                    points.push_back( RefPoint(4.65, 1, 0, 1) ); // H2O
            }
            else if ( options.GetRefSignals() == PROTON_NAA_CR_CHO )
            {
                points.push_back( RefPoint(2.01, 1, 0, 1) ); // NAA
                points.push_back( RefPoint(3.22, 1, 0, 1) ); // TCho
                points.push_back( RefPoint(3.03, 1, 0, 1) ); // Cr
                if ( pre_ws_shift )
                    points.push_back( RefPoint(4.65, 1, 0, 1) ); // H2O
            }
            else if ( options.GetRefSignals() == PROTON_CR_CHO )
            {
                points.push_back( RefPoint(3.22, 1, 0, 1) ); // TCho
                points.push_back( RefPoint(3.03, 1, 0, 1) ); // Cr
                if ( pre_ws_shift )
                    points.push_back( RefPoint(4.65, 1, 0, 1) ); // H2O
            }
            else if ( options.GetRefSignals() == PROTON_NAA_CHO )
            {
                points.push_back( RefPoint(2.01, 1, 0, 1) ); // NAA
                points.push_back( RefPoint(3.22, 1, 0, 1) ); // TCho
                if ( pre_ws_shift )
                    points.push_back( RefPoint(4.65, 1, 0, 1) ); // H2O
            }
            else if ( options.GetRefSignals() == PROTON_NAA )
            {
                points.push_back( RefPoint(2.01, 1, 0, 1) ); // NAA
                if ( pre_ws_shift )
                    points.push_back( RefPoint(4.65, 1, 0, 1) ); // H2O
            }
            else if ( options.GetRefSignals() == PROTON_CR )
            {
                points.push_back( RefPoint(3.03, 1, 0, 1) ); // CR 
                if ( pre_ws_shift )
                    points.push_back( RefPoint(4.65, 1, 0, 1) ); // H2O
            }
            else if ( options.GetRefSignals() == PROTON_CHO )
            {
                points.push_back( RefPoint(3.22, 1, 0, 1) ); // TCho 
                if ( pre_ws_shift )
                    points.push_back( RefPoint(4.65, 1, 0, 1) ); // H2O
            }
            else if ( options.GetRefSignals() == PROTON_H2O )
            {
                points.push_back( RefPoint(4.65, 1, 0, 1) ); // H2O
            }
            else if ( options.GetRefSignals() == PROTON_LIP )
            {
                points.push_back( RefPoint(1.28, 1, 0, 0) ); // Lip
                if ( pre_ws_shift )
                    points.push_back( RefPoint(4.65, 1, 0, 1) ); // H2O
            }
            else if ( options.GetRefSignals() == PHOSPH_PCR )
            {
                points.push_back( RefPoint(0, 1, 0, 1) ); // PCR
            }
            else if ( options.GetRefSignals() == PHOSPH_PCR_GAMMAATP )
            {
                points.push_back( RefPoint(0, 1, 0, 1) ); // PCR
                points.push_back( RefPoint(-2.67, 1, 16.1, 1) ); // GAMMAAPT
            }
            else
            {
                bos_log.LogMessage(LOG_ERROR, "Error Ref Mode not found.");
                return false;
            }
        }
        else
        {
            // otherwise read from file
            CCSVFile file;
            if( false == file.load(options.GetRefFile() ) )
                return false;

            std::vector<std::vector<double> >& ref_mat = file.getDoubleMatrix();

            for( std::vector<std::vector<double> >::iterator itP = ref_mat.begin(); itP != ref_mat.end(); itP++ ) 
            {
                //std::cout << (*itP)[0] << std::endl;
                //std::cout << (*itP)[1] << std::endl;
                //std::cout << (*itP)[2] << std::endl;
                points.push_back( RefPoint( (*itP)[0], (*itP)[1], (*itP)[2], (*itP)[3] ) );
            }
        }
    }
    else
    {
        if ( options.GetDynRefSignals() == PROTON_NAA_CR_CHO_LIP )
        {
            points.push_back( RefPoint(2.01, 1, 0, 1) ); // NAA
            points.push_back( RefPoint(3.22, 1, 0, 1) ); // TCho
            points.push_back( RefPoint(3.03, 1, 0, 1) ); // Cr
            points.push_back( RefPoint(1.28, 1, 0, 0) ); // Lip
            points.push_back( RefPoint(0.9, 1, 0, 0) ); // Lip
        }
        else if ( options.GetDynRefSignals() == PROTON_NAA_CR_CHO )
        {
            points.push_back( RefPoint(2.01, 1, 0, 1) ); // NAA
            points.push_back( RefPoint(3.22, 1, 0, 1) ); // TCho
            points.push_back( RefPoint(3.03, 1, 0, 1) ); // Cr
        }
        else if ( options.GetDynRefSignals() == PROTON_CR_CHO )
        {
            points.push_back( RefPoint(3.22, 1, 0, 1) ); // TCho
            points.push_back( RefPoint(3.03, 1, 0, 1) ); // Cr
        }
        else if ( options.GetDynRefSignals() == PROTON_NAA_CHO )
        {
            points.push_back( RefPoint(2.01, 1, 0, 1) ); // NAA
            points.push_back( RefPoint(3.22, 1, 0, 1) ); // TCho
        }
        else if ( options.GetDynRefSignals() == PROTON_NAA )
        {
            points.push_back( RefPoint(2.01, 1, 0, 1) ); // NAA
        }
        else if ( options.GetDynRefSignals() == PROTON_CR )
        {
            points.push_back( RefPoint(3.03, 1, 0, 1) ); // CR 
        }
        else if ( options.GetDynRefSignals() == PROTON_CHO )
        {
            points.push_back( RefPoint(3.22, 1, 0, 1) ); // TCho 
        }
        else if ( options.GetDynRefSignals() == PROTON_H2O )
        {
            points.push_back( RefPoint(4.65, 1, 0, 1) ); // H2O
        }
        else if ( options.GetDynRefSignals() == PHOSPH_PCR )
        {
            points.push_back( RefPoint(0, 1, 0, 1) ); // PCR
        }
        else if ( options.GetDynRefSignals() == PHOSPH_PCR_GAMMAATP )
        {
            points.push_back( RefPoint(0, 1, 0, 1) ); // PCR
            points.push_back( RefPoint(-2.67, 1, 16.1, 1) ); // GAMMAAPT
        }
        else
        {
            bos_log.LogMessage(LOG_ERROR, "Error Ref Mode not found.");
            return false;
        }
    }

    // maximum shift from the starting position allowed (ppm)
    treal max_shift = options.GetMaxDRef();

    // read in the fid	
	cvm::cvector y = fid.GetVectorFID(proc_coord);
    integer zf = 8;
    // apply zero filling
	integer N = y.size() * zf;
    y.resize(N);
	treal fs = fid.GetSamplingFrequency();
	treal ft = fid.GetTransmitterFrequency();
	// DFT is at this frequency resolution
	treal df = fs / N;
    // starting point for the ref
	treal ref_guess = fid.GetPPMRef(proc_coord);
	bos_log.LogMessage(LOG_INFO, "Ref start = %f ppm", ref_guess);

	cvm::cvector Y(N);
	fft(y, Y);
	Y = fftshift(Y);

	// get the magnitude spectrum
	cvm::cvector magY(N);
    for ( int n=1; n < N+1; n++ )
    {
        magY(n) = pow(pow(Y(n).real(),2) + pow(Y(n).imag(),2),0.5);
        //magY(n) = Y(n).real();
    }

    //plot(Y);
    //plot(magY);

    // generate a reference signal for all reference points
	cvm::cvector ref_sig(N);
	std::complex<double> j(0,1); 
    ref_sig.set(0);
	for( std::vector<RefPoint>::iterator itP = points.begin(); itP != points.end(); itP++ ) 
    {
        for (int n = 1; n < N+1; n++)
        {
            // convert ppm to Hz
            treal c_freq = (ref_guess-itP->ppm)*ft/pow(10,6.0);

	        std::complex<double> jwp(0,2*M_PI*(c_freq+itP->split_hz)*(n-1)/fs); 
	        std::complex<double> jwm(0,2*M_PI*(c_freq-itP->split_hz)*(n-1)/fs); 
            ref_sig(n) += itP->amp/2*exp(jwp);
            ref_sig(n) += itP->amp/2*exp(jwm);
        }
    }

    // generate a reference signal for just metabolites
	cvm::cvector ref_sig_metabs(N);
    ref_sig_metabs.set(0);
	for( std::vector<RefPoint>::iterator itP = points.begin(); itP != points.end(); itP++ ) 
    {
        // meet some condition to be a metabolite
        if ( itP->ml == 1 )
        {
            for (int n = 1; n < N+1; n++)
            {
                // convert ppm to Hz
                treal c_freq = (ref_guess-itP->ppm)*ft/pow(10,6.0);

                std::complex<double> jwp(0,2*M_PI*(c_freq+itP->split_hz)*(n-1)/fs); 
                std::complex<double> jwm(0,2*M_PI*(c_freq-itP->split_hz)*(n-1)/fs); 
                ref_sig_metabs(n) += itP->amp/2*exp(jwp);
                ref_sig_metabs(n) += itP->amp/2*exp(jwm);
            }
        }
    }

    // do the fft, fftshift and find the magnitude
    cvm::cvector ref_spec(N);
	fft(ref_sig, ref_spec);
	ref_spec = fftshift(ref_spec);

    cvm::cvector ref_spec_metabs(N);
	fft(ref_sig_metabs, ref_spec_metabs);
	ref_spec_metabs = fftshift(ref_spec_metabs);

	// get the magnitude spectrum
	cvm::cvector mag_ref_spec(N);
    for ( int n=1; n < N+1; n++ )
        mag_ref_spec(n) = pow(pow(ref_spec(n).real(),2) + pow(ref_spec(n).imag(),2),0.5);
    
    cvm::cvector mag_ref_spec_metabs(N);
    for ( int n=1; n < N+1; n++ )
        mag_ref_spec_metabs(n) = pow(pow(ref_spec_metabs(n).real(),2) + pow(ref_spec_metabs(n).imag(),2),0.5);
    
	//plot(mag_ref_spec);

    magY.resize(2*N);
    mag_ref_spec.resize(2*N);
    mag_ref_spec_metabs.resize(2*N);
    cvm::cvector ft_magY(2*N);
    cvm::cvector ft_mag_ref_spec(2*N);
    cvm::cvector ft_mag_ref_spec_metabs(2*N);
   
    // do the fft 
    fft(magY, ft_magY);
    fft(mag_ref_spec, ft_mag_ref_spec);
    // changed 26 Jun 12 from old line, possilbe bug?
    fft(mag_ref_spec_metabs, ft_mag_ref_spec_metabs);

    // old code
    // fft(mag_ref_spec, ft_mag_ref_spec_metabs);

    // conjugate the target
    ft_magY = ft_magY.conj();

    // multiply
    cvm::cvector prod(2*N);
    cvm::cvector prod_metabs(2*N);
    for ( int n = 1; n < 2*N+1; n++ )
    {
        prod(n) = ft_magY(n)*ft_mag_ref_spec(n);
        prod_metabs(n) = ft_magY(n)*ft_mag_ref_spec_metabs(n);
    }

	//prod(n) = ft_mag_ref_spec(n)*ft_mag_ref_spec(n);

    // do the ifft
    ifft(prod, prod);
    ifft(prod_metabs, prod_metabs);

    // fftshift
    cvm::cvector prod_shift(2*N);
    prod_shift = fftshift(prod);
    cvm::cvector prod_shift_metabs(2*N);
    prod_shift_metabs = fftshift(prod_metabs);
    
    // create an index
    cvm::rvector shift_index(2*N);
    for ( int n = 1; n < 2*N+1; n++ )
        shift_index(n) = -N+n-1;

    //std::cout << std::endl << max_shift << "Max PPM" << std::endl;
    //std::cout << max_shift*ft*(1e-6) << "Max Hz" << std::endl;
    //std::cout << static_cast<int>(max_shift*ft*(1e-6)/df) << "Max points" << std::endl;

    // trim the shift_index and prod_shift to the max shift allowed
    int max_pts = static_cast<int>(max_shift*ft*(1e-6)/df);
    cvm::cvector prod_shift_cut(max_pts*2+1);
    cvm::cvector prod_shift_cut_metabs(max_pts*2+1);
    cvm::rvector shift_index_cut(max_pts*2+1);
    for ( int n = 1; n < max_pts*2+2; n++ )
    {
        prod_shift_cut(n) = prod_shift(N-max_pts+n);
        prod_shift_cut_metabs(n) = prod_shift_metabs(N-max_pts+n);
        shift_index_cut(n) = shift_index(N-max_pts+n);
    }
    
    // and plot for debugging
    //plot(shift_index, prod_shift_metabs);
    //plot(shift_index_cut, prod_shift_cut);
    //plot(prod_shift_cut);
    //plot(shift_index_cut, prod_shift_cut);
    //plot(shift_index_cut, prod_shift_cut_metabs);
    //plot(shift_index_cut, prod_shift_cut);

    int shift_pts = shift_index_cut(prod_shift_cut.indofmax());

	// now compute reference for the biggest peak in the set
	treal ref_freq_hz = df * (shift_pts);
	treal ref_freq_ppm = ref_guess - (1e6 * ref_freq_hz / ft);

	// transfer this new ref to the FID
    if ( options.GetAutoReference() || pre_ws_shift )
    {
        fid.SetPPMRef(proc_coord, ref_freq_ppm);	
        bos_log.LogMessage(LOG_INFO, "Ref new = %f ppm", ref_freq_ppm);
    }

    if ( dyn_mode )
        return true;

    if ( skip_beta_guess )
        return true;
    
    int metab_cut_max_pt = prod_shift_cut_metabs.indofmax();
    int metab_max_pt = metab_cut_max_pt - max_pts + N;
    
    /*
	bos_log.LogMessage(LOG_INFO, "Max ind 1 = %i", prod_shift_cut_metabs.indofmax());
	bos_log.LogMessage(LOG_INFO, "Max ind 2 = %i", prod_shift_metabs.indofmax());
	bos_log.LogMessage(LOG_INFO, "Max ind 3 = %i", metab_max_pt);
	bos_log.LogMessage(LOG_INFO, "N = %i", N);
	bos_log.LogMessage(LOG_INFO, "Max pts = %i", max_pts);
    */

    //plot(shift_index, prod_shift_metabs);
    //plot(prod_shift_metabs);
    
    // work on smoothed data to find the base of the peak
	cvm::rvector prod_shift_metabs_smo(prod_shift_metabs.size());
	td_conv_ws(prod_shift_metabs.real(), prod_shift_metabs_smo, 9, 10);

    // guess init beta based on the fwhm of prod_shift_cut_metabs
    // find the right base of the peak
    int right_pt = -1;
    for ( int n = metab_max_pt + 2; n < prod_shift_metabs.size(); n++ )
    {
        if ( prod_shift_metabs_smo(n) > prod_shift_metabs_smo(n-1) )
        {
            right_pt = n-1;
            break;
        }
        if ( n < max_pts*2+1 )
            right_pt = n-1;
    }

	//bos_log.LogMessage(LOG_INFO, "Right pt = %i", right_pt);
    
    // find the left base of the peak
    int left_pt = -1;
    for ( int n = metab_max_pt-2; n > 0; n-- )
    {
        if ( prod_shift_metabs_smo(n) > prod_shift_metabs_smo(n+1) )
        {
            left_pt = n+1;
            break;
        }
        if ( n == 1 )
            left_pt = n+1;
    }

	//bos_log.LogMessage(LOG_INFO, "Left pt = %i", left_pt);

    
    if ( left_pt == -1 || right_pt == -1 )
    {
        bos_log.DebugMessage(DEBUG_LEVEL_1, "Warning, init beta calc failed, carry on regardless.");
        options.AppendInitBetaUsed(options.GetInitBeta());
        return true;
    }

    //plot(prod_shift_metabs);
    
    // calculate the base level of the peak
    treal base = ( prod_shift_metabs(left_pt).real() + prod_shift_metabs(right_pt).real() ) / 2.0;
	//bos_log.LogMessage(LOG_INFO, "base = %f", base);
    //treal base = 0;
    treal top = prod_shift_metabs( metab_max_pt ).real();
	//bos_log.LogMessage(LOG_INFO, "top = %f", top);
    
    //plot(shift_index_cut, prod_shift_cut_metabs);
    //plot(shift_index, prod_shift_metabs);

    // find half height
    treal hh = ( base + top ) / 2;

    // find the right hh point
    treal right_hh_pt = -1;
    for ( int n = metab_max_pt; n < prod_shift_metabs.size(); n++ )
    {
        if ( prod_shift_metabs(n).real() >= hh && prod_shift_metabs(n+1).real() < hh )
        {
            treal m = ( prod_shift_metabs(n).real() - prod_shift_metabs(n+1).real() ) / ( n - (n+1) );
            treal c = prod_shift_metabs(n).real() - m * n;
            right_hh_pt = (hh - c) / m;
            break;
        }
    }

    // find the left hh point
    treal left_hh_pt = -1;
    for ( int n = metab_max_pt; n > 1; n-- )
    {
        if ( prod_shift_metabs(n).real() >= hh && prod_shift_metabs(n-1).real() < hh )
        {
            treal m = ( prod_shift_metabs(n).real() - prod_shift_metabs(n-1).real() ) / ( n - (n-1) );
            treal c = prod_shift_metabs(n).real() - m * n;
            left_hh_pt = (hh - c) / m;
            break;
        }
    }

    /* 
    bos_log.DebugMessage(DEBUG_LEVEL_1, "HH val        : %g", hh);
    bos_log.DebugMessage(DEBUG_LEVEL_1, "Max pt        : %i", prod_shift_cut_metabs.indofmax());
    bos_log.DebugMessage(DEBUG_LEVEL_1, "Left HH pt    : %.2f", left_hh_pt);
    bos_log.DebugMessage(DEBUG_LEVEL_1, "Right HH pt   : %.2f", right_hh_pt);
    */

    treal width_hz = ( right_hh_pt - left_hh_pt ) * fs / N * 2;

    // width of mag mode FWHM is 3 times bigger than absorbtion mode 
    // so correct for this
    
    width_hz = width_hz / 3;

    // subtract 1.0Hz width_hz if needed
    treal init_beta;
    treal hz_correction = 1.0;
    if ( (width_hz - hz_correction) > 0 )
        init_beta = pow((width_hz - hz_correction)*M_PI/2, 2)/log(0.5);
    else
        init_beta = -1;
    
    // cap the value at 200 for poor spectra
    if ( init_beta < -200 )
        init_beta = -200;

    // std::cout << std::endl << width_hz << std::endl;
    bos_log.DebugMessage(DEBUG_LEVEL_1, "Width pts     : %.2f", right_hh_pt - left_hh_pt);
    bos_log.DebugMessage(DEBUG_LEVEL_1, "Width Hz est  : %.2f", width_hz);
    bos_log.DebugMessage(DEBUG_LEVEL_1, "Init beta est : %.2f", -init_beta); 

    //plot(shift_index_cut, prod_shift_cut_metabs);
    //plot(prod_shift_metabs);
    //plot(shift_index, prod_shift_metabs);
    //plot(prod_shift_metabs_smo);

    // only update init beta if the options variable hasn't been modified from 0
    // probally need a better way of doing this as a starting value of 0 should
    // be possible for auto referencing, one soln is to add an auto_beta option
    if ( options.GetInitBeta() == 0 )
    {
        if ( init_beta < 0 )
            options.AppendInitBetaUsed(-init_beta);
        else
            options.AppendInitBetaUsed(0);
    }
    else
    {
        options.AppendInitBetaUsed(options.GetInitBeta());
    }

    //options.SetInitBeta(100);

	return true;
}


//! Return a reference in ppm.
bool tarquin::AutoReference(const coord& proc_coord, CFID& fid, CBoswell& log)
{
	std::vector<PeakRef> peaks;

	// the peaks we will try and reference to (biggest peak is that which is referenced to)
	peaks.push_back( PeakRef(1.91, 2.11, 2.01) ); // NAA
	peaks.push_back( PeakRef( 3.12, 3.32, 3.22) ); // PC
	peaks.push_back( PeakRef(2.93, 3.13, 3.03) ); // Creatine

	const cvm::cvector& y = fid.GetVectorFID(proc_coord);
	integer N = y.size();
	treal fs = fid.GetSamplingFrequency();
	treal ft = fid.GetTransmitterFrequency();

	// bodge alert!!
	treal ref_guess = fid.GetPPMRef(proc_coord);
	log.LogMessage(LOG_INFO, "Ref start = %f ppm", ref_guess);


	cvm::cvector Y(N);
	fft(y, Y);
	Y = fftshift(Y);
	// get the real part of the spectrum
	cvm::rvector spec = Y.real();

	// DFT is at this frequency resolution
	treal df = fs / N;

	treal max_height_global = -std::numeric_limits<treal>::infinity();
	integer nMaxGlobal = 0;
	treal ref_ppm = 0;

	// for each possible peak to reference to...
	for( std::vector<PeakRef>::iterator itP = peaks.begin(); itP != peaks.end(); itP++ ) 
	{
		// convert ppm to frequency
		treal fstart = (ref_guess-itP->freq_ppm_higher)*ft/pow(10.0,6.0);
		treal fend = (ref_guess-itP->freq_ppm_lower)*ft/pow(10.0,6.0);
		// convert frequency to a point in the fft (how are we rounding here?)
		integer nRegionEnd = 1 + integer((fend + fs/2.0)/(fs/N));
		integer nRegionStart = 1 + integer((fstart + fs/2.0)/(fs/N));

		assert( nRegionEnd >= nRegionStart );

		// find biggest peak in this region 
		// set max_height_local to be -inf
		treal max_height_local = -std::numeric_limits<treal>::infinity();
		integer nMaxLocal = 0;

		for( integer n = nRegionStart; n <= nRegionEnd; n++ ) 
		{
			if( spec[n] > max_height_local ) 
			{
				max_height_local = spec[n];
				nMaxLocal = n;
			}
		}

		assert( 0 != nMaxLocal );

		// is this the biggest peak we have found so far?
		if( max_height_local > max_height_global ) 
		{
			max_height_global = max_height_local;
			nMaxGlobal = nMaxLocal;
			ref_ppm = itP->ref_ppm;
		}
	}

	// commented by MW, in odd cases ref_ppm could be both 0 and valid
	//assert( ref_ppm != 0 );
	assert( nMaxGlobal != 0);

	// now compute reference for the biggest peak in the set
	treal ref_freq_hz = -fs/2.0 + df * (nMaxGlobal-1);
	treal ref_freq_ppm = ref_ppm + (1e6 * ref_freq_hz / ft);

	// transfer this data to the FID
	fid.SetPPMRef(proc_coord, ref_freq_ppm);	
	log.LogMessage(LOG_INFO, "Ref new = %f ppm", ref_freq_ppm);

	return true;
}

namespace tarquin
{
	struct SResRefParams 
	{
		SResRefParams(
				integer nStart, 
				treal fs, 
				const cvm::cvector& y,
				cvm::rvector& shift_freq_range, 
				CBoswell& log
				) : 
			m_dt(1.0 / fs),
			m_y(y), 
			m_nStart(nStart),
			m_shift_freq_range(shift_freq_range), 
			m_log(log)
		{
			m_a.resize(1);
		}

		//! Time step.
		treal m_dt;

		//! The original whole FID vector.
		const cvm::cvector& m_y;

		//! The vector we are fitting to.
		cvm::rvector m_yActive;

		//! The offset in number of time samples.
		integer m_nStart;

		//! The basis vector.
		cvm::cvector m_q;

		//! The DFT of the basis vector.
		cvm::cvector m_Q;

		//! Subvector over which we are fitting [real; imag].
		cvm::rvector m_qsub;

		//! Frequency vector in the range [-fs/2, ..., 0, ... fs/2).
		cvm::rvector& m_shift_freq_range;

		//! Amplitude estimate (complex and size 1).
		cvm::cvector m_a;

		//! The logging object.
		CBoswell& m_log;
	};
}

/*
 * This is the objective function for fitting a single peak. 
 * This is called by the optimiser repeatedly to generate f(p), 
 * the levmar code will minimise ||x-f(p)||, where x is the observed signal, in this case
 * the water reference FID.
 *
 * Note: amplitude and phi0 are worked out implicitly in amplitude estimate.
 */
void tarquin::residual_objective_ref(treal* pp, treal* pyhat, integer nParams, integer activeNdbl, void* pParams)
{
	// parameters other than those being optimised over
	SResRefParams& params = *static_cast<SResRefParams*>(pParams);

	// some variables we use locally
	integer nStart     = params.m_nStart;
	integer nActive    = params.m_yActive.size() / 2;
	CBoswell& log      = params.m_log;
	treal dt           = params.m_dt;
	cvm::cvector& q    = params.m_q;
	cvm::cvector& Q    = params.m_Q;
	//cvm::rvector& qsub = params.m_qsub;
	cvm::cvector& a    = params.m_a;
	//cvm::rvector& shift_freq_range = params.m_shift_freq_range;

	// vector of parameters (shares memory)
	assert( 4 == nParams );
	cvm::rvector vParams(pp, nParams);

	// parameter we are determining
	treal freq  = vParams(1);
	treal alpha = vParams(2);
	treal beta  = vParams(3);
	treal phi1  = vParams(4);

	tcomplex jomega(0, 2.0 * M_PI * freq);

	// make the model with these parameters (all except phi1)
	for( integer n = 0; n < q.size(); n++ ) 
	{
		// time starts at the offset because of nStart
		treal t = dt*n;

		// the modified sample of the model signal
		q(n+1) = exp(-t*alpha -t*t*beta + t*jomega);
	}

	// take DFT of model signal
	fft(q, Q);

	// apply phi1 parameter
	for( integer n = 0; n < Q.size(); n++ ) 
	{
		treal freq = n / (Q.size() * dt);

		Q[n+1]  = Q[n+1] * exp( tcomplex(0, phi1*2.0*M_PI*freq) );
	}

	ifft(Q, q);

	cvm::cmatrix A( activeNdbl / 2, 1);
	cvm::cvector y( activeNdbl / 2 );

	for( integer n = 0; n < A.msize(); n++ ) 
	{
		A(n+1, 1) = q(n+nStart);
		y(n+1)    = params.m_y(n+nStart);
	}


	cvm::cmatrix X = A.pinv(NUMERICAL_TOL_NNLS);
	a = X*y;

	assert( a.size() == 1 );
	//assert( qsub.size() == params.m_yActive.size() );
	assert( activeNdbl == params.m_yActive.size() );

	// estimate signal - yhat (shares memory with pyhat provided by levmar)
	cvm::rvector yhat( pyhat, activeNdbl );

	//cvm::cmatrix f(yhat.size()/2, 2);

	// construct [real; imag] result vector
	for( integer n = 0; n < nActive; n++ ) 
	{
		tcomplex z = q(n+nStart) * a[1];

		yhat(n+1)         = z.real();
		yhat(n+1+nActive) = z.imag();

		//   f(n+1, 1) = z;
		//  f(n+1, 2) = params.m_y(n+nStart);
	}

	log.UpdateTask(".");

	//cvm::cmatrix F(f.msize(), f.nsize());
	//F = fft(f);

	//cvm::cmatrix Z = fftshift(F);

	//plot(Z);
}

tarquin::treal tarquin::ComputeWaterNormalisation(const coord& proc_coord, const CFID& fid, CBoswell& log)
{
	// some local copies for quick ref.
	const cvm::cvector& y = fid.GetVectorFID(proc_coord);

	cvm::cvector yy = y;
    //std::cout << std::endl << yy.size() << std::endl;
	//plot(yy);

    // fitting method from mathworld
    // http://mathworld.wolfram.com/LeastSquaresFittingExponential.html

	cvm::rvector t = fid.GetTimeScale();

	int start_pt = 10;
	int end_pt = 50;

	cvm::rvector x_sec(end_pt-start_pt+1);
	cvm::rvector x_sec_sq(end_pt-start_pt+1);
	cvm::rvector y_sec(end_pt-start_pt+1);
	cvm::rvector y_sec_sq(end_pt-start_pt+1);
	cvm::rvector y_sec_ln(end_pt-start_pt+1);

	for(int n = start_pt; n<= end_pt; n++) 
	{
		x_sec(n-start_pt+1) = t(n);
		x_sec_sq(n-start_pt+1) = t(n)*t(n);
		y_sec(n-start_pt+1) = pow(pow(y(n).real(),2) + pow(y(n).imag(),2),0.5);
		y_sec_sq(n-start_pt+1) = pow(pow(y(n).real(),2) + pow(y(n).imag(),2),0.5)*pow(pow(y(n).real(),2) + pow(y(n).imag(),2),0.5);
		y_sec_ln(n-start_pt+1) = std::log(y_sec(n-start_pt+1));
	}
    
	treal a;
	treal sum1 = 0;
	treal sum2 = 0;
	treal sum3 = 0;
	treal sum4 = 0;
	treal sum5 = 0;
	treal sum6 = 0;
	treal sum7 = 0;

	for(int n = 1; n<= end_pt-start_pt+1; n++) 
	{
		sum1 += x_sec_sq(n) * y_sec(n);
		sum2 += y_sec(n) * y_sec_ln(n);
		sum3 += x_sec(n) * y_sec(n);
		sum4 += x_sec(n) * y_sec(n) * y_sec_ln(n);
		sum5 += y_sec(n);
		sum6 += x_sec_sq(n) * y_sec(n);
	}
	sum7 = sum3*sum3;

	a = ( sum1 * sum2 - sum3 * sum4 ) / ( sum5 * sum6 - sum7 );

	//log.LogMessage(LOG_INFO, "Water amplitude is: %.2f", exp(a));
	return exp(a);
}
