#include <iostream>
#include "CBasis.hpp"
#include "Workspace.hpp"
#include "tarquin.hpp"
#include "cvm_util.hpp"
#include "preprocess.hpp"
#include "console.hpp"
#include "fitviewer.hpp"
#include "CBoswell.hpp"
#include "td_conv_ws.hpp"
#include <algorithm>
#include "proto_buf.hpp"
#include "export_data.hpp"
#include "exception.hpp"
#include "fidio/CFIDReader.hpp"

using namespace tarquin;

int main(int argc, char* argv[])
{
    // this contains all final and intermediate results
    Workspace workspace;

    // this is where all messages go
    CBoswell log(LOG_STDOUT);

    // this is the FID we will be preprocessing, then fitting to
    CFID& fidraw = workspace.GetFIDRaw();

    // options from the command line
    Options& options = workspace.GetOptions();

	if( false == ParseCommandLine(argc, argv, options, fidraw) ) 
		return -1;

	try
	{
		//std::cout << fidraw.GetPPMRef() << std::endl << std::flush;

		// are we passing control to the fitviewer code?
		if( options.GetViewFile() != "" ) 
		{
			//std::cout << options.GetPlotSigs() << std::endl;
			return FitView(options.GetViewFile(), options.GetPlotSigs(), 
					options.GetPPMstart(), options.GetPPMend(), options.Getlb(), options.GetPause());
		}

        // save some paras before they get corrupted
        CFID orig_fid = fidraw;

        // try and guess the format if not set
        if ( options.GetFormat() == tarquin::NOTSET )
        {
            std::string fn = options.GetFilename();
            // try to load as DICOM
            try
            {
                options.SetFormat(tarquin::DCM);
                fidraw.Load(fn, options, workspace, log);
            }
            catch( const tarquin::Exception& /*e*/ )
            {
                // remove any bad data fields from the previous load attempt
                fidraw = orig_fid;

                // doesn't look like DICOM so guess from filename
                if ( (fn.substr(fn.size()-5, fn.size()) == ".sdat") 
                        || (fn.substr(fn.size()-5, fn.size()) == ".SDAT")
                        || (fn.substr(fn.size()-5, fn.size()) == ".spar")
                        || (fn.substr(fn.size()-5, fn.size()) == ".SPAR" ) )
                    options.SetFormat(tarquin::PHILIPS);
                else if ( fn.substr(fn.size()-4, fn.size()) == ".shf" )
                    options.SetFormat(tarquin::SHF);
                else if ( fn.substr(fn.size()-2, fn.size()) == ".7" )
                    options.SetFormat(tarquin::GE);
                else if ( fn.substr(fn.size()-4, fn.size()) == ".rda"
                        || fn.substr(fn.size()-4, fn.size()) == ".RDA" )
                    options.SetFormat(tarquin::RDA);
                else
                    options.SetFormat(tarquin::SIEMENS);
                
                // try again
                fidraw.Load(fn, options, workspace, log);
            }
        }
        else
        {
		    fidraw.Load(options.GetFilename(), options, workspace, log);
        }


        // are we passing control to the export paras code?
        if( options.GetPrintParas() == true ) 
        {
            std::cout.precision(0);
            std::cout << std::endl;
            std::cout << "Data format                 : ";
            if( tarquin::SIEMENS == options.GetFormat() )
                std::cout << "Siemens IMA";
            else if( tarquin::PHILIPS == options.GetFormat() )
                std::cout << "Philips spar/sdat";
            else if( tarquin::SHF == options.GetFormat() )
                std::cout << "GE SHF";
            else if( tarquin::GE == options.GetFormat() )
                std::cout << "GE p-file";
            else if( tarquin::RDA == options.GetFormat() )
                std::cout << "Siemens RDA";
            else if( tarquin::DCM == options.GetFormat() )
                std::cout << "DICOM MRS";
            else if( tarquin::DANGER == options.GetFormat() )
                std::cout << "Dangerplot";
            else if( tarquin::VARIAN == options.GetFormat() )
                std::cout << "Varian";

            std::cout << std::endl;
            std::cout << "Transmitter frequency (Hz)  : " << std::fixed << fidraw.GetTransmitterFrequency() << std::endl;
            std::cout.precision(0);
            std::cout << "Sampling frequency (Hz)     : " << fidraw.GetSamplingFrequency() << std::endl;
            std::cout.precision(3);
            std::cout << "Echo time (s)               : " << fidraw.GetEchoTime() << std::endl;
            std::cout.precision(3);
            std::cout << "Number of points            : " << fidraw.GetNumberOfPoints() << std::endl;
            std::cout << "Rows                        : " << fidraw.GetRows() << std::endl;
            std::cout << "Columns                     : " << fidraw.GetCols() << std::endl;
            std::cout << "Slices                      : " << fidraw.GetSlices() << std::endl;
            std::cout << "Voxel dimensions            : ";

            if ( fidraw.IsKnownVoxelDim() )
            {
                std::vector<double> voxel_dims = fidraw.GetVoxelDim();
                std::cout << voxel_dims[0] << " " << voxel_dims[1] << " " << voxel_dims[2] << std::endl;
            }
            else
                std::cout << "Unknown" << std::endl;

            std::cout << "VOI dimensions              : ";

            if ( fidraw.IsKnownVoiDim() )
            {
                std::vector<double> voi_dims = fidraw.GetVoiDim();
                std::cout << voi_dims[0] << " " << voi_dims[1] << " " << voi_dims[2] << std::endl;
            }
            else
                std::cout << "Unknown" << std::endl;

            return 0;
        }
        
        // looks like the main algorithm is being run so display the splash
        DisplaySplash();

		// if ref has been defined from the c/l override the value from the fid
		if( options.GetRefSpec() ) 
		{
			fidraw.SetPPMRef(options.GetRef());
		}
		// otherwise use the default of 4.7 if the file format is not dpt
		else 
		{
			if( options.GetFormat() != tarquin::DANGER ) 
			{
				fidraw.SetPPMRef(options.GetRef());
			}
		}
		//std::cout << fidraw.GetPPMRef() << std::endl;
		// this logic will cause the ref to be taken from the dpt file unless overridden by the c/l option

		// are we writing the original file out as dangerplot format?
		if( options.GetOutRawFile() != "" ) 
		{

			log.Out(LOG_INFO) << "Writing raw FID to file \"" << options.GetOutRawFile() << "\".";

			if( !fidraw.SaveToFile(options.GetOutRawFile()) ) 
			{
				log.Out(LOG_ERROR) << "failed to write file.";
				return -1;
			}

			log.Out(LOG_INFO) << "Wrote file: \"" << options.GetOutRawFile() << "\"";
		}
        
        // are we writing the original file out as lcm format?
		if( options.GetOutRawFileLcm() != "" ) 
		{

			log.Out(LOG_INFO) << "Writing raw FID to file \"" << options.GetOutRawFileLcm() << "\".";

			if( !fidraw.SaveToFileLCM(options.GetOutRawFileLcm()) ) 
			{
				log.Out(LOG_ERROR) << "failed to write file.";
				return -1;
			}

			log.Out(LOG_INFO) << "Wrote file: \"" << options.GetOutRawFileLcm() << "\"";
		}


		// are we loading a reference FID?
		if( options.GetFilenameWater() != "" || fidraw.GetCWF() == true ) 
		{
			CFID& fidWater = workspace.GetFIDWater();

			// initialise the water FID to have the same parameters as the WS FID
			fidWater.SetParametersFromFID(fidraw);

            if ( fidraw.GetCWF() == true )
            {
                options.SetFilenameWater(options.GetFilename());
                fidWater.LoadW(options.GetFilenameWater(), options, log);
            }
            else
            {
                fidWater.Load(options.GetFilenameWater(), options, workspace, log);
            }
			
            if( options.GetOutRawFileW() != "" ) 
            {

                log.Out(LOG_INFO) << "Writing raw FID to file \"" << options.GetOutRawFileW() << "\".";

                if( !fidWater.SaveToFile(options.GetOutRawFileW()) ) 
                {
                    log.Out(LOG_ERROR) << "failed to write file.";
                    return -1;
                }

                log.Out(LOG_INFO) << "Wrote file: \"" << options.GetOutRawFileW() << "\"";
            }
            
            if( options.GetOutRawFileLcmW() != "" ) 
            {

                log.Out(LOG_INFO) << "Writing raw FID to file \"" << options.GetOutRawFileLcmW() << "\".";

                if( !fidWater.SaveToFileLCM(options.GetOutRawFileLcmW()) ) 
                {
                    log.Out(LOG_ERROR) << "failed to write file.";
                    return -1;
                }

                log.Out(LOG_INFO) << "Wrote file: \"" << options.GetOutRawFileLcmW() << "\"";
            }

		}


		//plot(fidWater)

		// output summary message of parameters (before preprocessing)
		log.Out(LOG_INFO) << "\nParameters Before Preprocessing";
		log.Out(LOG_INFO) << "\n-------------------------------";
		log.Out(LOG_INFO) << workspace;
		log.Flush();

		// check whether sufficient parameters are known
		if( false == workspace.CheckParameters() ) 
		{
			log.Out(LOG_ERROR) << "Not all necessary parameters are known.";
			log.Out(LOG_ERROR) << "Please check above messages to see which ones are missing.";
			return -1;
		}

		// this is the FID after some preprocessing, the FID we will actually fit to
		CFID& fidproc = workspace.GetFIDProc();

		Preprocessor preprocess(workspace, log);

		// preprocess the raw fid
		preprocess(); 
    
		// are we writing the original file out as dangerplot format?
		if( options.GetOutPreFile() != "" ) 
		{
			log.Out(LOG_INFO) << "\nWriting preprocessed FID to file \"" << options.GetOutPreFile() << "\".";

			if( false == fidproc.SaveToFile(options.GetOutPreFile()) ) 
			{
				log.Out(LOG_ERROR) << "failed to write file.";
				return -1;
			}

			log.Out(LOG_INFO) << "Wrote file: \"" << options.GetOutPreFile() << "\"";
		}

		// output summary message of parameters (after preprocessing)
		log.Out(LOG_INFO) << "\n\nParameters After Preprocessing";
		log.Out(LOG_INFO) << "\n------------------------------";
		log.Out(LOG_INFO) << workspace;
		log.Flush();

		// are we going to display the preprocessed FID?
		if( true == options.GetShowPreprocessed() )
			PlotFIDs(fidraw, fidproc);

		// the basis we will adapt to fit to the signal
		CBasis& basis = workspace.GetBasis();

		// load or simulate a basis to match this FID
		if( options.GetUsePrecompiled() ) 
		{
            int len = options.GetBasisPath().size();
            std::string str_end;
            if ( len > 6 )
                str_end = options.GetBasisPath().substr(len-6,6);
            else
                str_end = "";


            if (( str_end == ".basis" ) || ( str_end == ".BASIS" ))
            {
                if( false == basis.ReadLCMBasis(options.GetBasisPath(), fidraw, options, log) )
                {
                    log.Out(LOG_ERROR) << "failed to load basis file: \"" << options.GetBasisPath() << "\"";
                    return -1;
                }
            }
            else
            {
                if( false == load_basis(options.GetBasisPath(), basis) ) 
                {
                    log.Out(LOG_ERROR) << "failed to load basis file: \"" << options.GetBasisPath() << "\"";
                    return -1;
                }
            }

		}
		else 
		{

			// simulate an basis from csv
            if ( options.GetBasisPath() != "" )
            {
                if( false == basis.Simulate(options.GetBasisPath(), fidproc, options, log) ) 
                {
                    log.Out(LOG_ERROR) << "simulation failed, no basis available so quitting";
                    return -1;
                }
            }
            else
            {
                // internal basis set
                if( false == basis.Simulate(fidproc, options, log) ) 
                {
                    log.Out(LOG_ERROR) << "internal simulation failed, no basis available so quitting";
                    return -1;
                }
            }

			//basis.ApplyRef(fidproc.GetPPMRef());	

			// did the user specify to save it?
			if( options.GetBasisSaveFile() != "" ) 
			{
				basis.Initialise(4.65f);

				log.Out(LOG_INFO) << "\nSaving basis to \"" << options.GetBasisSaveFile() << "\".";

				if( false == save_basis(options.GetBasisSaveFile().c_str(), basis) ) 
				{
					log.Out(LOG_ERROR) << "failed to write basis file, quitting";
					return -1;
				}

				log.Out(LOG_INFO) << "\nWrote \"" << options.GetBasisSaveFile() << "\".";
			}
    
		}
        
        if( options.GetBasisSaveFileCSV() != "" ) 
        {
            basis.SaveTxt(options.GetBasisSaveFileCSV().c_str(),fidraw);
        }

        if( options.GetBasisSaveFileLCM() != "" ) 
        {

            if( options.GetUsePrecompiled() ) 
            {
                basis.SaveLCM(options.GetBasisSaveFileLCM().c_str(),fidproc);
            }
            else
            {
                CFID lcmfid = fidproc;
                CBasis lcm_basis;
                lcmfid.SetNumberOfPoints(options.GetLCMBasisPts());

                // resimulate with different N before writing file
                if ( options.GetBasisPath() != "" )
                {
                    if( false == lcm_basis.Simulate(options.GetBasisPath(), lcmfid, options, log) ) 
                    {
                        log.Out(LOG_ERROR) << "simulation failed, no basis available so quitting";
                        return -1;
                    }
                }
                else
                {
                    // internal basis set
                    if( false == lcm_basis.Simulate(lcmfid, options, log) ) 
                    {
                        log.Out(LOG_ERROR) << "internal simulation failed, no basis available so quitting";
                        return -1;
                    }
                }
                lcm_basis.SaveLCM(options.GetBasisSaveFileLCM().c_str(),lcmfid);
            }
        }

        // check that basis and fidraw match
        if ( false == basis.check(fidraw, log) )
        {
		    log.Out(LOG_ERROR) << "Basis set and fid parameter mismatch, please resimulate\n";
			return -1;
        }

        if ( options.GetNoFit() == true )
        {
            return 0;
        }

		// do the analysis
		if( false == RunTARQUIN(workspace, log) ) 
		{
			log.Out(LOG_WARNING) << "TARQUIN could not fit this data (low quality?).";
			return -1;
		}

		// serialise results 
		std::string OutputXMLPath = options.GetOutputXMLPath();

		// since this now has a default, it must always be set to something
		//assert( OutputXMLPath.size() > 0 );

		// try and save the results
        if ( OutputXMLPath.size() != 0 )
        {
            if( false == save_results( OutputXMLPath.c_str(), workspace ) ) 
            {
                log.Out(LOG_ERROR) << "failed to write results file to:\"" << OutputXMLPath << "\"";
                return -1;
            }
        }

		// are we outputting a dpt file of the signal with TARQUIN phasing?
		if( options.GetOutPostFile() != "" ) 
		{
			log.Out(LOG_INFO) << "\nWriting postprocessed FID to file \"" << 
				options.GetOutPostFile() << "\".";

			CFID yfid = workspace.GetFID();
			if( false == yfid.SaveToFile(options.GetOutPostFile()) ) 
			{
				log.Out(LOG_ERROR) << "failed to write file";
				return -1;
			}

			log.Out(LOG_INFO) << "\nWrote \"" << options.GetOutPostFile() << "\".";
		}


		// TODO this function should be depreciated
		if( options.GetFilenameImag() != "" ) 
		{
			// TODO: Shift all this to a function somewhere

			CFID yfid = workspace.GetFID();
			cvm::cvector y = yfid.GetVectorFID(0); //TODO
			cvec_stdvec yhat_vec = workspace.GetSignalEstimate();
			cvm::cvector yhat = yhat_vec[0];
            
            int zf = 1;
            // zero fill if less than 4096 points
            if ( y.size() < 4096 )
                zf = int(round(4096/y.size()));

            y.resize(y.size()*zf);
            yhat.resize(yhat.size()*zf);

            // copy last pts points of real fid to end of zfilled fid
            int pts = 5;
            if ( zf > 1 )
            {
                for ( int n = 1 ; n < pts + 1; n++)
                {
                    y(y.size()-pts+n) = y(y.size()/zf-pts+n);
                    y(y.size()/zf-pts+n) = 0;
                }
            }

            cvm::cvector Y(y.size());
            cvm::cvector YHAT(yhat.size());

			fft(y, Y);
			fft(yhat, YHAT);

			//cvm::cvector Y = fft(y);
			//cvm::cvector YHAT = fft(yhat);

			Y = fftshift(Y);
			YHAT = fftshift(YHAT);

			cvm::cvector residual = y-yhat;
			cvm::cvector baseline = residual;

			/* old method
			// set baseline points not considered in fitting to be zero	
			for ( int n = options.GetRangeStart(); n < baseline.size()+1-11; n ++ )
			baseline(n) = 0;

			cvm::cvector BASELINE = fft(baseline);
			BASELINE = fftshift(BASELINE);
			*/
			//cvm::cvector BASELINE = fft(baseline);
			//BASELINE = fftshift(BASELINE);
			//BASELINE,

			cvm::cvector RESIDUAL = fft(residual);
			RESIDUAL = fftshift(RESIDUAL);

			cvm::cvector BASELINE;
			td_conv_ws( RESIDUAL, BASELINE, 100, 10);	

			cmat_stdvec s_vec = workspace.GetBasisMatrix();
			cvm::cmatrix s = s_vec[0];
            s.resize(s.msize()*zf, s.nsize());
            
			//cvm::cmatrix& s = workspace.GetGroupMatrix();
			rvec_stdvec ahat_vec = workspace.GetAmplitudes();
			cvm::rvector ahat = ahat_vec[0];

			for ( int n = 1; (n < ahat.size() + 1); n++) 
				for ( int m = 1; (m < s.msize() + 1); m++)
					s(m,n) = s(m,n)*ahat(n);
				

			cvm::cmatrix S = fft(s);
			S = fftshift(S);
			coord voxel(1, 1, 1); // TODO
			cvm::rvector freq_scale = yfid.GetPPMScale(voxel,zf);

			cvm::cmatrix all(S.msize(),S.nsize()+7);
			all.assign(1,8,S);

			// find points corresponding to 0.2 and 4.0 ppm
			int left = 1, right = 1;
			for ( int n = 1; n < (freq_scale.size()); n++ ) 
			{
				if ( ( freq_scale(n) > 4 ) && ( freq_scale(n+1) < 4 ) )
					left = n;
				if ( ( freq_scale(n) > 0.2 ) && ( freq_scale(n+1) < 0.2 ) )
					right = n;
			}
			// std::cout << left << std::endl;
			// std::cout << right << std::endl;

			double Ymax = 0;
			// find max of y for plotting residual	
			for ( int n = left; n < right+1; n++ ) 
			{
				if ( Ymax < Y(n).real() )
					Ymax = Y(n).real();
				if ( Ymax < YHAT(n).real()+BASELINE(n).real() )
					Ymax = YHAT(n).real()+BASELINE(n).real();
			}

			double res_min = 0;
			// find min of residual for plotting residual	
			for ( int n = left; n < (right+1); n++ ) 
				if ( res_min > RESIDUAL(n).real() - BASELINE(n).real() )
					res_min = RESIDUAL(n).real() - BASELINE(n).real() ;

			double res_max = 0;
			// find max of residual for plotting residual	
			for ( int n = left; n < (right+1); n++ ) 
				if ( res_max < RESIDUAL(n).real() - BASELINE(n).real() )
					res_max = RESIDUAL(n).real() - BASELINE(n).real() ;

			for ( int m = 1; m < (S.msize()+1); m++ ) 
			{
				all(m,1) = Y(m);
				all(m,2) = YHAT(m)+BASELINE(m);
				all(m,3) = RESIDUAL(m) - BASELINE(m) + Ymax - res_min;
				all(m,4) = Ymax-res_min;
				all(m,5) = Ymax;
				all(m,6) = Ymax-res_min+res_max;
				all(m,7) = BASELINE(m);
				for ( int n = 8; n < all.nsize() + 1 ; n++ )
					all(m,n) = all(m,n) + BASELINE(m);
			}

			// and plot results	
			savefitimag(freq_scale, all, options.GetFilenameImag());

		}
    
    	if( options.GetFilenamePdf() != "" ) 
            ExportPdfResults(options.GetFilenamePdf(), workspace);

		// are we outputting a txt file of results?
		if( options.GetFilenameTxt() != "" )
			ExportTxtResults(options.GetFilenameTxt(), workspace);

		if( options.GetFilenameCSV() != "" ) 
			ExportCsvResults(options.GetFilenameCSV(), workspace);
    
        if( options.GetFilenameCSVFit() != "" ) 
			ExportCsvFit(options.GetFilenameCSVFit(), workspace);

        if( options.GetFilenameCSVSpectraAligned() != "" ) 
			ExportCsvSpectraAligned(options.GetFilenameCSVSpectraAligned(), workspace);

        if( options.GetFilenameCSVSpectraAlignedMag() != "" ) 
			ExportCsvSpectraAligned(options.GetFilenameCSVSpectraAlignedMag(), workspace, true);

		// output some convergence information
        // TODO fix [0]s below for CSI data
		std::vector<std::vector < double > > info = workspace.GetLMinfo();
		std::cout << "\n";
		std::cout << "\nOptimisation details";
		std::cout << "\n--------------------";
		std::cout << "\nl2 norm of error at initial p  = " << info[0][0];
		std::cout << "\nl2 norm of error at final p    = " << info[0][1];
		std::cout << "\nl2 norm of J.'*e at final p    = " << info[0][2];
		std::cout << "\nl2 norm of D*p at final p      = " << info[0][3];
		std::cout << "\nnumber of iterations           = " << info[0][5];
		std::cout << "\nnumber of function evaluations = " << info[0][6];
		std::cout << "\nnumber of Jacobian evaluations = " << info[0][7];

		if( 1 == info[0][6] )
			std::cout << "\nstopped by small gradient\n";
		else if( 2 == info[0][6] )
			std::cout << "\nstopped by small D*p\n";
		else if( 3 == info[0][6] )
			std::cout << "\nstopped by iteration limit\n";
		else if( 4 == info[0][6] )
			std::cout << "\nstopped by singular matrix - restart with bigger mu";
		else if( 5 == info[0][6] )
			std::cout << "\nno further reduction possible - restart with bigger mu\n";
		else if( 6 == info[0][6] )
			std::cout << "\nstopping by small l2 norm of error\n";

		std::cout << "\n\nTARQUIN has finished. Please use \"--view_fit\" to inspect the results.";
		std::cout << "\n" << std::endl;
	}
	catch( const tarquin::CFIDReader::Exception& e )
	{
		std::cerr << "\nerror loading FID: " << e.what() << std::endl;
		return -1;
	}
	catch( const tarquin::Exception& e )
	{
		std::cerr << "\ngeneral error: " << e.what() << std::endl;
		return -1;
	}
	catch( ... )
	{
		std::cerr << "\nunknown exception (FAILURE!)" << std::endl;
		return -1;
	}

	return 0;
}



