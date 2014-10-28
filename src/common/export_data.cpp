#include <iostream>
#include "CBasis.hpp"
#include "Workspace.hpp"
#include "tarquin.hpp"
#include "cvm_util.hpp"
#include <algorithm>
#include <string>
#include "td_conv_ws.hpp"
#include <sstream>
#include <iomanip>
#include "version/version.h"

using namespace tarquin; 

void ExportCsvSpectraAligned(const std::string& strFilename, const Workspace& workspace, bool mag)
{
	// create a new options 
	Options options = workspace.GetOptions();

	CFID yfid = workspace.GetFID();
	coord_vec fit_list = options.GetFitList();

    std::ofstream fout(strFilename.c_str(), std::ios::out);

    int spectra = fit_list.size();
    cvm::rvector freq_scale = yfid.GetPPMScale(fit_list[0], options.GetZF());
    int N = freq_scale.size();
    
    std::vector<cvm::cvector> spectra_vec;

    for ( size_t p = 0; p < spectra; p++ )
    {
        yfid.ShiftRef(options.GetRef(), yfid.vox2ind(fit_list[p]));
        cvm::cvector y = yfid.GetVectorFID(fit_list[p]);
        lb(y,workspace.GetFIDRaw(),options.Getlb());
        ZeroPad(y,options.GetZF());
        cvm::cvector Y = fft(y);
        Y = fftshift(Y);
        spectra_vec.push_back(Y);
    }

    fout << "PPMScale";
    for ( size_t n = 0; n < spectra; n++ )
        fout << "," << fit_list[n].row  << "-" << fit_list[n].col << "-" << fit_list[n].slice;

    fout << std::endl;

    for ( size_t row = 1; row < N+1; row++ )
    {
        fout << freq_scale(row);
        for ( size_t col = 1; col < spectra+1; col++ )
        {
            if ( mag == false )
                fout << "," << spectra_vec[col-1](row).real();
            else
                fout << "," << std::abs(spectra_vec[col-1](row));

        }
        fout << std::endl;
    }
}

/* OLDER FUNCTION, bit like the above but transposed so difficult to read in other software
void ExportCsvSpectraAligned(const std::string& strFilename, const Workspace& workspace, bool mag)
{
	// create a new options 
	Options options = workspace.GetOptions();

	CFID yfid = workspace.GetFID();
	coord_vec fit_list = options.GetFitList();

    std::ofstream fout(strFilename.c_str(), std::ios::out);
    
    for ( size_t p = 0; p < fit_list.size(); p++ )
    {
        // align data
        yfid.ShiftRef(options.GetRef(), yfid.vox2ind(fit_list[p]));

        cvm::cvector y = yfid.GetVectorFID(fit_list[p]);

        lb(y,workspace.GetFIDRaw(),options.Getlb());
        ZeroPad(y,options.GetZF());

        cvm::cvector Y = fft(y);
        Y = fftshift(Y);
        cvm::rvector freq_scale = yfid.GetPPMScale(fit_list[p], options.GetZF());

        // write ppm scale first
        if ( p == 0 )
        {
            fout << ",,,,,PPMScale";
            for ( int m = 1; ( m < Y.size() + 1 ); m++) 
            {
                fout << "," << freq_scale(m);
            }
            fout << std::endl;
        }

        fout << "Row : " << fit_list[p].row << ", Col : " << fit_list[p].col << ", Slice : " << fit_list[p].slice;

        // write the data points
        for ( int m = 1; ( m < Y.size() + 1 ); m++) 
        {
            if ( mag == false )
                fout << "," << Y(m).real();
            else
                fout << "," << std::abs(Y(m));
        }
        fout << std::endl;
    }
}*/

void ExportCsvSpectrum(const std::string& strFilename, const Workspace& workspace)
{
	// create a new options 
	Options options = workspace.GetOptions();

	CFID yfid = workspace.GetFID();
	coord_vec fit_list = options.GetFitList();

    std::ofstream fout(strFilename.c_str(), std::ios::out);
    
    for ( size_t p = 0; p < fit_list.size(); p++ )
    {

        fout << "Row : " << fit_list[p].row << ", Col : " << fit_list[p].col << ", Slice : " << fit_list[p].slice << std::endl;

        cvm::cvector y = yfid.GetVectorFID(fit_list[p]);

        lb(y,workspace.GetFIDRaw(),options.Getlb());
        ZeroPad(y,options.GetZF());

        cvm::cvector Y = fft(y);
        Y = fftshift(Y);
        cvm::rvector freq_scale = yfid.GetPPMScale(fit_list[p], options.GetZF());

        // export freq_scale, Y, YHAT, BASELINE and S as a csv file
        // write the first line with column headings
        fout << "PPMScale,Data";

        fout << std::endl;

        // write the data points
        for ( int m = 1; ( m < Y.size() + 1 ); m++) 
        {
            fout << freq_scale(m) << ",";
            fout << Y(m).real();
            //fout << Y(m).imag();
            //fout << abs(Y(m));
            fout << std::endl;
        }
    }
}

void ExportCsvFit(const std::string& strFilename, const Workspace& workspace)
{
	// create a new options 
	Options options = workspace.GetOptions();

	const CFID& yfid = workspace.GetFID();
	coord_vec fit_list = options.GetFitList();
	const cvec_stdvec& yhat_vec = workspace.GetSignalEstimate();

    std::ofstream fout(strFilename.c_str(), std::ios::out);
    
    for ( size_t p = 0; p < yhat_vec.size(); p++ )
    {

        fout << "Row : " << fit_list[p].row << ", Col : " << fit_list[p].col << ", Slice : " << fit_list[p].slice << std::endl;

        cvm::cvector y = yfid.GetVectorFID(fit_list[p]);
        cvm::cvector yhat = yhat_vec[p];

        lb(y,workspace.GetFIDRaw(),options.Getlb());
        ZeroPad(y,options.GetZF());

        lb(yhat,workspace.GetFIDRaw(),options.Getlb());
        ZeroPad(yhat,options.GetZF());

        cvm::cvector Y = fft(y);
        cvm::cvector YHAT = fft(yhat);

        Y = fftshift(Y);
        YHAT = fftshift(YHAT);

        cvm::cvector residual = y-yhat;
        cvm::cvector baseline = residual;

        cvm::cvector RESIDUAL = fft(residual);
        RESIDUAL = fftshift(RESIDUAL);

        cvm::cvector BASELINE;
        td_conv_ws( RESIDUAL, BASELINE, options.GetBL()*options.GetZF(), 10);	

        const cmat_stdvec& s_vec = workspace.GetBasisMatrix(); // TODO should be a reference?
        cvm::cmatrix s = s_vec[p];

        CBasis basis = workspace.GetBasis();
        std::vector<std::string> sig_names = basis.GetSignalNames();
        
        rvec_stdvec ahat_vec = workspace.GetAmplitudes();
        cvm::rvector ahat = ahat_vec[p];

        cvm::cmatrix s_proc(s.msize()*options.GetZF(), s.nsize());
        for ( int n = 1; (n < ahat.size() + 1); n++) 
        {
            cvm::cvector s_temp = s(n);
            lb(s_temp,workspace.GetFIDRaw(),options.Getlb());
            ZeroPad(s_temp,options.GetZF());
            s_proc(n) = s_temp*ahat(n);
        }

        cvm::cmatrix S = fft(s_proc);

        S = fftshift(S);
        cvm::rvector freq_scale = yfid.GetPPMScale(fit_list[p], options.GetZF());

        // export freq_scale, Y, YHAT, BASELINE and S as a csv file
        // write the first line with column headings
        fout << "PPMScale,Data,Fit,Baseline";

        if ( options.GetExtCSVFit() )
        {
            fout << ",";
            for ( size_t n = 0; n < sig_names.size() - 1; n++ )
                fout << sig_names[n] << ",";
            
            // last one without a ","
            fout << sig_names[sig_names.size()-1];
        }

        fout << std::endl;

        // write the data points
        for ( int m = 1; ( m < Y.size() + 1 ); m++) 
        {
            fout << freq_scale(m) << ",";
            fout << Y(m).real() << ",";
            fout << YHAT(m).real() << ",";
            fout << BASELINE(m).real();
            
            if ( options.GetExtCSVFit() )
            {
                fout << ",";
                for ( int n = 1; ( n < s.nsize() ); n++) 
                    fout << S(m,n).real() << ",";
                // last one without a ","
                fout << S(m, s.nsize()).real();
            }

            fout << std::endl;
        }
    }
}

void ExportCsvResults(const std::string& strFilename, const Workspace& workspace)
{
    const CBasis& basis     = workspace.GetBasis();
	const rvec_stdvec ahat  = workspace.GetAmplitudesNormalised();
	const rvec_stdvec crlb  = workspace.GetCRLBsNormalised();
	const Options& options  = workspace.GetOptions();
    const CFID& yfid        = workspace.GetFID();

	std::ofstream fout(strFilename.c_str(), std::ios::out);
	fout << std::showpoint;
	fout << std::scientific;

	std::vector < std::string > signal_names = basis.GetSignalNames();	
	std::vector < std::string > signal_names_sorted = signal_names;

	std::sort(signal_names_sorted.begin(), signal_names_sorted.end()); 
	int width = ahat[0].size();

	fout << "Signal amplitudes" << std::endl;
	// output the vol location column headings
	fout << "Row" << "," << "Col" << "," << "Slice" << ","; 
	// output the signal names	
	for(integer n = 1; n < width; n++) 
		fout << signal_names_sorted[n-1] << ",";
	
    if ( workspace.GetAmplitudesNormalisedComb().size() > 0 )
    {
        fout << signal_names_sorted[width-1] << ",";
        for(integer n = 1; n < workspace.GetMetabNamesComb().size(); n++)
            fout << workspace.GetMetabNamesComb()[n-1] << ",";
        
        fout << workspace.GetMetabNamesComb()[workspace.GetMetabNamesComb().size()-1] << std::endl;
    }
    else
    {
        // last one without the comma	
        fout << signal_names_sorted[width-1] << std::endl;
    }
	
	// output the signal amplitudes
    const coord_vec& fit_list = options.GetFitList();
	coord_vec::const_iterator i_coord = fit_list.begin();

    size_t fit = 0;
	for( rvec_stdvec::const_iterator i = ahat.begin(); i != ahat.end(); ++i )	
	{
		// write the voxel info
		fout << (*i_coord).row << "," << (*i_coord).col << "," << (*i_coord).slice << ","; 
		double tempamp = 0;
		for(integer n = 1; n < width; n++) {
			// find index of ahat now the list is sorted
			for(integer m = 1; m < width+1; m++) 
			{
				if ( signal_names_sorted[n-1] == signal_names[m-1] ){
					tempamp = (*i)(m);
					break;
				}
			}
			fout << tempamp << ",";	
		}
        
        if ( workspace.GetAmplitudesNormalisedComb().size() > 0 )
        {
            for(integer m = 1; m < width+1; m++) 
            {
                if ( signal_names_sorted[width-1] == signal_names[m-1] )
                {
                    tempamp = (*i)(m);
                    break;
                }
            }
            fout << tempamp << ",";

            for(integer n = 1; n < workspace.GetAmplitudesNormalisedComb()[fit].size(); n++)
            {
                fout << workspace.GetAmplitudesNormalisedComb()[fit][n] << ",";
            }
            fout << workspace.GetAmplitudesNormalisedComb()[fit][workspace.GetAmplitudesNormalisedComb()[fit].size()] << std::endl;
        }
        else
        {
            for(integer m = 1; m < width+1; m++) 
            {
                if ( signal_names_sorted[width-1] == signal_names[m-1] )
                {
                    tempamp = (*i)(m);
                    break;
                }
            }
            fout << tempamp << std::endl;
        }

		// move to the next voxel
		i_coord++;
        fit++;
	}

	fout << "CRLBs (standard deviation)" << std::endl;
	fout << "Row" << "," << "Col" << "," << "Slice" << ","; 
	// output the signal names	
	for(integer n = 1; n < width; n++) 
	{
		fout << signal_names_sorted[n-1] << ",";
	}
    
    if ( workspace.GetAmplitudesNormalisedComb().size() > 0 )
    {
        fout << signal_names_sorted[width-1] << ",";
        for(integer n = 1; n < workspace.GetMetabNamesComb().size(); n++)
            fout << workspace.GetMetabNamesComb()[n-1] << ",";
        
        fout << workspace.GetMetabNamesComb()[workspace.GetMetabNamesComb().size()-1] << std::endl;
    }
    else
    {
        // last one without the comma	
        fout << signal_names_sorted[width-1] << std::endl;
    }

	// output the signal amplitudes
	i_coord = fit_list.begin();
    fit = 0;
	for( rvec_stdvec::const_iterator i = crlb.begin(); i != crlb.end(); ++i )	
	{
		// write the voxel info
		fout << (*i_coord).row << "," << (*i_coord).col << "," << (*i_coord).slice << ","; 
		double tempamp = 0;
		for(integer n = 1; n < width; n++) {
			// find index of ahat now the list is sorted
			for(integer m = 1; m < width+1; m++) 
			{
				if ( signal_names_sorted[n-1] == signal_names[m-1] ){
					tempamp = (*i)(m);
					break;
				}
			}
			fout << tempamp << ",";	
		}

        if ( workspace.GetAmplitudesNormalisedComb().size() > 0 )
        {
            for(integer m = 1; m < width+1; m++) 
            {
                if ( signal_names_sorted[width-1] == signal_names[m-1] )
                {
                    tempamp = (*i)(m);
                    break;
                }
            }
            fout << tempamp << ",";

            for(integer n = 1; n < workspace.GetCRLBsNormalisedComb()[fit].size(); n++)
            {
                fout << workspace.GetCRLBsNormalisedComb()[fit][n] << ",";
            }
            fout << workspace.GetCRLBsNormalisedComb()[fit][workspace.GetCRLBsNormalisedComb()[fit].size()] << std::endl;
        }
        else
        {
            for(integer m = 1; m < width+1; m++) 
            {
                if ( signal_names_sorted[width-1] == signal_names[m-1] )
                {
                    tempamp = (*i)(m);
                    break;
                }
            }
            fout << tempamp << std::endl;
        }


		// move to the next voxel
		i_coord++;
        fit++;
	}
    
    //fout << std::showpoint;
    //fout << std::fixed;
    //fout << std::setprecision(6);

    fout << "Fit diagnostics" << std::endl;
	fout << "Row" << "," << "Col" << "," << "Slice" << ","; 
	// output the row column names
	fout << "Q,max res,metab ratio,peak metab ratio,metab FWHM (PPM),metab FWHM (Hz),SNR,SNR max,SNR metab,spec noise,td noise,ref,init beta,final beta,final beta (PPM),phi0 (deg),phi1 (deg/PPM),water amp,water FWHM (Hz),water FWHM (PPM),water freq (Hz),baseline dev,max baseline,min baseline,baseline shape,initial residual,final residual,stopping reason" << std::endl;
	
	// output fitting info
    // Q
    const std::vector<double>& Q_vec = workspace.GetQ();
    std::vector<double>::const_iterator i_Q = Q_vec.begin();

    const std::vector<double>& Q_vec_rel = workspace.GetQ_rel();
    std::vector<double>::const_iterator i_Q_rel = Q_vec_rel.begin();

    const std::vector<double>& metab_rat = workspace.GetMetabRat();
    std::vector<double>::const_iterator i_metab_rat = metab_rat.begin();

    const std::vector<double>& peak_metab_rat = workspace.GetPeakMetabRat();
    std::vector<double>::const_iterator i_peak_metab_rat = peak_metab_rat.begin();

    const std::vector<double>& metab_fwhm_vec = workspace.GetMetabFWHM();
    std::vector<double>::const_iterator i_metab_fwhm = metab_fwhm_vec.begin();

    const std::vector<double>& BLV_vec = workspace.GetBLV();
    std::vector<double>::const_iterator i_BLV = BLV_vec.begin();

    const std::vector<double>& max_bl_vec = workspace.GetMaxBL();
    std::vector<double>::const_iterator i_max_bl = max_bl_vec.begin();

    const std::vector<double>& min_bl_vec = workspace.GetMinBL();
    std::vector<double>::const_iterator i_min_bl = min_bl_vec.begin();

    const std::vector<double>& BLS_vec = workspace.GetBLS();
    std::vector<double>::const_iterator i_BLS = BLS_vec.begin();

    const std::vector<double>& SpecNoise_vec = workspace.GetSpecNoise();
    std::vector<double>::const_iterator i_spec_noise = SpecNoise_vec.begin();

    const std::vector<double>& TdNoise_vec = workspace.GetTdNoise();
    std::vector<double>::const_iterator i_td_noise = TdNoise_vec.begin();
    
    const std::vector<double>& MetabSNR_vec = workspace.GetMetabSNR();
    std::vector<double>::const_iterator i_metab_snr = MetabSNR_vec.begin();

    // SNR
    const pair_vec& snr = yfid.GetSNR();
    pair_vec::const_iterator i_snr = snr.begin();
    
    // optimiser results
    std::vector<std::vector < double > > info = workspace.GetLMinfo();

	integer M = basis.GetGroupMatrix().nsize();
	std::size_t nIdxBeta   = 2*M+1;
    
    fit = 0;
	for( coord_vec::const_iterator i_coord = options.GetFitList().begin(); i_coord != options.GetFitList().end(); ++i_coord )	
	{
        std::string stop;
        // find out the reason for stopping
        if( 1 == info[fit][6] )
            stop = "stopped by small gradient";
        else if( 2 == info[fit][6] )
            stop = "stopped by small D*p";
        else if( 3 == info[fit][6] )
            stop = "stopped by iteration limit";
        else if( 4 == info[fit][6] )
            stop = "stopped by singular matrix - restart with bigger mu";
        else if( 5 == info[fit][6] )
            stop = "no further reduction possible - restart with bigger mu";
        else if( 6 == info[fit][6] )
            stop = "stopped by small l2 norm of error";

		// write the voxel info
		fout << (*i_coord).row << "," << (*i_coord).col << "," << (*i_coord).slice << ","; 
        fout << *i_Q << "," << *i_Q_rel << "," << *i_metab_rat << "," << *i_peak_metab_rat << "," << *i_metab_fwhm << "," << *i_metab_fwhm * (yfid.GetTransmitterFrequency() / 1.0e6) << "," << (*i_snr).first << "," << *i_Q*(*i_snr).first << "," <<  *i_metab_snr << "," << *i_spec_noise << "," << *i_td_noise << "," << yfid.GetPPMRef(*i_coord) << "," << options.GetInitBetaUsed(fit) << "," << workspace.GetParas(fit)(nIdxBeta) << "," << pow(-workspace.GetParas(fit)(nIdxBeta)*log(0.5),0.5)*2.0/M_PI/(yfid.GetTransmitterFrequency()/1.0e6) << "," << yfid.GetPhi0(*i_coord)*180/M_PI << "," << -yfid.GetPhi1(*i_coord) * 180/M_PI * (yfid.GetTransmitterFrequency() / 1.0e6) * 2.0 * M_PI  << "," << workspace.GetNormalisationValue()[fit] << "," << workspace.GetWaterWidth(fit) << "," << workspace.GetWaterWidth(fit)/(yfid.GetTransmitterFrequency()/1.0e6) << "," << workspace.GetWaterFreq(fit) << "," << *i_BLV << ","<< *i_max_bl << "," << *i_min_bl << "," << *i_BLS << "," << info[fit][0] << "," << info[fit][1] << "," << stop << std::endl;

        // advance to the next voxel
        ++i_Q;
        ++i_Q_rel;
        ++i_metab_rat;
        ++i_peak_metab_rat;
        ++i_metab_fwhm;
        ++i_BLV;
        ++i_max_bl;
        ++i_min_bl;
        ++i_BLS;
        ++i_snr;
        ++i_spec_noise;
        ++i_td_noise;
        ++i_metab_snr;
        ++fit;
	}
    
    fout << "Basis shifts and dampings" << std::endl;
	fout << "Row" << "," << "Col" << "," << "Slice" << ","; 
	// output the row column names
	fout << "name,group,shift (PPM),damping (Hz)" << std::endl;

    fit = 0;
    for( coord_vec::const_iterator i_coord = options.GetFitList().begin(); i_coord != options.GetFitList().end(); ++i_coord )	
    {
        int basis_cnt = 1;
        int group_cnt = 0;
        for( integer c = 1; c <= M; c++ ) 
        {
            group_cnt++;
            std::size_t nIdxShift = c;
            std::size_t nIdxAlpha = c+M;
            // some output
            integer i = workspace.GetBasis().GetBasisFromGroup(c-1);
            if ( i != basis_cnt )
            {
                basis_cnt = i;
                group_cnt = 1;
            }
            fout << (*i_coord).row << "," << (*i_coord).col << "," << (*i_coord).slice << ","; 
            fout << signal_names[i-1] << "," << group_cnt << "," << -workspace.GetParas(fit)(nIdxShift) / yfid.GetTransmitterFrequency() * 1e6 << "," << workspace.GetParas(fit)(nIdxAlpha)/M_PI << std::endl;
        }
        fit++;
    }
    
    if ( workspace.GetDynShiftSize() > 0 )
    {
        fout << "Dynamic frequency corrections" << std::endl;
        fout << "index,shift (Hz)" << std::endl;
        for ( size_t n = 0; n < workspace.GetDynShiftSize(); n++ )
            fout << n+1 << "," << workspace.GetDynShift(n) << std::endl;
    }

    fout << "Geometry parameters" << std::endl;
    fout << "Rows," << yfid.GetRows() << std::endl;
    fout << "Cols," << yfid.GetCols() << std::endl;
    fout << "Slices," << yfid.GetSlices() << std::endl;

    if ( yfid.IsKnownVoxelDim() )
    {
        const std::vector<double>& voxel_dim = yfid.GetVoxelDim();
        fout << "PixelSpacing," << voxel_dim[0] << "\\" << voxel_dim[1] << std::endl;
        fout << "SliceThickness," << voxel_dim[2] << std::endl;
    }
    else
    {
        fout << "PixelSpacing,Unknown" << std::endl;
        fout << "SliceThickness,Unknown" << std::endl;
    }
    
    if ( yfid.IsKnownRowDirn() && yfid.IsKnownColDirn() )
    {
        fout << "ImageOrientationPatient,";
        std::string row_ori_str;
        rvec2str(yfid.GetRowDirn(),row_ori_str);
        std::string col_ori_str;
        rvec2str(yfid.GetColDirn(),col_ori_str);
        fout << row_ori_str + "\\" + col_ori_str << std::endl;
    }
    else
    {
        fout << "ImageOrientationPatient,Unknown" << std::endl;
    }


    if ( yfid.IsKnownPos() )
    {
        fout << "ImagePositionPatient,";
        std::string pos_str;
        rvec2str(yfid.GetPos(), pos_str);
        fout << pos_str << std::endl;
    }
    else
    {
        fout << "ImagePositionPatient,Unknown" << std::endl;
    }

    /*
    if ( yfid.IsKnownVoiDim() )
    {
        const std::vector<double>& voi_dim = yfid.GetVoiDim();
        fout << "VoiDimRow," << voi_dim[0] << std::endl;
        fout << "VoiDimCol," << voi_dim[1] << std::endl;
        fout << "VoiDimSlice," << voi_dim[2] << std::endl;
    }
    else
    {
        fout << "VoiDimRow,Unknown" << std::endl;
        fout << "VoiDimCol,Unknown" << std::endl;
        fout << "VoiDimSlice,Unknown" << std::endl;
    }
    */

    if ( yfid.IsKnownVoiDim() )
    {
        fout << "VoiDim,";
        std::string voi_str;
        rvec2str(yfid.GetVoiDim(), voi_str);
        fout << voi_str << std::endl;
    }
    else
    {
        fout << "VoiDim,Unknown" << std::endl;
    }

    fout << "Data parameters" << std::endl;
    fout << "Fs (Hz)," << yfid.GetSamplingFrequency() << std::endl;
    fout << "Ft (Hz)," << yfid.GetTransmitterFrequency() << std::endl;
    fout << "N," << yfid.GetNumberOfPoints() << std::endl;
    fout << "TE (s)," << yfid.GetEchoTime() << std::endl;

    std::string format_str;
    if( tarquin::NOTSET == options.GetFormat() ) 
        format_str = "not set";
    else if( tarquin::DANGER == options.GetFormat() ) 
        format_str = "dpt";
    else if( tarquin::SIEMENS == options.GetFormat() ) 
        format_str = "siemens";
    else if( tarquin::DCM == options.GetFormat() ) 
        format_str = "dcm";
    else if( tarquin::RDA == options.GetFormat() ) 
        format_str = "rda";
    else if( tarquin::PHILIPS == options.GetFormat() ) 
        format_str = "philips";
    else if( tarquin::PHILIPS_DCM == options.GetFormat() ) 
        format_str = "philips_dcm";
    else if( tarquin::GE == options.GetFormat() ) 
        format_str = "ge";
    else if( tarquin::SHF == options.GetFormat() ) 
        format_str = "shf";
    else if( tarquin::LCM == options.GetFormat() ) 
        format_str = "lcm";
    else if( tarquin::VARIAN == options.GetFormat() ) 
        format_str = "varian";
    else if( tarquin::BRUKER == options.GetFormat() ) 
        format_str = "bruker";
    else if( tarquin::JMRUI_TXT == options.GetFormat() ) 
        format_str = "jmrui_txt";

    fout << "Format," << format_str << std::endl;
    
    fout << std::endl;

    fout << "TARQUIN version," << version::version_string() << std::endl;

}

void ExportTxtResults(const std::string& strFilename, const Workspace& workspace)
{
    const CBasis& basis      = workspace.GetBasis();
	const Options& options  = workspace.GetOptions();
    const CFID& yfid = workspace.GetFID();
    const coord_vec& fit_list = options.GetFitList();

	rvec_stdvec ahat = workspace.GetAmplitudesNormalised();
	rvec_stdvec crlbnorm    = workspace.GetCRLBsNormalised();

	std::ofstream fout(strFilename.c_str(), std::ios::out);

	integer M = basis.GetGroupMatrix().nsize();
	std::size_t nIdxBeta   = 2*M+1;
	
    //for ( size_t fit = 0; fit < fit_list.size(); fit++ )

    size_t fit = 0;
	for ( coord_vec::const_iterator i_coord = fit_list.begin(); i_coord != fit_list.end(); i_coord++ )
    {
        // header line
        if ( fit_list.size() > 1 )
            fout << "Row : " << (*i_coord).row << ", Col : " << (*i_coord).col << ", Slice : " << (*i_coord).slice << std::endl << std::endl;

        size_t max_string_sz = 0;
        std::string temp_str;
        for( integer n = 1; n < ahat[fit].size(); n++ ) {
            temp_str = basis.GetSignalName(n);
            if ( temp_str.size() > max_string_sz )
                max_string_sz = temp_str.size();
        }

        if ( options.GetFilenameWater() == "" )
        {
            fout << std::setw(10) << std::left << "Signal";
            fout << std::setw(10) << std::right << "Conc (au)";
            fout << std::setw(10) << std::right << "%SD";
            fout << std::setw(10) << std::right << "SD (au)" << std::endl;
        }
        else
        {
            fout << std::setw(10) << std::left << "Signal";
            fout << std::setw(10) << std::right << "Conc (mM)";
            fout << std::setw(10) << std::right << "%SD";
            fout << std::setw(10) << std::right << "SD (mM)" << std::endl;
        }

        for(size_t m = 0; m < 40; m++)
            fout << "-";

        fout << std::endl;

        std::vector < std::string > signal_names = basis.GetSignalNames();	
        std::vector < std::string > signal_names_sorted = signal_names;

        std::sort(signal_names_sorted.begin(), signal_names_sorted.end()); 

        /*	
        // unsorted version kept for testing
        for(integer n = 1; n < ahat.size()+1; n++) {
        fout << signal_names[n-1]; 
        temp_str = signal_names[n-1];
        // pad with spaces
        for(size_t m = temp_str.size(); m < max_string_sz; m++)
        fout << " ";
        fout << "   " << ahat(n) << endl;
        }*/

        fout << std::showpoint;
        fout << std::setprecision(4);
        //fout << std::scientific;

        // sorted version
        double tempamp = 0;
        double tempcrlb = 0;
        for(integer n = 1; n < ahat[fit].size()+1; n++) {
            // find index of ahat now the list is sorted
            for(integer m = 1; m < ahat[fit].size()+1; m++) {
                if ( signal_names_sorted[n-1] == signal_names[m-1] ){
                    tempamp = ahat[fit](m);
                    tempcrlb = crlbnorm[fit](m);
                    break;
                }
            }
            temp_str = signal_names_sorted[n-1];
            fout << std::setw(10) << std::left << temp_str;
            fout << std::setw(10) << std::right << tempamp;
            fout << std::setw(10) << std::right << tempcrlb/tempamp*100;
            fout << std::setw(10) << std::right << tempcrlb << std::endl;
        }

        // combined metabolites
        
        if ( workspace.GetAmplitudesNormalisedComb().size() > 0 )
        {
            for(integer n = 1; n < workspace.GetAmplitudesNormalisedComb()[fit].size()+1; n++) {
                fout << std::setw(10) << std::left << workspace.GetMetabNamesComb()[n-1];
                fout << std::setw(10) << std::right << workspace.GetAmplitudesNormalisedComb()[fit][n];
                fout << std::setw(10) << std::right << workspace.GetCRLBsNormalisedComb()[fit][n]/workspace.GetAmplitudesNormalisedComb()[fit][n]*100;
                fout << std::setw(10) << std::right << workspace.GetCRLBsNormalisedComb()[fit][n] << std::endl;
            }
        }

        // output some other useful information
        fout << std::endl;
        fout << "Fit quality" << std::endl;
        fout << "-----------" << std::endl;
        const std::vector<double> metab_fwhm = workspace.GetMetabFWHM();
        fout << "Metab FWHM (PPM) : " << metab_fwhm[fit] << std::endl;
        fout << "Metab FWHM (Hz)  : " << metab_fwhm[fit]*(yfid.GetTransmitterFrequency()/1.0e6) << std::endl;
        const pair_vec snr = yfid.GetSNR();
        fout << "SNR residual     : " << snr[fit].first << std::endl;
        const std::vector<double> Q = workspace.GetQ();
        fout << "SNR max          : " << snr[fit].first * Q[fit] << std::endl;
        fout << "Q                : " << Q[fit] << std::endl;
        fout << std::endl;
        fout << "Fitting parameters" << std::endl;
        fout << "------------------" << std::endl;
        fout << "Initial beta : " << options.GetInitBetaUsed(fit) << std::endl;
        fout << "Final beta   : " << workspace.GetParas(fit)(nIdxBeta) << std::endl;
        //fout << "Initial phi0 : " << options.GetTypicalPhi0() << std::endl;
        //fout << "Final phi0   : " << std::endl;
        //fout << "Initial phi1 : " << options.GetTypicalPhi1() << std::endl;
        //fout << "Final phi1   : " << std::endl;
        fout << "Max iters    : " << options.GetMaxIters() << std::endl;
        fout << "Start point  : " << options.GetRangeStart() << std::endl;
        fout << "End point    : " << options.GetRangeEnd() << std::endl;
        const pair_vec ref = yfid.GetPPMRef();
        fout << "Ref (ppm)    : " << ref[fit].first << std::endl;
        fout << std::endl;
        fout << "Concentration scaling" << std::endl;
        fout << "---------------------" << std::endl;
        if ( options.GetFilenameWater() == "" )
            fout << "Water scaling    : false" << std::endl;
        else
            fout << "Water scaling    : true" << std::endl;
        fout << "Water amp        : " << workspace.GetNormalisationValue()[fit] << std::endl;
        fout << "Water conc       : " << options.GetWConc() << std::endl;
        fout << "Water att        : " << options.GetWAtt() << std::endl;
        if ( options.GetFilenameWater() == "" )
        {
            fout << "Water FWHM (Hz)  : NA" << std::endl; 
            fout << "Water FWHM (PPM) : NA" << std::endl; 
            fout << "Water freq (Hz)  : NA" << std::endl; 
        }
        else
        {
            fout << "Water FWHM (Hz)  : " << workspace.GetWaterWidth(fit) << std::endl; 
            fout << "Water FWHM (PPM) : " << workspace.GetWaterWidth(fit)/(yfid.GetTransmitterFrequency()/1.0e6) << std::endl; 
            fout << "Water freq (Hz)  : " << workspace.GetWaterFreq(fit) << std::endl; 
        }

        fout << std::endl;
        fout << "Data parameters" << std::endl;
        fout << "---------------" << std::endl;
        fout << "Fs (Hz) : " << yfid.GetSamplingFrequency() << std::endl;
        fout << "Ft (Hz) : " << yfid.GetTransmitterFrequency() << std::endl;
        fout << "N       : " << yfid.GetNumberOfPoints() << std::endl;
        fout << "TE (s)  : " << yfid.GetEchoTime() << std::endl;

        std::string format_str;
        if( tarquin::NOTSET == options.GetFormat() ) 
            format_str = "not set";
        else if( tarquin::DANGER == options.GetFormat() ) 
            format_str = "dpt";
        else if( tarquin::SIEMENS == options.GetFormat() ) 
            format_str = "siemens";
        else if( tarquin::DCM == options.GetFormat() ) 
            format_str = "dcm";
        else if( tarquin::RDA == options.GetFormat() ) 
            format_str = "rda";
        else if( tarquin::PHILIPS == options.GetFormat() ) 
            format_str = "philips";
        else if( tarquin::PHILIPS_DCM == options.GetFormat() ) 
            format_str = "philips_dcm";
        else if( tarquin::GE == options.GetFormat() ) 
            format_str = "ge";
        else if( tarquin::SHF == options.GetFormat() ) 
            format_str = "shf";
        else if( tarquin::LCM == options.GetFormat() ) 
            format_str = "lcm";
        else if( tarquin::VARIAN == options.GetFormat() ) 
            format_str = "varian";
        else if( tarquin::BRUKER == options.GetFormat() ) 
            format_str = "bruker";
        else if( tarquin::JMRUI_TXT == options.GetFormat() ) 
            format_str = "jmrui_txt";

        fout << "Format  : " << format_str << std::endl;

        std::vector<std::vector < double > > info = workspace.GetLMinfo();
        fout << std::endl;
        fout << "Optimisation details" << std::endl;
        fout << "--------------------" << std::endl;
        fout << "l2 norm of error at initial p  : " << info[fit][0] << std::endl;
        fout << "l2 norm of error at final p    : " << info[fit][1] << std::endl;
        fout << "l2 norm of J.'*e at final p    : " << info[fit][2] << std::endl;
        fout << "l2 norm of D*p at final p      : " << info[fit][3] << std::endl;
        fout << "number of iterations           : " << info[fit][5] << std::endl;
        fout << "number of function evaluations : " << info[fit][6] << std::endl;
        fout << "number of Jacobian evaluations : " << info[fit][7] << std::endl;
        if( 1 == info[fit][6] )
            fout << "stopped by small gradient" << std::endl;
        else if( 2 == info[fit][6] )
            fout << "stopped by small D*p" << std::endl;
        else if( 3 == info[fit][6] )
            fout << "stopped by iteration limit" << std::endl;
        else if( 4 == info[fit][6] )
            fout << "stopped by singular matrix - restart with bigger mu" << std::endl;
        else if( 5 == info[fit][6] )
            fout << "no further reduction possible - restart with bigger mu" << std::endl;
        else if( 6 == info[fit][6] )
            fout << "stopped by small l2 norm of error" << std::endl;
        fout << std::endl;
        fit++;
    }
    fout << "TARQUIN version " << version::version_string() << std::endl;
}

void GetTable(std::ostringstream& fout, const Workspace& workspace)
{
    const CBasis& basis      = workspace.GetBasis();
	const Options& options  = workspace.GetOptions();

	rvec_stdvec ahat = workspace.GetAmplitudesNormalised();
	rvec_stdvec crlbnorm    = workspace.GetCRLBsNormalised();

	size_t max_string_sz = 0;
	std::string temp_str;
	for( integer n = 1; n < ahat[0].size(); n++ ) { // TODO
		temp_str = basis.GetSignalName(n);
		if ( temp_str.size() > max_string_sz )
			max_string_sz = temp_str.size();
	}

	if ( options.GetFilenameWater() == "" )
    {
        fout << std::setw(10) << std::left << "Signal";
        fout << std::setw(10) << std::right << "Conc (au)";
        fout << std::setw(10) << std::right << "%SD";
        fout << std::setw(10) << std::right << "SD (au)";
        fout << std::setw(10) << "\\n";
    }
    else
    {
        fout << std::setw(10) << std::left << "Signal";
        fout << std::setw(10) << std::right << "Conc (mM)";
        fout << std::setw(10) << std::right << "%SD";
        fout << std::setw(10) << std::right << "SD (mM)";
        fout << std::setw(10) << "\\n";
    }
	
    for(size_t m = 0; m < 40; m++)
		fout << "-";

	fout << "\\n";

	std::vector < std::string > signal_names = basis.GetSignalNames();	
	std::vector < std::string > signal_names_sorted = signal_names;

	std::sort(signal_names_sorted.begin(), signal_names_sorted.end()); 

    fout << std::showpoint;
    fout << std::setprecision(3);
	//fout << std::scientific;

	// sorted version
	double tempamp = 0;
	double tempcrlb = 0;
	for(integer n = 1; n < ahat[0].size()+1; n++) { // TODO
		// find index of ahat now the list is sorted
		for(integer m = 1; m < ahat[0].size()+1; m++) { // TODO
			if ( signal_names_sorted[n-1] == signal_names[m-1] ){
				tempamp = ahat[0](m); // TODO
                tempcrlb = crlbnorm[0](m); // TODO
				break;
			}
		}
		temp_str = signal_names_sorted[n-1];
        fout << std::setw(10) << std::left << temp_str;
        fout << std::setw(10) << std::right << tempamp;
        fout << std::setw(10) << std::right << tempcrlb/tempamp*100;
        fout << std::setw(10) << std::right << tempcrlb << "\\n";
	}
        
    // combined metabolites
    if ( workspace.GetAmplitudesNormalisedComb().size() > 0 )
    {
        for(integer n = 1; n < workspace.GetAmplitudesNormalisedComb()[0].size()+1; n++) { // TODO
            fout << std::setw(10) << std::left << workspace.GetMetabNamesComb()[n-1];
            fout << std::setw(10) << std::right << workspace.GetAmplitudesNormalisedComb()[0][n]; // TODO
            fout << std::setw(10) << std::right << workspace.GetCRLBsNormalisedComb()[0][n]/workspace.GetAmplitudesNormalisedComb()[0][n]*100; //TODO
            fout << std::setw(10) << std::right << workspace.GetCRLBsNormalisedComb()[0][n] << "\\n"; // TODO
        }
    }


	CFID yfid = workspace.GetFID();
    
    // SNR
    const pair_vec& snr = yfid.GetSNR();

    // Q
    const std::vector<double>& Q_vec = workspace.GetQ();
    
    // FWHM
    const std::vector<double>& metab_fwhm_vec = workspace.GetMetabFWHM();
    
    // baseline var
    const std::vector<double>& BLV_vec = workspace.GetBLV();

	integer M = basis.GetGroupMatrix().nsize();
	std::size_t nIdxBeta   = 2*M+1;
    
    for(size_t m = 0; m < 40; m++)
		fout << "-";

    fout << "\\n";
    fout << "QC INFORMATION\\n";
    
    bool qc_state = true;
    fout << "Metab FWHM (PPM) = " << metab_fwhm_vec[0] <<  " : ";
    
    double metab_fwhm = metab_fwhm_vec[0];
    if ( (metab_fwhm > 0.1) || (metab_fwhm == -1) )
    {
        fout << "FAIL\\n";
        qc_state = false;
    }
    else if ( ( metab_fwhm <= 0.1 ) && ( metab_fwhm > 0.075 ) )
        fout << "PASS (borderline)\\n";
    else if ( ( metab_fwhm <= 0.075 ) && ( metab_fwhm > 0.05 ) )
        fout << "PASS (average)\\n";
    else if ( ( metab_fwhm <= 0.05 ) )
        fout << "PASS (good)\\n";

    double snr_qc = snr[0].first;
    fout << "SNR              = " << snr_qc << "   : " ;
    if (snr_qc < 4)
    {
        fout << "FAIL\\n";
        qc_state = false;
    }
    else if ( ( snr_qc >= 4 ) && ( snr_qc < 10 ) )
        fout << "PASS (borderline)\\n";
    else if ( ( snr_qc >= 10 ) && ( snr_qc < 20 ) )
        fout << "PASS (average)\\n";
    else if ( ( snr_qc >= 20 ) )
        fout << "PASS (good)\\n";

    if ( qc_state )
        fout << "Overall QC                : PASS\\n";
    else
        fout << "Overall QC                : FAIL\\n";
    
    for(size_t m = 0; m < 40; m++)
		fout << "-";
    fout << "\\n";

    fout << "DIAGNOSTICS\\n";
    //fout << "Metab FWHM (Hz)   = " << metab_fwhm_vec[0]*(yfid.GetTransmitterFrequency()/1.0e6) << "\\n";
    fout << "SNR max           = " << snr[0].first * Q_vec[0] << "\\n";
    //fout << "SNR residul       = " << snr[0].first << "\\n";
    fout << "Q                 = " << Q_vec[0] << "\\n";
    
    fout << "Water FWHM (PPM)  = ";
    if ( options.GetFilenameWater() == "" )
        fout << "NA\\n";
    else
        fout << workspace.GetWaterWidth(0)/(yfid.GetTransmitterFrequency()/1.0e6) << "\\n";

    /*fout << "Water FWHM (Hz)   = ";
    if ( options.GetFilenameWater() == "" )
        fout << "NA\\n";
    else
        fout << workspace.GetWaterWidth(0) << "\\n";
        */

    /*fout << "Water freq (Hz)   = ";
    if ( options.GetFilenameWater() == "" )
        fout << "NA\\n";
    else
        fout << workspace.GetWaterFreq(0) << "\\n";
        */

    fout << "Init beta         = " << options.GetInitBetaUsed(0) << "\\n";
    fout << "Final beta        = " << workspace.GetParas(0)(nIdxBeta) << "\\n";
    fout << "Final beta (PPM)  = " << pow(-workspace.GetParas(0)(nIdxBeta)*log(0.5),0.5)*2.0/M_PI/(yfid.GetTransmitterFrequency()/1.0e6) << "\\n";
    //fout << "Phi0 (deg)        = " << yfid.GetPhi0()[0].first * 180/M_PI << "\\n";
    //fout << "Phi1 (deg/PPM)    = " << -yfid.GetPhi1()[0].first * 180/M_PI * (yfid.GetTransmitterFrequency() / 1.0e6) * 2.0 * M_PI << "\\n";
    //fout << "Baseline var      = " << BLV_vec[0] << "\\n";
    fout << "Start point       = " << options.GetRangeStart() << "\\n";
    fout << "End point         = " << options.GetRangeEnd() << "\\n";
}

void ExportPdfResults(const std::string& strFilename, const Workspace& workspace, int fit_num)
{
    const Options& options = workspace.GetOptions();
    CFID yfid = workspace.GetFID();
	coord_vec fit_list = options.GetFitList();
    cvm::cvector y = yfid.GetVectorFID(fit_list[fit_num]);
    cvec_stdvec yhat_vec = workspace.GetSignalEstimate();
    cvm::cvector yhat = yhat_vec[fit_num];
	
    std::string title = options.GetTitle();
    
    // test code
    /*
    cvm::cvector linesh(yhat.size());
    treal frac;
    //for ( int n = 1; n < 500 + 1; n++ ) 
    for ( int n = 1; n < yhat.size() + 1; n++ ) 
    {
        frac = (treal) n/yhat.size();
        linesh(n) = y(n)/(yhat(n)*0.999+y(n)*0.001);
        //linesh(n) = y(n)/(yhat(n)*(1-frac*frac*frac)+y(n)*frac*frac*frac);
        //linesh(n) = y(n)/yhat(n);
        //linesh(n) = linesh(n).real();
    }
    
    cvm::cvector linesh_smo(yhat.size());
    td_conv_ws( linesh, linesh_smo, yfid.GetSamplingFrequency()/20, 10);	
    //plot(linesh);
    //plot(linesh_smo);

    cvm::cvector LIN(yhat.size());
    fft(linesh_smo, LIN);
    LIN = fftshift(LIN);
    //plot(LIN);
        
    for ( int n = 1; n < yhat.size() + 1; n++ ) 
        yhat(n) = yhat(n) * linesh_smo(n);
    */

    lb(y,workspace.GetFIDRaw(),options.Getlb());
    ZeroPad(y,options.GetZF());
    lb(yhat,workspace.GetFIDRaw(),options.Getlb());
    ZeroPad(yhat,options.GetZF());

    cvm::cvector Y(y.size());
    cvm::cvector YHAT(yhat.size());
    
    fft(y, Y);
    fft(yhat, YHAT);

    Y = fftshift(Y);
    YHAT = fftshift(YHAT);
    
    cvm::cvector residual = y-yhat;
    cvm::cvector baseline = residual;
  
    cvm::cvector RESIDUAL = fft(residual);
    RESIDUAL = fftshift(RESIDUAL);

    cvm::cvector BASELINE;
    td_conv_ws( RESIDUAL, BASELINE, options.GetBL()*options.GetZF(), 10);	

	const cmat_stdvec& s_vec = workspace.GetBasisMatrix(); // TODO should be a const reference?
	cvm::cmatrix s = s_vec[fit_num];
    
    cvm::cmatrix s_proc(s.msize()*options.GetZF(), s.nsize());

	rvec_stdvec ahat_vec = workspace.GetAmplitudes();
	cvm::rvector ahat = ahat_vec[fit_num];
    for ( int n = 1; (n < ahat.size() + 1); n++) 
    {
        cvm::cvector s_temp = s(n);
        lb(s_temp,workspace.GetFIDRaw(),options.Getlb());
        ZeroPad(s_temp,options.GetZF());
        s_proc(n) = s_temp*ahat(n);
    }

    cvm::cmatrix S = fft(s_proc);
    S = fftshift(S);

	coord voxel = fit_list[fit_num];
    cvm::rvector freq_scale = yfid.GetPPMScale(voxel, options.GetZF());

    cvm::cmatrix all(S.msize(),S.nsize()+7);
    all.assign(1,8,S);

    // find points corresponding to ppm start and ppm end
    int left = 1, right = 1;
    for ( int n = 1; n < (freq_scale.size()); n++ ) 
    {
        if ( ( freq_scale(n) > options.GetPPMend() ) && ( freq_scale(n+1) < options.GetPPMend() ) )
            left = n;
        if ( ( freq_scale(n) > options.GetPPMstart() ) && ( freq_scale(n+1) < options.GetPPMstart() ) )
            right = n;
    }

    //std::cout << freq_scale(left) << std::endl;
    //std::cout << freq_scale(right) << std::endl;
    
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

    /*double noise_min = std::numeric_limits<double>::infinity();
    double noise_temp = 0;
    int block_size = 200;
    for ( int n = 1; n < (RESIDUAL.size()+1) - block_size; n = n + block_size ) 
    {
        noise_temp = stdev(RESIDUAL.real(),n,n+block_size-1);
        if ( noise_temp < noise_min )
            noise_min = noise_temp;
    }

    std::cout << "stdev of res = " << noise_min << std::endl;
    std::cout << "stdev of fit res = " << stdev(RESIDUAL.real() - BASELINE.real(),left,right) << std::endl;
    std::cout << "Q = " << stdev(RESIDUAL.real() - BASELINE.real(),left,right) /noise_min << std::endl;*/

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
   
    // generate a table
    std::ostringstream table;
    GetTable(table, workspace);

	std::vector < std::string > signal_names = workspace.GetBasis().GetSignalNames();	

    bool ext_output = options.GetPdfExt();

    // and plot results	
    savepdffit(freq_scale, all, strFilename, table, title, signal_names, ext_output, options.GetPPMstart(), options.GetPPMend(), options.GetGnuplotCex());
}


