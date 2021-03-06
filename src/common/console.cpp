#include "console.hpp"
#include "version/version.h"
#include <boost/function.hpp>



void tarquin::DisplaySplash()
{
	std::cout << "\nTARQUIN Version " << version::version_string() << "\n";
	std::cout << "\n**** TARQUIN is a research tool and is NOT for clinical use. **** \n";
}

void tarquin::DisplayUsage()
{
	std::cout << "\nCommand line arguments are:";
	std::cout << "\n\t--help              displays this message";
	std::cout << "\n\t--para_file         read parameters from a specified txt file";
	std::cout << "\n\t--input             fidfile (water suppressed)";
	std::cout << "\n\t--input_w           fidfile (water reference)";
	std::cout << "\n\t--format            {siemens | philips | ge | dpt | rda | lcm | shf | varian | bruker | jmrui_txt}";
	//std::cout << "\n\t--ge_offset         offset_to_data";
	//std::cout << "\n\t--ge_wsframes       number of ws averages";
	std::cout << "\n\t--ge_wframes        number of w averages";
	//std::cout << "\n\t--ge_coils          number of coils";
	//std::cout << "\n\t--ge_samples        number of samples";
	std::cout << "\n\t--start_pnt         starting sample (1-based)";
	std::cout << "\n\t--end_pnt           ending sample (1-based)";
	std::cout << "\n\t--max_iters         maximum number of iterations to perform";
	std::cout << "\n\t--water_width       width of hsvd in hz";
	std::cout << "\n\t--water_width_df    width of downfield hsvd in hz (default is -inf)";
	std::cout << "\n\t--hsvd_comp         number of components in the hsvd decomp";
	std::cout << "\n\t--max_hsvd_pts      max number of pts used in the hsvd decomp";
	std::cout << "\n\t--conv_width        width of convolution window in points";
	std::cout << "\n\t--lipid_filter      {true | false} remove signals upfield of 1.8ppm in hsvd";
	std::cout << "\n\t--lipid_filter_freq Set ppm of lipid filter";
	std::cout << "\n\t--swap_row_col      {true | false} swap CSI rows and cols";
	std::cout << "\n\t--av_list           CSV file containing voxels to be averaged prior to fitting";
	std::cout << "\n\t--auto_phase        {true | false}";
	std::cout << "\n\t--auto_ref          {true | false}";
	std::cout << "\n\t--crlb_optim        calculate more optimistic CRLBs in output {true | false}";
	std::cout << "\n\t--max_dref          the max deviation from ref allowed by auto_ref";
	std::cout << "\n\t--max_phi0          the value of phi0 in rads";
	std::cout << "\n\t--max_phi1          the value of phi1_max/fs/2";
	std::cout << "\n\t--zfill_kspace      factor to zerofill the row, col direction of MRSI data";
	std::cout << "\n\t--filter_kspace     {true | false } apply a Hamming filter to MRSI k-space in row, col direction";
	std::cout << "\n\t--ref_freq          frequency of a single peak to be used for auto referencing (ppm)";
	std::cout << "\n\t--ref_file          CSV file containing reference peak list";
	std::cout << "\n\t--ref_signals       {1h_naa_cr_cho_lip | 1h_naa_cho | 1h_naa_cr_cho | 1h_cr_cho | 1h_naa | 1h_cr | 1h_cho | 1h_h2o | 1h_lip | 31p_pcr | 31p_pcr_gammaapt}";
	std::cout << "\n\t--dref_signals      reference signals for dynamic frequency correction, options as above";
	std::cout << "\n\t--fs                sampling frequency in Hz";
	std::cout << "\n\t--ft                transmitter frequency in Hz";
	std::cout << "\n\t--basis_csv         path to basis (CSV files)";
	std::cout << "\n\t--basis_xml         path to basis (precompiled XML files)";
	std::cout << "\n\t--basis_lcm         path to basis (LCModel .basis format)";
	std::cout << "\n\t--int_basis         {1h_brain | 1h_brain_exmm | h1_brain_glth | 1h_brain_gly_glth | 1h_brain_gly_cit_glth | 1h_brain_full | 1h_brain_le | 1h_brain_no_pcr | 1h_brain_metab_only | megapress_gaba | braino | 31p_brain}";
	std::cout << "\n\t--echo              echo time in seconds";
	std::cout << "\n\t--te1               te1 time in seconds for PRESS sequence";
	std::cout << "\n\t--tm                tm time in seconds for STEAM sequence";
	std::cout << "\n\t--acq_delay         acquistion delay time for pulse acquire seq in seconds";
	std::cout << "\n\t--cpmg_pulses       number of pulses for CPMG sequence";
	std::cout << "\n\t--pul_seq           {press | steam | slaser | laser | pulse_acq | cpmg | se | mega_press}";
	std::cout << "\n\t--dyn_av            {default | none | all | subtract | odd | even} water sup. dyn averaging scheme";
	std::cout << "\n\t--dyn_av_w          {default | none | all | subtract | odd | even} water dyn averaging scheme";
	std::cout << "\n\t--inv_even_paris    {true | false} invert even pairs of dynamic scans (useful for GE data in no_add mode_";
	std::cout << "\n\t--dyn_freq_corr     {true | false}";
	std::cout << "\n\t--output_dyn_shifts CSV file to save the dynamic shifts in Hz.";
	std::cout << "\n\t--pdfc              Pair-wise dynamic frequency correction {true | false}";
	std::cout << "\n\t--ref               reference offset in ppm";
	std::cout << "\n\t--water_eddy        {true | false}";
	std::cout << "\n\t--w_conc            NMR visable water concentration (35880)";
	std::cout << "\n\t--w_att             Water attenuation (0.7)";
	std::cout << "\n\t--basis_comp        {true | false}";
	std::cout << "\n\t--show_pre          {true | false}";
	std::cout << "\n\t--svs_only          {true | false}";
	std::cout << "\n\t--output_xml        output_file.xml";
	std::cout << "\n\t--view_fit          output_file.xml";
	std::cout << "\n\t--plot_sig          whitespace seperated list of signals";
	std::cout << "\n\t--gnuplot           path to gnuplot";
	std::cout << "\n\t--gnuplot_cex       gnuplot font expansion factor";
	std::cout << "\n\t--gnuplot_xtic      gnuplot xtic spacing (0.2ppm)";
	std::cout << "\n\t--output_basis      basis XML file";
	std::cout << "\n\t--output_basis_lcm  LCModel basis file";
	std::cout << "\n\t--output_image      plot of the fit in pdf format";
	std::cout << "\n\t--output_pdf        A4 pdf results page for printing";
	std::cout << "\n\t--ext_pdf           Extended output {true | false}";
    std::cout << "\n\t--stack_pdf         Stacked plot {true | false}";
	std::cout << "\n\t--title             title of results page";
	std::cout << "\n\t--output_txt        txt output of results";
	std::cout << "\n\t--output_csv        csv output of results";
	std::cout << "\n\t--output_fit        csv output of fit";
	std::cout << "\n\t--output_fit_m      csv output of fit (magnitude)";
	std::cout << "\n\t--output_spec       csv output of processed spectra";
	std::cout << "\n\t--output_spec_m     csv output of processed spectra (magnitude)";
	std::cout << "\n\t--write_raw         DPT fidfile of the raw data";
	std::cout << "\n\t--write_raw_w       DPT fidfile of the raw data";
	std::cout << "\n\t--write_raw_lcm     LCM fidfile of the raw data";
	std::cout << "\n\t--write_raw_lcm_w   LCM fidfile of the raw data";
	std::cout << "\n\t--write_pre         DPT fidfile of the preprocessed data";
	std::cout << "\n\t--write_post        DPT fidfile of the postprocessed data";

	std::cout << "\nExample:\n";
	std::cout << "\n\ttarquin --input MYFILE.DCM --format siemens --fs 1000 -ft 63e6 --basis_csv ../basis";
	std::cout << "\n\n";
}

namespace
{
	tarquin::fid_format_e parse_format_type(std::string strVal)
	{
		if( strVal == "siemens" ) 
			return tarquin::SIEMENS;

		else if( strVal == "philips" ) 
			return tarquin::PHILIPS;

		else if( strVal == "ge" ) 
			return tarquin::GE;

		else if( strVal == "dpt" ) 
			return tarquin::DANGER;

		else if( strVal == "rda" )
			return tarquin::RDA;

		else if( strVal == "lcm" )
			return tarquin::LCM;

		else if( strVal == "shf" )
			return tarquin::SHF;

		else if( strVal == "varian" )
			return tarquin::VARIAN;

		else if( strVal == "bruker" )
			return tarquin::BRUKER;

		else if( strVal == "philips_dcm" )
			return tarquin::PHILIPS_DCM;

		else if( strVal == "dcm" )
			return tarquin::DCM;

		else if( strVal == "jmrui_txt" )
			return tarquin::JMRUI_TXT;

		return tarquin::NOTSET;
	}

	bool parse_binary(std::string strVal)
    {
        if( strVal == "true" )
            return true;
        else
            return false;
    }
}

typedef boost::function<bool (std::string, std::string, tarquin::Options&)> handler_type;

bool lb_parser(std::string strKey, std::string strVal, tarquin::Options& options)
{
    std::istringstream iss(strVal, std::istringstream::in);
	iss >> options.m_lb;
    if( iss.fail() ) {
        std::cerr << "\nerror: couldn't recognise '" << strVal << "' as a number" << std::endl;
        return false;
    }
    return true;
}
bool lb_ref_parser(std::string strKey, std::string strVal, tarquin::Options& options)
{
    std::istringstream iss(strVal, std::istringstream::in);
	iss >> options.m_lb_ref;
    if( iss.fail() ) {
        std::cerr << "\nerror: couldn't recognise '" << strVal << "' as a number" << std::endl;
        return false;
    }
    return true;
}
bool init_beta_parser(std::string strKey, std::string strVal, tarquin::Options& options)
{
    std::istringstream iss(strVal, std::istringstream::in);
	iss >> options.m_init_beta;
    if( iss.fail() ) {
        std::cerr << "\nerror: couldn't recognise '" << strVal << "' as a number" << std::endl;
        return false;
    }
    return true;
}
bool max_beta_parser(std::string strKey, std::string strVal, tarquin::Options& options)
{
    std::istringstream iss(strVal, std::istringstream::in);
	iss >> options.m_max_beta;
    if( iss.fail() ) {
        std::cerr << "\nerror: couldn't recognise '" << strVal << "' as a number" << std::endl;
        return false;
    }
    return true;
}
bool conv_width_parser(std::string strKey, std::string strVal, tarquin::Options& options)
{
    std::istringstream iss(strVal, std::istringstream::in);
	iss >> options.m_conv_window_width;
    if( iss.fail() ) {
        std::cerr << "\nerror: couldn't recognise '" << strVal << "' as a number" << std::endl;
        return false;
    }
    return true;
}
bool water_width_parser(std::string strKey, std::string strVal, tarquin::Options& options)
{
    std::istringstream iss(strVal, std::istringstream::in);
	iss >> options.m_water_window;
    if( iss.fail() ) {
        std::cerr << "\nerror: couldn't recognise '" << strVal << "' as a number" << std::endl;
        return false;
    }
    return true;
}
bool start_pnt_parser(std::string strKey, std::string strVal, tarquin::Options& options)
{
    std::istringstream iss(strVal, std::istringstream::in);
	iss >> options.m_nStart;
    if( iss.fail() ) {
        std::cerr << "\nerror: couldn't recognise '" << strVal << "' as a number" << std::endl;
        return false;
    }
    return true;
}
bool end_pnt_parser(std::string strKey, std::string strVal, tarquin::Options& options)
{
    std::istringstream iss(strVal, std::istringstream::in);
	iss >> options.m_nEnd;
    if( iss.fail() ) {
        std::cerr << "\nerror: couldn't recognise '" << strVal << "' as a number" << std::endl;
        return false;
    }
    return true;
}
bool ppm_left_parser(std::string strKey, std::string strVal, tarquin::Options& options)
{
    std::istringstream iss(strVal, std::istringstream::in);
	iss >> options.m_ppm_end;
    if( iss.fail() ) {
        std::cerr << "\nerror: couldn't recognise '" << strVal << "' as a number" << std::endl;
        return false;
    }
    return true;
}
bool ppm_right_parser(std::string strKey, std::string strVal, tarquin::Options& options)
{
    std::istringstream iss(strVal, std::istringstream::in);
	iss >> options.m_ppm_start;
    if( iss.fail() ) {
        std::cerr << "\nerror: couldn't recognise '" << strVal << "' as a number" << std::endl;
        return false;
    }
    return true;
}
bool min_metab_alpha_parser(std::string strKey, std::string strVal, tarquin::Options& options)
{
    std::istringstream iss(strVal, std::istringstream::in);
	iss >> options.m_alpha_metab_lower;
    if( iss.fail() ) {
        std::cerr << "\nerror: couldn't recognise '" << strVal << "' as a number" << std::endl;
        return false;
    }
    return true;
}
bool max_metab_alpha_parser(std::string strKey, std::string strVal, tarquin::Options& options)
{
    std::istringstream iss(strVal, std::istringstream::in);
	iss >> options.m_alpha_metab_upper;
    if( iss.fail() ) {
        std::cerr << "\nerror: couldn't recognise '" << strVal << "' as a number" << std::endl;
        return false;
    }
    return true;
}
bool typ_metab_alpha_parser(std::string strKey, std::string strVal, tarquin::Options& options)
{
    std::istringstream iss(strVal, std::istringstream::in);
	iss >> options.m_alpha_metab_typ;
    if( iss.fail() ) {
        std::cerr << "\nerror: couldn't recognise '" << strVal << "' as a number" << std::endl;
        return false;
    }
    return true;
}
bool min_broad_alpha_parser(std::string strKey, std::string strVal, tarquin::Options& options)
{
    std::istringstream iss(strVal, std::istringstream::in);
	iss >> options.m_alpha_broad_lower;
    if( iss.fail() ) {
        std::cerr << "\nerror: couldn't recognise '" << strVal << "' as a number" << std::endl;
        return false;
    }
    return true;
}
bool max_broad_alpha_parser(std::string strKey, std::string strVal, tarquin::Options& options)
{
    std::istringstream iss(strVal, std::istringstream::in);
	iss >> options.m_alpha_broad_upper;
    if( iss.fail() ) {
        std::cerr << "\nerror: couldn't recognise '" << strVal << "' as a number" << std::endl;
        return false;
    }
    return true;
}
bool typ_broad_alpha_parser(std::string strKey, std::string strVal, tarquin::Options& options)
{
    std::istringstream iss(strVal, std::istringstream::in);
	iss >> options.m_alpha_broad_typ;
    if( iss.fail() ) {
        std::cerr << "\nerror: couldn't recognise '" << strVal << "' as a number" << std::endl;
        return false;
    }
    return true;
}

bool tarquin::ParseCommandLine(int argc, char* argv[], Options& options, CFID& fid)
{

	// create an empty vector of strings
	std::vector <std::string> args;
	// convert array of pointers to a vector of strings
	for(int n = 0; n < argc; n++ ) {
		args.push_back(argv[n]);
	}

	// search for a para_file before processing options
	for(size_t n = 1; n < args.size(); n+=2 ) 
	{
		std::string strKey = std::string(args[n]);

		if( strKey == "--help" ) {
			DisplayUsage();
			return true;
		}

		if( n == args.size()-1 ) {
			std::cerr << "\nerror: expected value after key '" << strKey << "'" << std::endl;
			return false;
		}

		std::string strVal = std::string(args[n+1]);
		if( strKey == "--para_file" ) {
			std::string line;
			std::ifstream para_file(strVal.c_str());
			if ( para_file.is_open() ) {
				while ( ! para_file.eof() )
				{
					getline(para_file,line);
					// strip leading and trailing whitespace from line
					size_t startpos = line.find_first_not_of(" \t");
					size_t endpos = line.find_last_not_of(" \t");
					if (( std::string::npos == startpos ) || ( std::string::npos == endpos))
					{
						line = "";
					}
					else
                    {
						line = line.substr(startpos, endpos-startpos+1);
                    }

					// split line into two parts dictated by the first found space char
					std::string::size_type pos = line.find(" ");
					// if not a comment, append line to the args vector
					if ( line.substr(0,1) != "#" ) {
						args.push_back("--"+line.substr(0,pos));
						args.push_back(line.substr(pos+1,line.size()-1));
					}

					//std::cout << "'--" + line.substr(0,pos) << "'" << std::endl;
					//std::cout << "'" << line.substr(pos+1,line.size()-1) << "'" << std::endl;
				}
				para_file.close();
			}
			else std::cerr << "\nUnable to open parameter file : " << strVal << "." << std::endl;
		}

	}



	// for each command line argument
	for(size_t n = 1; n < args.size(); n+=2 ) 
	{
		std::string strKey = std::string(args[n]);

		if( n == args.size()-1 ) {
			std::cerr << "\nerror: expected value after key '" << strKey << "'" << std::endl;
			return false;
		}

		std::string strVal = std::string(args[n+1]);

        //std::cout << strKey << std::endl;
        //std::cout << strVal << std::endl;

        bool handler_key_found = false;
    
        std::map<std::string, handler_type> handlers;
        handlers.insert( std::make_pair("--lb", &lb_parser) );
        handlers.insert( std::make_pair("--lb_ref", &lb_ref_parser) );
        handlers.insert( std::make_pair("--init_beta", &init_beta_parser) );
        handlers.insert( std::make_pair("--max_beta", &max_beta_parser) );
        handlers.insert( std::make_pair("--conv_width", &conv_width_parser) );
        handlers.insert( std::make_pair("--water_width", &water_width_parser) );
        handlers.insert( std::make_pair("--start_pnt", &start_pnt_parser) );
        handlers.insert( std::make_pair("--end_pnt", &end_pnt_parser) );
        handlers.insert( std::make_pair("--ppm_right", &ppm_right_parser) );
        handlers.insert( std::make_pair("--ppm_left", &ppm_left_parser) );
        handlers.insert( std::make_pair("--min_metab_alpha", &min_metab_alpha_parser) );
        handlers.insert( std::make_pair("--max_metab_alpha", &max_metab_alpha_parser) );
        handlers.insert( std::make_pair("--typ_metab_alpha", &typ_metab_alpha_parser) );
        handlers.insert( std::make_pair("--min_broad_alpha", &min_broad_alpha_parser) );
        handlers.insert( std::make_pair("--max_broad_alpha", &max_broad_alpha_parser) );
        handlers.insert( std::make_pair("--typ_broad_alpha", &typ_broad_alpha_parser) );

        for( std::map<std::string, handler_type>::iterator i = handlers.begin(); i != handlers.end(); ++i )
        {
            if( i->first == strKey )
            {
                bool ok = (i->second)(strKey, strVal, options);
                if (!ok)
                    return false;
                handler_key_found = true;
            }
        }

		// name of FID file
		if( strKey == "--input" ) 
        {
			options.m_strFile = strVal;
        }

		// name of water reference fid file
		else if( strKey == "--input_w" )
			options.m_strFileWater = strVal;

		// format of FID object
		else if( strKey == "--format" ) 
		{
			options.m_format = parse_format_type(strVal);

			if( options.m_format == NOTSET )
			{
				std::cerr << "\nerror: unrecognised format '" << strVal << "'" << std::endl;
				return false;
			}

		}

		// format of FID object
		else if( strKey == "--ref_signals" ) 
		{
			if( strVal == "1h_naa_cr_cho_lip" ) 
				options.m_ref_signals = tarquin::PROTON_NAA_CR_CHO_LIP;

			else if( strVal == "1h_naa_cr_cho" ) 
				options.m_ref_signals = tarquin::PROTON_NAA_CR_CHO;

			else if( strVal == "1h_naa_cho" ) 
				options.m_ref_signals = tarquin::PROTON_NAA_CHO;
            
            else if( strVal == "1h_cr_cho" ) 
				options.m_ref_signals = tarquin::PROTON_CR_CHO;

			else if( strVal == "1h_h2o" ) 
				options.m_ref_signals = tarquin::PROTON_H2O;

			else if( strVal == "1h_lip" ) 
				options.m_ref_signals = tarquin::PROTON_LIP;

			else if( strVal == "1h_naa" ) 
				options.m_ref_signals = tarquin::PROTON_NAA;

			else if( strVal == "1h_cr" ) 
				options.m_ref_signals = tarquin::PROTON_CR;

			else if( strVal == "1h_cho" ) 
				options.m_ref_signals = tarquin::PROTON_CHO;
			
            else if( strVal == "31p_pcr" ) 
				options.m_ref_signals = tarquin::PHOSPH_PCR;

            else if( strVal == "31p_pcr_gammaatp" ) 
				options.m_ref_signals = tarquin::PHOSPH_PCR_GAMMAATP;

			else 
			{
				std::cerr << "\nerror: unrecognised ref mode '" << strVal << "'" << std::endl;
				return false;
			}
        }

		else if( strKey == "--dref_signals" ) 
		{
			if( strVal == "1h_naa_cr_cho_lip" ) 
				options.m_dyn_ref_signals = tarquin::PROTON_NAA_CR_CHO_LIP;

			else if( strVal == "1h_naa_cr_cho" ) 
				options.m_dyn_ref_signals = tarquin::PROTON_NAA_CR_CHO;

			else if( strVal == "1h_naa_cho" ) 
				options.m_dyn_ref_signals = tarquin::PROTON_NAA_CHO;
            
            else if( strVal == "1h_cr_cho" ) 
				options.m_dyn_ref_signals = tarquin::PROTON_CR_CHO;

			else if( strVal == "1h_h2o" ) 
				options.m_dyn_ref_signals = tarquin::PROTON_H2O;
            
            else if( strVal == "1h_lip" ) 
				options.m_dyn_ref_signals = tarquin::PROTON_LIP;

			else if( strVal == "1h_naa" ) 
				options.m_dyn_ref_signals = tarquin::PROTON_NAA;

			else if( strVal == "1h_cr" ) 
				options.m_dyn_ref_signals = tarquin::PROTON_CR;

			else if( strVal == "1h_cho" ) 
				options.m_dyn_ref_signals = tarquin::PROTON_CHO;
			
            else if( strVal == "31p_pcr" ) 
				options.m_dyn_ref_signals = tarquin::PHOSPH_PCR;

            else if( strVal == "31p_pcr_gammaatp" ) 
				options.m_dyn_ref_signals = tarquin::PHOSPH_PCR_GAMMAATP;

			else {
				std::cerr << "\nerror: unrecognised ref mode '" << strVal << "'" << std::endl;
				return false;
			}
        }

		// pulse seq name
		// default to press	
		else if( strKey == "--pul_seq" ) 
		{
			if( strVal == "press" ) 
				options.SetPulSeq(PRESS);
			else if( strVal == "steam" ) 
				options.SetPulSeq(STEAM);
			else if( strVal == "slaser" ) 
				options.SetPulSeq(SEMI_LASER);
			else if( strVal == "laser" ) 
				options.SetPulSeq(LASER);
			else if( strVal == "pulse_acq" ) 
				options.SetPulSeq(PULSE_ACQUIRE);
			else if( strVal == "cpmg" ) 
				options.SetPulSeq(CPMG);
            else if( strVal == "se" ) 
				options.SetPulSeq(SPIN_ECHO);
			else if( strVal == "mega_press" ) 
				options.SetPulSeq(MEGA_PRESS);
            else if( strVal == "profile_press" ) 
				options.SetPulSeq(PROFILE_PRESS);
            else if( strVal == "shaped_press" ) 
				options.SetPulSeq(SHAPED_PRESS);
			else {
				std::cerr << "\nerror: unrecognised pulse sequence '" << strVal << "'" << std::endl;
				return false;
			}
		}

		// dynamic averaging
		// default to none
		else if( strKey == "--dyn_av" ) 
		{
			if( strVal == "default" ) 
				options.SetDynAv(DEFAULT);
			else if ( strVal == "none" ) 
				options.SetDynAv(NONE);
			else if( strVal == "all" ) 
				options.SetDynAv(ALL);
			else if( strVal == "subtract" ) 
				options.SetDynAv(SUBTRACT);
			else if( strVal == "odd" ) 
				options.SetDynAv(ODD);
			else if( strVal == "even" ) 
				options.SetDynAv(EVEN);
			else {
				std::cerr << "\nerror: unrecognised dynamic averaging scheme '" << strVal << "'" << std::endl;
				return false;
			}
		}

		else if( strKey == "--dyn_av_w" ) 
		{
			if( strVal == "default" ) 
				options.SetDynAvW(DEFAULT);
			else if ( strVal == "none" ) 
				options.SetDynAvW(NONE);
			else if( strVal == "all" ) 
				options.SetDynAvW(ALL);
			else if( strVal == "subtract" ) 
				options.SetDynAvW(SUBTRACT);
			else if( strVal == "odd" ) 
				options.SetDynAvW(ODD);
			else if( strVal == "even" ) 
				options.SetDynAvW(EVEN);
			else {
				std::cerr << "\nerror: unrecognised dynamic averaging scheme '" << strVal << "'" << std::endl;
				return false;
			}
		}
        
        // center the water peak before water removal?
        else if( strKey == "--pre_ws_shift" )
            options.m_pre_ws_shift = parse_binary(strVal);

        // keep the above shift?
        else if( strKey == "--keep_pre_ws_shift" )
            options.m_keep_pre_ws_shift = parse_binary(strVal);

        // pre dyn av freq correction
        else if( strKey == "--dyn_freq_corr" )
            options.m_dyn_freq_corr = parse_binary(strVal);

		// name of dynamic shift file
		else if( strKey == "--output_dyn_shifts" ) {
			// store the path for loading later
			options.m_dyn_shift_file = strVal;
		}

        // pre dyn av freq correction
        else if( strKey == "--rescale_lcm_basis" )
            options.m_rescale_lcm_basis = parse_binary(strVal);

        else if( strKey == "--inv_even_pairs" )
            options.m_invert_even_pairs = parse_binary(strVal);

        // pairwise dyn av freq correction
        else if( strKey == "--pdfc" )
            options.m_pdfc = parse_binary(strVal);

        // td or fd for noise estimate
        else if( strKey == "--crlb_td" )
            options.m_crlb_td = parse_binary(strVal);
        
        // optimistic CRLBs
        else if( strKey == "--crlb_optim" )
            options.m_crlb_optim = parse_binary(strVal);

        // allow negative amps
        else if( strKey == "--nnls" )
            options.m_nnls = parse_binary(strVal);

        // soft cons
        else if( strKey == "--soft_cons" )
            options.m_soft_cons = parse_binary(strVal);

        // lineshape correction
        //else if( strKey == "--ls_corr" )
		//	options.m_lineshape_corr = parse_binary(strVal);

        else if( strKey == "--replace_fp" )
            options.m_replace_fp = parse_binary(strVal);

        else if( strKey == "--prepend_pts" ) {
            int pts;
			std::istringstream iss(strVal, std::istringstream::in);
			iss >> pts;
			if( iss.fail() ) {
				std::cerr << "\nerror: couldn't recognise '" << strVal << "' as a number" << std::endl;
				return false;
			}
			options.m_prepend_pts = pts;
        }

        else if( strKey == "--trunc_pts" ) {
            int pts;
			std::istringstream iss(strVal, std::istringstream::in);
			iss >> pts;
			if( iss.fail() ) {
				std::cerr << "\nerror: couldn't recognise '" << strVal << "' as a number" << std::endl;
				return false;
			}
			options.m_truncate_pts = pts;
        }
        
        // interal basis
		// default to none
		else if( strKey == "--int_basis" ) 
		{
			if( strVal == "1h_brain" ) 
				options.SetIntBasisSet(PROTON_BRAIN);
			else if( strVal == "1h_brain_exmm" ) 
				options.SetIntBasisSet(PROTON_BRAIN_MMEXP);
			else if( strVal == "1h_brain_metab_only" ) 
				options.SetIntBasisSet(PROTON_BRAIN_METAB_ONLY);
			else if( strVal == "1h_brain_glth" ) 
				options.SetIntBasisSet(PROTON_BRAIN_GLTH);
			else if( strVal == "1h_brain_gly_glth" ) 
				options.SetIntBasisSet(PROTON_BRAIN_GLY_GLTH);
			else if( strVal == "1h_brain_gly_cit_glth" ) 
				options.SetIntBasisSet(PROTON_BRAIN_GLY_CIT_GLTH);
			else if( strVal == "1h_brain_full" ) 
				options.SetIntBasisSet(PROTON_BRAIN_FULL);
			else if( strVal == "1h_brain_le" ) 
				options.SetIntBasisSet(PROTON_BRAIN_LE);
			else if( strVal == "1h_brain_no_pcr" ) 
				options.SetIntBasisSet(PROTON_BRAIN_NO_PCR);
			else if( strVal == "1h_brain_lcm" ) 
				options.SetIntBasisSet(PROTON_BRAIN_LCM);
			else if( strVal == "megapress_gaba" ) 
				options.SetIntBasisSet(PROTON_MEGAPRESS_GABA);
			else if( strVal == "braino" ) 
				options.SetIntBasisSet(PROTON_BRAINO);
			else if( strVal == "31p_brain" ) 
				options.SetIntBasisSet(PHOSPH_BRAIN_DECOUP);
			else {
				std::cerr << "\nerror: unrecognised internal basis '" << strVal << "'" << std::endl;
				return false;
			}
		}


		// basis directory parameter 
		else if( strKey == "--output_xml" ) {
			// store the path for loading later
			options.m_strOutputXMLPath = strVal;
		}

		// basis directory parameter 
		else if( strKey == "--basis_csv" || strKey == "--basis_xml" || strKey == "--basis_lcm" ) {
			// make sure that we have only created one kind of basis
			if( options.m_strBasisPath.size() != 0 ) {
				std::cerr << "\nerror: you can only specify one of {--basis_csv, --basis_xml, basis_lcm}" <<
					std::endl;

				return false;
			}

			// store the path for loading later
			options.m_strBasisPath = strVal;

			// a single precompiled XML file?
			if( strKey == "--basis_xml" || strKey == "--basis_lcm" )
				options.m_bUsePrecompiled = true;
			// no, CSV files, do we need to append a trailing /?
			else {
				if( options.m_strBasisPath[options.m_strBasisPath.size()-1] != filesep[0] )
					options.m_strBasisPath += filesep;
			}
		}

		// echo time parameter
		else if( strKey == "--echo" ) {

			// convert string to number and check for errors
			treal echo_time;
			std::istringstream iss(strVal, std::istringstream::in);
			iss >> echo_time;
			//
			if( iss.fail() ) {
				std::cerr << "\nerror: couldn't recognise '" << strVal << "' as a number" << std::endl;
				return false;
			}
			fid.SetEchoTime(echo_time);
		}
    
        // te1 press parameter
		else if( strKey == "--te1" ) {

			// convert string to number and check for errors
			treal te1;
			std::istringstream iss(strVal, std::istringstream::in);
			iss >> te1;
			//
			if( iss.fail() ) {
				std::cerr << "\nerror: couldn't recognise '" << strVal << "' as a number" << std::endl;
				return false;
			}
			options.SetPRESS_TE1(te1);
		}

        // tm steam parameter
		else if( strKey == "--tm" ) {

			// convert string to number and check for errors
			treal tm;
			std::istringstream iss(strVal, std::istringstream::in);
			iss >> tm;
			//
			if( iss.fail() ) {
				std::cerr << "\nerror: couldn't recognise '" << strVal << "' as a number" << std::endl;
				return false;
			}
			options.SetSTEAM_TM(tm);
		}

        // tm steam parameter
		else if( strKey == "--acq_delay" ) {

			// convert string to number and check for errors
			treal acq_delay;
			std::istringstream iss(strVal, std::istringstream::in);
			iss >> acq_delay;
			//
			if( iss.fail() ) {
				std::cerr << "\nerror: couldn't recognise '" << strVal << "' as a number" << std::endl;
				return false;
			}
			options.m_acq_delay = acq_delay;
		}
        
        // tm steam parameter
		else if( strKey == "--cpmg_pulses" ) {

			// convert string to number and check for errors
			int cpmg_pulses;
			std::istringstream iss(strVal, std::istringstream::in);
			iss >> cpmg_pulses;
			//
			if( iss.fail() ) {
				std::cerr << "\nerror: couldn't recognise '" << strVal << "' as a number" << std::endl;
				return false;
			}
			options.SetCPMG_N(cpmg_pulses);
		}

        // gnuplot font expansion
		else if( strKey == "--gnuplot_cex" ) {

			// convert string to number and check for errors
			treal cex;
			std::istringstream iss(strVal, std::istringstream::in);
			iss >> cex;
			//
			if( iss.fail() ) {
				std::cerr << "\nerror: couldn't recognise '" << strVal << "' as a number" << std::endl;
				return false;
			}
			options.m_gnuplot_cex = cex;
		}

        // gnuplot font expansion
		else if( strKey == "--gnuplot_xtic" ) {

			// convert string to number and check for errors
			treal cex;
			std::istringstream iss(strVal, std::istringstream::in);
			iss >> cex;
			//
			if( iss.fail() ) {
				std::cerr << "\nerror: couldn't recognise '" << strVal << "' as a number" << std::endl;
				return false;
			}
			options.m_gnuplot_xtic = cex;
		}

		/*else if( strKey == "--ppm_right" ) {
			// convert string to number and check for errors
			treal temp;
			std::istringstream iss(strVal, std::istringstream::in);
			iss >> temp;
			//
			if( iss.fail() ) {
				std::cerr << "\nerror: couldn't recognise '" << strVal << "' as a number" << std::endl;
				return false;
			}
			options.m_ppm_start = temp;
		}

		else if( strKey == "--ppm_left" ) {
			// convert string to number and check for errors
			treal temp;
			std::istringstream iss(strVal, std::istringstream::in);
			iss >> temp;
			//
			if( iss.fail() ) {
				std::cerr << "\nerror: couldn't recognise '" << strVal << "' as a number" << std::endl;
				return false;
			}
			options.m_ppm_end = temp;
		}*/
        
        /*else if( strKey == "--min_metab_alpha" ) {
			// convert string to number and check for errors
			treal temp;
			std::istringstream iss(strVal, std::istringstream::in);
			iss >> temp;
			//
			if( iss.fail() ) {
				std::cerr << "\nerror: couldn't recognise '" << strVal << "' as a number" << std::endl;
				return false;
			}
			options.m_alpha_metab_lower = temp;
		}

        else if( strKey == "--max_metab_alpha" ) {
			// convert string to number and check for errors
			treal temp;
			std::istringstream iss(strVal, std::istringstream::in);
			iss >> temp;
			//
			if( iss.fail() ) {
				std::cerr << "\nerror: couldn't recognise '" << strVal << "' as a number" << std::endl;
				return false;
			}
			options.m_alpha_metab_upper = temp;
		}*/
        /*
        else if( strKey == "--typ_metab_alpha" ) {
			// convert string to number and check for errors
			treal temp;
			std::istringstream iss(strVal, std::istringstream::in);
			iss >> temp;
			//
			if( iss.fail() ) {
				std::cerr << "\nerror: couldn't recognise '" << strVal << "' as a number" << std::endl;
				return false;
			}
			options.m_alpha_metab_typ = temp;
		}

        else if( strKey == "--min_broad_alpha" ) {
			// convert string to number and check for errors
			treal temp;
			std::istringstream iss(strVal, std::istringstream::in);
			iss >> temp;
			//
			if( iss.fail() ) {
				std::cerr << "\nerror: couldn't recognise '" << strVal << "' as a number" << std::endl;
				return false;
			}
			options.m_alpha_broad_lower = temp;
		}

        else if( strKey == "--max_broad_alpha" ) {
			// convert string to number and check for errors
			treal temp;
			std::istringstream iss(strVal, std::istringstream::in);
			iss >> temp;
			//
			if( iss.fail() ) {
				std::cerr << "\nerror: couldn't recognise '" << strVal << "' as a number" << std::endl;
				return false;
			}
			options.m_alpha_broad_upper = temp;
		}

        else if( strKey == "--typ_broad_alpha" ) {
			// convert string to number and check for errors
			treal temp;
			std::istringstream iss(strVal, std::istringstream::in);
			iss >> temp;
			//
			if( iss.fail() ) {
				std::cerr << "\nerror: couldn't recognise '" << strVal << "' as a number" << std::endl;
				return false;
			}
			options.m_alpha_broad_typ = temp;
		} */

        /*else if( strKey == "--beta_scale" ) {
			// convert string to number and check for errors
			treal temp;
			std::istringstream iss(strVal, std::istringstream::in);
			iss >> temp;
			//
			if( iss.fail() ) {
				std::cerr << "\nerror: couldn't recognise '" << strVal << "' as a number" << std::endl;
				return false;
			}
			options.m_beta_scale = temp;
		}*/

        else if( strKey == "--max_phi1" ) {
			// convert string to number and check for errors
			treal temp;
			std::istringstream iss(strVal, std::istringstream::in);
			iss >> temp;
			//
			if( iss.fail() ) {
				std::cerr << "\nerror: couldn't recognise '" << strVal << "' as a number" << std::endl;
				return false;
			}
			options.m_phi1_lower = -temp;
			options.m_phi1_upper = +temp;
		}

        else if( strKey == "--max_phi0" ) {
			// convert string to number and check for errors
			treal temp;
			std::istringstream iss(strVal, std::istringstream::in);
			iss >> temp;
			//
			if( iss.fail() ) {
				std::cerr << "\nerror: couldn't recognise '" << strVal << "' as a number" << std::endl;
				return false;
			}
			options.m_phi0_lower = -temp;
			options.m_phi0_upper = +temp;
		}

		// sampling frequency parameter
		else if( strKey == "--fs") {

			// convert string to number and check for errors
			treal fs;
			std::istringstream iss(strVal, std::istringstream::in);
			iss >> fs;

			if( iss.fail() ) {
				std::cerr << "\nerror: couldn't recognise '" << strVal << "' as a number" << std::endl;
				return false;
			}

			fid.SetSamplingFrequency(fs);
		}

		// transmitter frequency parameter
		else if( strKey == "--ft") {

			// convert string to number and check for errors
			treal ft;
			std::istringstream iss(strVal, std::istringstream::in);
			iss >> ft;

			if( iss.fail() ) {
				std::cerr << "\nerror: couldn't recognise '" << strVal << "' as a number" << std::endl;
				return false;
			}

			fid.SetTransmitterFrequency(ft);
		}

		// ending point of fid
        /*
		else if( strKey == "--end_pnt" ) {

			// convert string to number and check for errors
			std::istringstream iss(strVal, std::istringstream::in);
			iss >> options.m_nEnd;

			if( iss.fail() ) {
				std::cerr << "\nerror: couldn't recognise '" << strVal << "' as a number" << std::endl;
				return false;
			}
		}

		// starting point of fid
		else if( strKey == "--start_pnt" ) {
			// convert string to number and check for errors
			std::istringstream iss(strVal, std::istringstream::in);
			iss >> options.m_nStart;

			if( iss.fail() ) {
				std::cerr << "\nerror: couldn't recognise '" << strVal << "' as a number" << std::endl;
				return false;
			}
		}*/

		// min starting point of fid
		else if( strKey == "--min_start_pnt" ) {
			// convert string to number and check for errors
			std::istringstream iss(strVal, std::istringstream::in);
			iss >> options.m_min_range;

			if( iss.fail() ) {
				std::cerr << "\nerror: couldn't recognise '" << strVal << "' as a number" << std::endl;
				return false;
			}
		}

		else if( strKey == "--start_time" ) {
			// convert string to number and check for errors
			std::istringstream iss(strVal, std::istringstream::in);
			iss >> options.m_nStartTime;

			if( iss.fail() ) {
				std::cerr << "\nerror: couldn't recognise '" << strVal << "' as a number" << std::endl;
				return false;
			}
		}
        
        else if( strKey == "--threads" ) {
			// convert string to number and check for errors
			std::istringstream iss(strVal, std::istringstream::in);
			iss >> options.m_threads;

			if( iss.fail() ) {
				std::cerr << "\nerror: couldn't recognise '" << strVal << "' as a number" << std::endl;
				return false;
			}
		}

		// adaptive starting point of fid
		/*else if( strKey == "--adapt_sp" )
			options.m_bAdaptSp = parse_binary(strVal);*/

		// pre-hsvd
		else if( strKey == "--pre_hsvd" )
			options.m_pre_hsvd = parse_binary(strVal);

		// old phi1 limits
		else if( strKey == "--old_phase" )
			options.m_old_phase = parse_binary(strVal);

		// adaptive ending point of fid
		/*else if( strKey == "--adapt_ep" )
            options.m_bAdaptEp = parse_binary(strVal);*/

        // CSI rows to be fit
		else if( strKey == "--rows" ) {
			// convert string to number and check for errors
			std::istringstream iss(strVal, std::istringstream::in);
			iss >> options.m_fit_rows;

			if( iss.fail() ) {
				std::cerr << "\nerror: couldn't recognise '" << strVal << "' as a number" << std::endl;
				return false;
			}
		}

        // CSI cols to be fit
		else if( strKey == "--cols" ) {
			// convert string to number and check for errors
			std::istringstream iss(strVal, std::istringstream::in);
			iss >> options.m_fit_cols;

			if( iss.fail() ) {
				std::cerr << "\nerror: couldn't recognise '" << strVal << "' as a number" << std::endl;
				return false;
			}
		}
        
        // CSI slices to be fit
		else if( strKey == "--slices" ) {
			// convert string to number and check for errors
			std::istringstream iss(strVal, std::istringstream::in);
			iss >> options.m_fit_slices;

			if( iss.fail() ) {
				std::cerr << "\nerror: couldn't recognise '" << strVal << "' as a number" << std::endl;
				return false;
			}
		}

		// maximum number of iterations to perform
		else if( strKey == "--max_iters" ) {

			// convert string to number and check for errors
			std::istringstream iss(strVal, std::istringstream::in);
			iss >> options.m_max_iters;

			if( iss.fail() ) {
				std::cerr << "\nerror: couldn't recognise '" << strVal << "' as a number" << std::endl;
				return false;
			}
		}

		else if( strKey == "--zfill" ) {

			// convert string to number and check for errors
			std::istringstream iss(strVal, std::istringstream::in);
			iss >> options.m_zero_fill;

			if( iss.fail() ) {
				std::cerr << "\nerror: couldn't recognise '" << strVal << "' as a number" << std::endl;
				return false;
			}
		}
        
        else if( strKey == "--bl" ) {

			// convert string to number and check for errors
			std::istringstream iss(strVal, std::istringstream::in);
			iss >> options.m_baseline;

			if( iss.fail() ) {
				std::cerr << "\nerror: couldn't recognise '" << strVal << "' as a number" << std::endl;
				return false;
			}
		}

		else if( strKey == "--zfill_kspace" ) {

			// convert string to number and check for errors
			std::istringstream iss(strVal, std::istringstream::in);
			iss >> options.m_zfill_kspace;

			if( iss.fail() ) {
				std::cerr << "\nerror: couldn't recognise '" << strVal << "' as a number" << std::endl;
				return false;
			}
		}

		else if( strKey == "--filter_kspace" )
            options.m_au_norm  = parse_binary(strVal);


		// gnuplot path
		else if( strKey == "--gnuplot" ) {

			g_strGnuPlot = strVal;
		}

		// path for fit plot
		else if( strKey == "--output_image" ) {

			options.m_strFileOutImag = strVal;
		}
		
        // path for fit plot
		else if( strKey == "--output_pdf" ) {

			options.m_strFileOutPdfImag = strVal;
		}

		// path for txt results 
		else if( strKey == "--output_txt" ) {

			options.m_strFileOutTxt = strVal;
		}

		// path for csv results 
		else if( strKey == "--output_csv" ) {

			options.m_strFileOutCSV = strVal;
		}

		// path for csv results 
		else if( strKey == "--output_csv_geom" ) {

			options.m_strFileOutCSVGeom = strVal;
		}
        
        // path for csv results 
		else if( strKey == "--output_fit" ) {

			options.m_strFileOutCSVFit = strVal;
		}

        // path for csv results mag mode
		else if( strKey == "--output_fit_m" ) {

			options.m_strFileOutCSVFitMag = strVal;
		}

		else if( strKey == "--output_spec" ) {

			options.m_strFileOutCSVSpectraAligned = strVal;
		}

		else if( strKey == "--output_spec_m" ) {

			options.m_strFileOutCSVSpectraAlignedMag = strVal;
		}

		// water concentraion 
		else if( strKey == "--w_conc" ) {

			// convert string to number and check for errors
			std::istringstream iss(strVal, std::istringstream::in);
			iss >> options.m_w_conc;

			if( iss.fail() ) {
				std::cerr << "\nerror: couldn't recognise '" << strVal << "' as a number" << std::endl;
				return false;
			}
		}

		// water attenuation
		else if( strKey == "--w_att" ) {

			// convert string to number and check for errors
			std::istringstream iss(strVal, std::istringstream::in);
			iss >> options.m_w_att;

			if( iss.fail() ) {
				std::cerr << "\nerror: couldn't recognise '" << strVal << "' as a number" << std::endl;
				return false;
			}
		}
        
        // Do the au scaling for SVS data?
		else if( strKey == "--au_norm" ) 
            options.m_au_norm  = parse_binary(strVal);

		// water removal window width
		/*else if( strKey == "--water_width" ) {

			// convert string to number and check for errors
			std::istringstream iss(strVal, std::istringstream::in);
			iss >> options.m_water_window;

			if( iss.fail() ) {
				std::cerr << "\nerror: couldn't recognise '" << strVal << "' as a number" << std::endl;
				return false;
			}
		}*/

		// downfield water removal window width
		else if( strKey == "--water_width_df" ) {

			// convert string to number and check for errors
			std::istringstream iss(strVal, std::istringstream::in);
			iss >> options.m_df_water_window;

			if( iss.fail() ) {
				std::cerr << "\nerror: couldn't recognise '" << strVal << "' as a number" << std::endl;
				return false;
			}
		}
		
        // number of hsvd comps
		else if( strKey == "--hsvd_comps" ) {

			// convert string to number and check for errors
			std::istringstream iss(strVal, std::istringstream::in);
			iss >> options.m_hsvd_comps;

			if( iss.fail() ) {
				std::cerr << "\nerror: couldn't recognise '" << strVal << "' as a number" << std::endl;
				return false;
			}
		}

		else if( strKey == "--max_hsvd_pts" ) {

			// convert string to number and check for errors
			std::istringstream iss(strVal, std::istringstream::in);
			iss >> options.m_max_hsvd_pts;

			if( iss.fail() ) {
				std::cerr << "\nerror: couldn't recognise '" << strVal << "' as a number" << std::endl;
				return false;
			}
		}
        
		// guess for the ppm value corresponding to 0Hz
		else if( strKey == "--ref" ) {

			// convert string to number and check for errors
			std::istringstream iss(strVal, std::istringstream::in);
			iss >> options.m_ref;
			options.m_ref_spec = true; 
			// this refers to the raw fid
			fid.SetPPMRef(options.m_ref);

			if( iss.fail() ) {
				std::cerr << "\nerror: couldn't recognise '" << strVal << "' as a number" << std::endl;
				return false;
			}
		}

		// auto-phasing
		else if( strKey == "--auto_phase" )
            options.m_bAutoPhase  = parse_binary(strVal);

		else if( strKey == "--combine_preproc" )
            options.m_bCombinePreproc  = parse_binary(strVal);

        else if( strKey == "--pre_fit_phase" )
            options.m_bPreFitPhase  = parse_binary(strVal);

        else if( strKey == "--pre_fit_shift" )
            options.m_bPreFitShift = parse_binary(strVal);

        else if( strKey == "--pre_fit_bl" )
            options.m_bPreFitBl = parse_binary(strVal);

        else if( strKey == "--append_lcm_basis" )
            options.m_bAppendLCMBasis = parse_binary(strVal);

        else if( strKey == "--append_neg_ref_basis" )
            options.m_bAppendNegRefBasis = parse_binary(strVal);
        
        else if( strKey == "--lipid_filter" )
            options.m_bLipidFilter = parse_binary(strVal);
        
        else if( strKey == "--lipid_filter_freq" ) 
		{

			// convert string to number and check for errors
			treal temp;
			std::istringstream iss(strVal, std::istringstream::in);
			iss >> temp;
			//
			if( iss.fail() ) {
				std::cerr << "\nerror: couldn't recognise '" << strVal << "' as a number" << std::endl;
				return false;
			}
			options.m_lipid_filt_freq = temp;
		}

        else if( strKey == "--swap_row_col" )
            options.m_bSwapRowCol = parse_binary(strVal);

        else if( strKey == "--full_echo" )
            options.m_bFullEcho = parse_binary(strVal);

		// auto-referencing
		else if( strKey == "--auto_ref" )
            options.m_bAutoRef = parse_binary(strVal);

		// max shift allowed
		else if( strKey == "--max_dref" ) 
		{

			// convert string to number and check for errors
			treal temp;
			std::istringstream iss(strVal, std::istringstream::in);
			iss >> temp;
			//
			if( iss.fail() ) {
				std::cerr << "\nerror: couldn't recognise '" << strVal << "' as a number" << std::endl;
				return false;
			}
			options.m_max_dref = temp;
		}

        // max metab shift allowed
		else if( strKey == "--max_metab_shift" ) 
		{

			// convert string to number and check for errors
			treal temp;
			std::istringstream iss(strVal, std::istringstream::in);
			iss >> temp;
			//
			if( iss.fail() ) {
				std::cerr << "\nerror: couldn't recognise '" << strVal << "' as a number" << std::endl;
				return false;
			}
			options.m_max_metab_shift = temp;
		}
        
        else if( strKey == "--nt_init_mu" ) 
		{

			// convert string to number and check for errors
			treal temp;
			std::istringstream iss(strVal, std::istringstream::in);
			iss >> temp;
			//
			if( iss.fail() ) {
				std::cerr << "\nerror: couldn't recognise '" << strVal << "' as a number" << std::endl;
				return false;
			}
			options.m_nt_init_mu = temp;
		}
        
        else if( strKey == "--nt_l2_dp" ) 
		{

			// convert string to number and check for errors
			treal temp;
			std::istringstream iss(strVal, std::istringstream::in);
			iss >> temp;
			//
			if( iss.fail() ) {
				std::cerr << "\nerror: couldn't recognise '" << strVal << "' as a number" << std::endl;
				return false;
			}
			options.m_nt_l2_dp = temp;
		}

        // max broad shift allowed
		else if( strKey == "--max_broad_shift" ) 
		{

			// convert string to number and check for errors
			treal temp;
			std::istringstream iss(strVal, std::istringstream::in);
			iss >> temp;
			//
			if( iss.fail() ) {
				std::cerr << "\nerror: couldn't recognise '" << strVal << "' as a number" << std::endl;
				return false;
			}
			options.m_max_broad_shift = temp;
		}

		else if( strKey == "--ref_file" )
		{
			options.m_strRefFile = strVal;
		}

		else if( strKey == "--av_list" )
		{
			options.m_strAvListFile = strVal;
		}

		else if( strKey == "--ref_freq" ) 
		{

			// convert string to number and check for errors
			treal temp;
			std::istringstream iss(strVal, std::istringstream::in);
			iss >> temp;
			//
			if( iss.fail() ) {
				std::cerr << "\nerror: couldn't recognise '" << strVal << "' as a number" << std::endl;
				return false;
			}
			options.m_ref_freq = temp;
		}

		// plot pausing 
		//else if( strKey == "--pause" )
        //    options.m_pause = parse_binary(strVal);

		// water-reference eddy current correction
		else if( strKey == "--water_eddy" )
            options.m_bWaterEddy = parse_binary(strVal);

		// LCModel mode option
		else if( strKey == "--basis_comp" )
            options.m_basis_comp = parse_binary(strVal);

		// LCModel mode option
		else if( strKey == "--ext_pdf" )
            options.m_bPdfExt = parse_binary(strVal);

        else if( strKey == "--stack_pdf" )
            options.m_bPdfStack = parse_binary(strVal);

        else if( strKey == "--ext_csv_fit" )
            options.m_ext_csv_fit = parse_binary(strVal);

		// fast fit?
		//else if( strKey == "--ff" )
         //   options.m_ff = parse_binary(strVal);

		// show preprocessed result
		else if( strKey == "--show_pre" )
            options.m_bShowPreprocessed = parse_binary(strVal);

		// only fit SVS data?
		else if( strKey == "--svs_only" )
            options.m_svs_only = parse_binary(strVal);

		// display an old results file
		else if( strKey == "--view_fit" ) {
			options.m_strViewFile = strVal;
		}

		// display FID paras
		else if( strKey == "--print_paras" )
            options.m_bPrintParas = parse_binary(strVal);

		else if( strKey == "--no_fit" )
            options.m_bNofit = parse_binary(strVal);

		else if( strKey == "--read_only" )
            options.m_bReadOnly = parse_binary(strVal);

		else if( strKey == "--rw_only" )
            options.m_bReadWriteOnly = parse_binary(strVal);

		// comma seperated list of signals to be plotted (no spaces allowed!)
		else if( strKey == "--plot_sigs" ) {
			options.m_strPlotSigs = strVal;
		}
        
        // title of pdf results
		else if( strKey == "--title" ) {
			options.m_title = strVal;
		}

        // pdf fit line col
		else if( strKey == "--fit_col" ) {
			options.m_fit_color = strVal;
		}

        // lambda
		else if( strKey == "--lambda" ) {
            std::istringstream iss(strVal, std::istringstream::in);
			iss >> options.m_lambda;
		}

		// where the basis XML will go
		else if( strKey == "--output_basis" ) {
			options.m_strBasisSaveFile = strVal;
		}
        
		// where the basis csv will go
		else if( strKey == "--output_basis_csv" ) {
			options.m_strBasisSaveFileCSV = strVal;
		}
        
        // where the basis lcm will go
		else if( strKey == "--output_basis_lcm" ) {
			options.m_strBasisSaveFileLCM = strVal;
		}

		// starting point of fid
		else if( strKey == "--lcm_basis_pts" ) {
			// convert string to number and check for errors
			std::istringstream iss(strVal, std::istringstream::in);
			iss >> options.m_lcm_basis_pts;

			if( iss.fail() ) {
				std::cerr << "\nerror: couldn't recognise '" << strVal << "' as a number" << std::endl;
				return false;
			}
		}

		// write an intermediate dangerplot file of the raw FID
		else if( strKey == "--write_raw" ) {
			options.m_strOutRaw = strVal;
		}

		// write an intermediate dangerplot file of the raw FID
		else if( strKey == "--write_raw_v3" ) {
			options.m_strOutRawV3 = strVal;
		}
        
        // write an intermediate dangerplot file of the raw water ref FID
		else if( strKey == "--write_raw_w" ) {
			options.m_strOutRawW = strVal;
		}

        // write an intermediate dangerplot file of the raw water ref FID
		else if( strKey == "--write_raw_w_v3" ) {
			options.m_strOutRawWV3 = strVal;
		}
        
        // write an intermediate dangerplot file of the raw FID
		else if( strKey == "--write_raw_lcm" ) {
			options.m_strOutRawLcm = strVal;
		}
        
        // write an intermediate dangerplot file of the raw water ref FID
		else if( strKey == "--write_raw_lcm_w" ) {
			options.m_strOutRawLcmW = strVal;
		}

		// write an intermediate dangerplot file of the preprocess FID
		else if( strKey == "--write_pre" ) {
			options.m_strOutPre = strVal;
		}

		// write an intermediate dangerplot file of the postprocessed FID
		else if( strKey == "--write_post" ) {
			options.m_strOutPost = strVal;
		}

		// this is a dirty hack to for GE to work
		/*else if( strKey == "--ge_offset" ) {

			std::istringstream iss(strVal, std::istringstream::in);
			iss >> options.m_geOptions.nOffset;
		}*/

		/*else if( strKey == "--ge_wsframes" ) {

			std::istringstream iss(strVal, std::istringstream::in);
			iss >> options.m_geOptions.nWSFrames;
		}*/

		/*else if( strKey == "--ge_wframes" ) {

			std::istringstream iss(strVal, std::istringstream::in);
			iss >> options.m_geOptions.nWaterFrames;
		}*/

		/*else if( strKey == "--ge_coils" ) {

			std::istringstream iss(strVal, std::istringstream::in);
			iss >> options.m_geOptions.nCoils;
		}*/

		/*else if( strKey == "--ge_samples" ) {

			std::istringstream iss(strVal, std::istringstream::in);
			iss >> options.m_geOptions.nFieldSize;
		}*/

		else if( strKey == "--para_file" ) {
            // skip over this one
        }

		else if( strKey == "--" ) {
            // skip over this one
        }
        
        else
        { 
            // if we got here there was some kind of problem
            if (!handler_key_found)
            {
                std::cerr << "\nerror: couldn't recognise '" << strKey << "' as a valid option" << std::endl;
                return false;
            }
        }


	}

	if( options.m_strFile == "" && options.m_strViewFile.size() == 0 && options.m_bPrintParas == false && std::string(argv[0]) != "tarquingui" ) 
	{
		std::cerr << "\nerror: you must specify at least:";
		std::cerr << "\n    --input";
		std::cerr << "\n" << std::endl;
		std::cerr << "\nuse --help to see a list of common options" << std::endl;

		return false;
	}
    
    /*
	// if the user didn't specify the output file, default to something sensible
	if( 0 == options.GetOutputXMLPath().size() && options.m_strViewFile.size() == 0 ) {

		//std::cout << "\nOutput file not specified, defaulting to \"TARQUINresults.xml\"." << std::endl;

		options.m_strOutputXMLPath = "results.xml";
		options.m_strFileOutCSV = "results.csv";
	}
    */

	// if the user didn't specify the output file, default to something sensible
	//if( 0 == options.GetFilenameTxt().size() && options.m_strViewFile.size() == 0 ) 
	    //options.m_strFileOutTxt = "results.txt";

	return true;
}

