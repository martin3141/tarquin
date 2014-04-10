#include "console.hpp"
#include "version/version.h"

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
	std::cout << "\n\t--ge_offset         offset_to_data";
	std::cout << "\n\t--ge_wsframes       number of ws averages";
	std::cout << "\n\t--ge_wframes        number of w averages";
	std::cout << "\n\t--ge_coils          number of coils";
	std::cout << "\n\t--ge_samples        number of samples";
	std::cout << "\n\t--start_pnt         starting sample (1-based)";
	std::cout << "\n\t--end_pnt           ending sample (1-based)";
	std::cout << "\n\t--max_iters         maximum number of iterations to perform";
	std::cout << "\n\t--water_width       width of hsvd in hz";
	std::cout << "\n\t--conv_width        width of convolution window in points";
	std::cout << "\n\t--lipid_filter      {true | false} remove signals upfield of 1.8ppm in hsvd";
	std::cout << "\n\t--swap_row_col      {true | false} swap CSI rows and cols";
	std::cout << "\n\t--auto_phase        {true | false}";
	std::cout << "\n\t--auto_ref          {true | false}";
	std::cout << "\n\t--max_dref          the max deviation from ref allowed by auto_ref";
	std::cout << "\n\t--max_phi1          the value of phi1_max/fs/2";
	std::cout << "\n\t--ref_file          CSV file containing reference peak list";
	std::cout << "\n\t--ref_signals       {1h_naa_cr_cho_lip | 1h_naa_cho | 1h_naa_cr_cho | 1h_cr_cho | 1h_naa | 1h_cr | 1h_cho | 1h_h2o | 1h_lip | 31p_pcr | 31p_pcr_gammaapt}";
	std::cout << "\n\t--dref_signals      reference signals for dynamic frequency correction, options as above";
	std::cout << "\n\t--fs                sampling frequency in Hz";
	std::cout << "\n\t--ft                transmitter frequency in Hz";
	std::cout << "\n\t--basis_csv         path to basis (CSV files)";
	std::cout << "\n\t--basis_xml         path to basis (precompiled XML files)";
	std::cout << "\n\t--basis_lcm         path to basis (LCModel .basis format)";
	std::cout << "\n\t--int_basis         {1h_brain | 1h_brain_gly_glth | 1h_brain_gly_cit_glth | 1h_brain_full | 1h_brain_le | megapress_gaba | 31p_brain}";
	std::cout << "\n\t--echo              echo time in seconds";
	std::cout << "\n\t--te1               te1 time in seconds for PRESS sequence";
	std::cout << "\n\t--tm                tm time in seconds for STEAM sequence";
	std::cout << "\n\t--cpmg_pulses       number of pulses for CPMG sequence";
	std::cout << "\n\t--pul_seq           {press | steam | slaser | laser | pulse_acq | cpmg | se | mega_press}";
	std::cout << "\n\t--dyn_av            {default | none | all | subtract | odd | even} water sup. dyn averaging scheme";
	std::cout << "\n\t--dyn_av_w          {default | none | all | subtract | odd | even} water dyn averaging scheme";
	std::cout << "\n\t--dyn_freq_corr     {true | false}";
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
	std::cout << "\n\t--output_basis      basis XML file";
	std::cout << "\n\t--output_basis_lcm  LCModel basis file";
	std::cout << "\n\t--output_image      plot of the fit in pdf format";
	std::cout << "\n\t--output_pdf        A4 pdf results page for printing";
	std::cout << "\n\t--ext_pdf           Extended output {true | false}";
	std::cout << "\n\t--title             title of results page";
	std::cout << "\n\t--output_txt        txt output of results";
	std::cout << "\n\t--output_csv        csv output of results";
	std::cout << "\n\t--output_fit        csv output of fit";
	std::cout << "\n\t--output_spec       csv output of processed spectra";
	std::cout << "\n\t--output_spec_m     csv output of processed spectra";
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

bool tarquin::ParseCommandLine(int argc, char* argv[], Options& options, CFID& fid)
{

	// create an empty vector of strings
	std::vector <std::string> args;
	// convert array of pointers to a vector of strings
	for(int n = 0; n < argc; n++ ) {
		args.push_back(argv[n]);
	}

	// search for a para_file before processing options
	for(size_t n = 1; n < args.size(); n+=2 ) {

		std::string strKey = std::string(args[n]);

		if( strKey == "--help" ) {
			DisplayUsage();
			return false;
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
	for(size_t n = 1; n < args.size(); n+=2 ) {

		std::string strKey = std::string(args[n]);

		if( n == args.size()-1 ) {
			std::cerr << "\nerror: expected value after key '" << strKey << "'" << std::endl;
			return false;
		}

		std::string strVal = std::string(args[n+1]);

        //std::cout << strKey << std::endl;
        //std::cout << strVal << std::endl;

		// name of FID file
		if( strKey == "--input" ) 
        {
			options.m_strFile = strVal;
        }

		// name of water reference fid file
		else if( strKey == "--input_w" )
			options.m_strFileWater = strVal;

		// format of FID object
		else if( strKey == "--format" ) {

			if( strVal == "siemens" ) 
				options.m_format = tarquin::SIEMENS;

			else if( strVal == "philips" ) 
				options.m_format = tarquin::PHILIPS;

			else if( strVal == "ge" ) 
				options.m_format = tarquin::GE;

			else if( strVal == "dpt" ) 
				options.m_format = tarquin::DANGER;

			else if( strVal == "rda" )
				options.m_format = tarquin::RDA;

			else if( strVal == "lcm" )
				options.m_format = tarquin::LCM;

			else if( strVal == "shf" )
				options.m_format = tarquin::SHF;

			else if( strVal == "varian" )
				options.m_format = tarquin::VARIAN;

			else if( strVal == "bruker" )
				options.m_format = tarquin::BRUKER;
            
            else if( strVal == "philips_dcm" )
				options.m_format = tarquin::PHILIPS_DCM;

            else if( strVal == "dcm" )
				options.m_format = tarquin::DCM;

            else if( strVal == "jmrui_txt" )
				options.m_format = tarquin::JMRUI_TXT;
			else {
				std::cerr << "\nerror: unrecognised format '" << strVal << "'" << std::endl;
				return false;
			}
		}

		// format of FID object
		else if( strKey == "--ref_signals" ) {
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

			else {
				std::cerr << "\nerror: unrecognised ref mode '" << strVal << "'" << std::endl;
				return false;
			}
        }

		else if( strKey == "--dref_signals" ) {
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
		else if( strKey == "--pul_seq" ) {
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
		else if( strKey == "--dyn_av" ) {
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

		else if( strKey == "--dyn_av_w" ) {
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
        
        // pre dyn av freq correction
        else if( strKey == "--dyn_freq_corr" ) {
			if( strVal == "true" )
				options.m_dyn_freq_corr = true;
			else
				options.m_dyn_freq_corr = false;
		}

        
        // interal basis
		// default to none
		else if( strKey == "--int_basis" ) {
			if( strVal == "1h_brain" ) 
				options.SetIntBasisSet(PROTON_BRAIN);
			else if( strVal == "1h_brain_gly_glth" ) 
				options.SetIntBasisSet(PROTON_BRAIN_GLY_GLTH);
			else if( strVal == "1h_brain_gly_cit_glth" ) 
				options.SetIntBasisSet(PROTON_BRAIN_GLY_CIT_GLTH);
			else if( strVal == "1h_brain_full" ) 
				options.SetIntBasisSet(PROTON_BRAIN_FULL);
			else if( strVal == "1h_brain_le" ) 
				options.SetIntBasisSet(PROTON_BRAIN_LE);
			else if( strVal == "megapress_gaba" ) 
				options.SetIntBasisSet(PROTON_MEGAPRESS_GABA);
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

		else if( strKey == "--ppm_right" ) {
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
		}

		else if( strKey == "--lb" ) {
			// convert string to number and check for errors
			treal temp;
			std::istringstream iss(strVal, std::istringstream::in);
			iss >> temp;
			//
			if( iss.fail() ) {
				std::cerr << "\nerror: couldn't recognise '" << strVal << "' as a number" << std::endl;
				return false;
			}
			options.m_lb = temp;
		}

		else if( strKey == "--init_beta" ) {
			// convert string to number and check for errors
			treal temp;
			std::istringstream iss(strVal, std::istringstream::in);
			iss >> temp;
			//
			if( iss.fail() ) {
				std::cerr << "\nerror: couldn't recognise '" << strVal << "' as a number" << std::endl;
				return false;
			}
			options.m_init_beta = temp;
		}
		
        else if( strKey == "--max_beta" ) {
			// convert string to number and check for errors
			treal temp;
			std::istringstream iss(strVal, std::istringstream::in);
			iss >> temp;
			//
			if( iss.fail() ) {
				std::cerr << "\nerror: couldn't recognise '" << strVal << "' as a number" << std::endl;
				return false;
			}
			options.m_max_beta = temp;
		}

        else if( strKey == "--beta_scale" ) {
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
		}

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
		}

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
		else if( strKey == "--adapt_sp" ) {
			if( strVal == "true" )
				options.m_bAdaptSp = true;
			else
				options.m_bAdaptSp = false;
		}

		// pre-hsvd
		else if( strKey == "--pre_hsvd" ) {
			if( strVal == "true" )
				options.m_pre_hsvd = true;
			else
				options.m_pre_hsvd = false;
		}

		// old phi1 limits
		else if( strKey == "--old_phase" ) {
			if( strVal == "true" )
				options.m_old_phase = true;
			else
				options.m_old_phase = false;
		}

		// adaptive ending point of fid
		else if( strKey == "--adapt_ep" ) {
			if( strVal == "true" )
				options.m_bAdaptEp = true;
			else
				options.m_bAdaptEp = false;
		}
        
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
		else if( strKey == "--output_fit" ) {

			options.m_strFileOutCSVFit = strVal;
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
        {
			if( strVal == "true" )
				options.m_au_norm = true;
			else
				options.m_au_norm = false;
		}

		// water removal window width
		else if( strKey == "--water_width" ) {

			// convert string to number and check for errors
			std::istringstream iss(strVal, std::istringstream::in);
			iss >> options.m_water_window;

			if( iss.fail() ) {
				std::cerr << "\nerror: couldn't recognise '" << strVal << "' as a number" << std::endl;
				return false;
			}
		}

		// water removal convolution window width
		else if( strKey == "--conv_width" ) {

			// convert string to number and check for errors
			std::istringstream iss(strVal, std::istringstream::in);
			iss >> options.m_conv_window_width;

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
		else if( strKey == "--auto_phase" ) {

			if( strVal == "true" )
				options.m_bAutoPhase = true;
			else
				options.m_bAutoPhase = false;
		}

		else if( strKey == "--combine_preproc" ) {

			if( strVal == "true" )
				options.m_bCombinePreproc = true;
			else
				options.m_bCombinePreproc = false;
		}
        
        else if( strKey == "--lipid_filter" ) {

			if( strVal == "true" )
				options.m_bLipidFilter = true;
			else
				options.m_bLipidFilter = false;
		}

        else if( strKey == "--swap_row_col" ) {

			if( strVal == "true" )
				options.m_bSwapRowCol = true;
			else
				options.m_bSwapRowCol = false;
		}

        else if( strKey == "--full_echo" ) {

			if( strVal == "true" )
				options.m_bFullEcho = true;
			else
				options.m_bFullEcho = false;
		}

		// auto-referencing
		else if( strKey == "--auto_ref" ) {

			if( strVal == "true" )
				options.m_bAutoRef = true;
			else
				options.m_bAutoRef = false;
		}

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

		// plot pausing 
		else if( strKey == "--pause" ) {

			if( strVal == "true" )
				options.m_pause = true;
			else
				options.m_pause = false;
		}

		// water-reference eddy current correction
		else if( strKey == "--water_eddy" ) {

			if( strVal == "true" )
				options.m_bWaterEddy = true;
			else
				options.m_bWaterEddy = false;
		}

		// LCModel mode option
		else if( strKey == "--basis_comp" ) {

			if( strVal == "true" )
				options.m_basis_comp = true;
			else
				options.m_basis_comp = false;
		}

		// LCModel mode option
		else if( strKey == "--ext_pdf" ) {

			if( strVal == "true" )
				options.m_bPdfExt = true;
			else
				options.m_bPdfExt = false;
		}

		// fast fit?
		else if( strKey == "--ff" ) {

			if( strVal == "true" )
				options.m_ff = true;
			else
				options.m_ff = false;
		}

		// show preprocessed result
		else if( strKey == "--show_pre" ) {

			if( strVal == "true" )
				options.m_bShowPreprocessed = true;
			else
				options.m_bShowPreprocessed = false;
		}

		// only fit SVS data?
		else if( strKey == "--svs_only" ) {

			if( strVal == "true" )
				options.m_svs_only = true;
			else
				options.m_svs_only = false;
		}

		// display an old results file
		else if( strKey == "--view_fit" ) {
			options.m_strViewFile = strVal;
		}

		// display FID paras
		else if( strKey == "--print_paras" ) {
			if( strVal == "true" )
				options.m_bPrintParas = true;
			else
				options.m_bPrintParas = false;
		}

		else if( strKey == "--no_fit" ) {
			if( strVal == "true" )
				options.m_bNofit = true;
			else
				options.m_bNofit = false;
		}

		// comma seperated list of signals to be plotted (no spaces allowed!)
		else if( strKey == "--plot_sigs" ) {
			options.m_strPlotSigs = strVal;
		}
        
        // title of pdf results
		else if( strKey == "--title" ) {
			options.m_title = strVal;
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
        
        // write an intermediate dangerplot file of the raw water ref FID
		else if( strKey == "--write_raw_w" ) {
			options.m_strOutRawW = strVal;
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
		else if( strKey == "--ge_offset" ) {

			std::istringstream iss(strVal, std::istringstream::in);
			iss >> options.m_geOptions.nOffset;
		}

		else if( strKey == "--ge_wsframes" ) {

			std::istringstream iss(strVal, std::istringstream::in);
			iss >> options.m_geOptions.nWSFrames;
		}

		else if( strKey == "--ge_wframes" ) {

			std::istringstream iss(strVal, std::istringstream::in);
			iss >> options.m_geOptions.nWaterFrames;
		}

		else if( strKey == "--ge_coils" ) {

			std::istringstream iss(strVal, std::istringstream::in);
			iss >> options.m_geOptions.nCoils;
		}

		else if( strKey == "--ge_samples" ) {

			std::istringstream iss(strVal, std::istringstream::in);
			iss >> options.m_geOptions.nFieldSize;
		}
		else if( strKey == "--para_file" ) {
            // skip over this one
        }

		else if( strKey == "--" ) {
            // skip over this one
        }
        
        else
        { 
            // if we got here there was some kind of problem
            std::cerr << "\nerror: couldn't recognise '" << strKey << "' as a valid option" << std::endl;
            return false;
        }
	}

	if( options.m_strFile == "" && options.m_strViewFile.size() == 0 && options.m_bPrintParas == false && std::string(argv[0]) != "tarquingui" ) 
	{
		std::cerr << "\nerror: you must specify at least:";
		std::cerr << "\n    --input";
		std::cerr << "\nOR:";
		std::cerr << "\n    --view_fit" << std::endl;
		std::cerr << "\nOR:";
		std::cerr << "\n    --input" << std::endl;
		std::cerr << "\n    --print_paras" << std::endl;

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
	if( 0 == options.GetFilenameTxt().size() && options.m_strViewFile.size() == 0 ) 
		options.m_strFileOutTxt = "results.txt";

	return true;
}

