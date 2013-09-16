#include "CFIDReaderDPT.hpp"
#include "CBoswell.hpp"
#include "CFID.hpp"

#include <string>
#include <sstream>

tarquin::CFIDReaderDPT::CFIDReaderDPT(CFID& fid, CBoswell& log) : 
	CFIDReader(fid, log)
{
}

void tarquin::CFIDReaderDPT::Load(std::string strFilename, const Options& opts, CBoswell& log)
{
	// use the boost tokenizer that keeps things like strings in quotes intact
	// middle argument is token separator 
	boost::escaped_list_separator<char> sep("\\", " \t\r", "\"");

	Lex(strFilename, sep); 

	EatTokens();
}

void tarquin::CFIDReaderDPT::EatTokens()
{
	std::istringstream ins;

	// this field is not stored explicitly, but is stored as size property of FID size
	size_t nSamples = 0;

	for( TokenList::iterator it = m_tokens.begin(); it != m_tokens.end(); it++ ) {

		// clear up for this conversion
		ins.clear();

		if( it->first == "Number_of_points" ) 
		{
			it++;
			ins.str(it->first);
			ins >> nSamples;
            if ( nSamples > 16384 )
                nSamples = 16384;
		}
		else if( it->first == "Sampling_frequency" ) 
		{
			it++;
			ins.str(it->first);
			treal fs = 0;
			ins >> fs;
			m_fid.SetSamplingFrequency(fs);
		}
		else if( it->first == "Transmitter_frequency" ) 
		{
			it++;
			ins.str(it->first);
			treal ft = 0;
			ins >> ft;
			m_fid.SetTransmitterFrequency(ft);
		}
		else if( it->first == "Phi0" ) 
		{
			it++;
			ins.str(it->first);
			treal phi0 = 0;
			ins >> phi0;
            pair_vec phi0_vec;
            phi0_vec.push_back(std::make_pair(phi0, true));
			m_fid.SetPhi0(phi0_vec);
		}
		else if( it->first == "Phi1" ) 
		{
			it++;
			ins.str(it->first);
			treal phi1 = 0;
			ins >> phi1;
            pair_vec phi1_vec;
            phi1_vec.push_back(std::make_pair(phi1, true));
			m_fid.SetPhi1(phi1_vec);
		}
		else if( it->first == "PPM_reference" ) 
		{
			it++;
			ins.str(it->first);
			treal ref = 0;
			ins >> ref;
			pair_vec ref_vec;
			ref_vec.push_back(std::make_pair(ref, true));
			m_fid.SetPPMRef(ref_vec);
		}
		else if( it->first == "Pulse_sequence") 
		{
			it++;
			ins.str(it->first);
			std::string sequence;
			ins >> sequence;
			m_fid.SetPulseSequence(sequence);
		}
		else if( it->first == "Averages" ) 
		{
			it++;
			ins.str(it->first);
			int averages = 0;
			ins >> averages;
			m_fid.SetAverages(averages);
		}
		else if( it->first == "Echo_time") 
		{
			it++;
			ins.str(it->first);
			treal tau = 0;
			ins >> tau;
			m_fid.SetEchoTime(tau);
		}
		else if( it->first == "Real_FID" ) 
		{
			it+=2; // loose "imaginary part" header
			return EatTokensFID(it, nSamples);
		}
		else 
		{
			// skip all tokens on this line
			size_t nLineCurrent = it->second;

			while( it->second == nLineCurrent && it != m_tokens.end() )
				it++;

			// skip one back because of outer increment
			it--;
		}
	}

	// if we haven't read any samples things didn't work well
	if( 0 == nSamples )
		throw Exception("failed to ready any samples from the file");
}

void tarquin::CFIDReaderDPT::EatTokensFID(TokenList::iterator it_token, std::size_t nSamples)
{
	assert( nSamples > 0 );

	std::istringstream ins;

	// allocate space for signal	
	cvm::cvector FID(nSamples);

	// read in points
	for(std::size_t n = 0; n < nSamples; n++) 
	{
		double real_part, imag_part;	

		ins.clear();
		ins.str(it_token->first);
		ins >> real_part;

		if( ins.fail() ) 
		{
			throw Exception("error reading FID tokens: '%s' should be a number.", it_token->first.c_str());
		}

		it_token++;

		ins.clear();
		ins.str(it_token->first);
		ins >> imag_part;

		if( ins.fail() ) 
		{
			throw Exception("error reading FID tokens: '%s' should be a number.", it_token->first.c_str());
		}

		it_token++;

		FID[n+1] = tcomplex(real_part, imag_part);
	}

    // truncate FID if long, prob need a warning here...
    if ( FID.size() > 16384 )
        FID.resize(16384);

	m_fid.AppendFromVector(FID);
}

