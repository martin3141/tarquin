#include "CFIDReaderJMRUITXT.hpp"
#include "CBoswell.hpp"
#include "CFID.hpp"

#include <string>
#include <sstream>

tarquin::CFIDReaderJMRUITXT::CFIDReaderJMRUITXT(CFID& fid, CBoswell& log) : 
	CFIDReader(fid, log)
{
}

void tarquin::CFIDReaderJMRUITXT::Load(std::string strFilename, const Options& opts, CBoswell& log)
{
	// use the boost tokenizer that keeps things like strings in quotes intact
	// middle argument is token separator 
	boost::escaped_list_separator<char> sep("\\", " \t\r", "\"");

	Lex(strFilename, sep); 

	EatTokens(opts);
}

void tarquin::CFIDReaderJMRUITXT::EatTokens(const Options& opts)
{
	std::istringstream ins;

	// this field is not stored explicitly, but is stored as size property of FID size
	size_t nSamples = 0;

	for( TokenList::iterator it = m_tokens.begin(); it != m_tokens.end(); it++ ) {

		// clear up for this conversion
		ins.clear();

		if( it->first == "PointsInDataset:" ) 
		{
			it++;
			ins.str(it->first);
			ins >> nSamples;
            if ( nSamples > 16384 )
                nSamples = 16384;
		}
		else if( it->first == "DatasetsInFile:" ) 
		{
			it++;
			ins.str(it->first);
			size_t cols = 0;
			ins >> cols;
			m_fid.SetCols(cols);
		}
		else if( it->first == "SamplingInterval:" ) 
		{
			it++;
			ins.str(it->first);
			treal fs = 0;
			ins >> fs;
			m_fid.SetSamplingFrequency(1/(fs*1e-3));
		}
		else if( it->first == "TransmitterFrequency:" ) 
		{
			it++;
			ins.str(it->first);
			treal ft = 0;
			ins >> ft;
			m_fid.SetTransmitterFrequency(ft);
		}
		else if( it->first == "Signal" ) 
		{
			it+=14;
			return EatTokensFID(it, nSamples, opts, m_tokens.end());
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

void tarquin::CFIDReaderJMRUITXT::EatTokensFID(TokenList::iterator it_token, std::size_t nSamples, const Options& opts, TokenList::iterator end_token)
{
	assert( nSamples > 0 );
    
    /*
    if ( opts.GetDynAv() == ALL || opts.GetDynAv() == ODD || opts.GetDynAv() == EVEN )
    {
        std::istringstream ins;

        // allocate space for signal	
        cvm::cvector FID(nSamples);

        for ( int m = 0; m < m_fid.GetVoxelCount(); m++ )
        {
            for(std::size_t n = 0; n < nSamples; n++) 
            {
                double real_part, imag_part;	
                ins.clear();
                if ( it_token > end_token )
                    throw Exception("error file ended unexpectedly");

                ins.str(it_token->first);

                ins >> real_part;

                if( ins.fail() ) 
                {
                    throw Exception("error reading FID tokens: '%s' should be a number.", it_token->first.c_str());
                }

                it_token++;

                ins.clear();
                if ( it_token > end_token )
                    throw Exception("error file ended unexpectedly");

                ins.str(it_token->first);
                ins >> imag_part;

                if( ins.fail() ) 
                {
                    throw Exception("error reading FID tokens: '%s' should be a number.", it_token->first.c_str());
                }

                it_token+=3;
                
                if ( ( m % 2 == 0 ) && ( opts.GetDynAv() == ODD ) )
                    FID[n+1] += tcomplex(real_part, -imag_part);
                else if ( ( m % 2 == 1 ) && ( opts.GetDynAv() == EVEN ) )
                    FID[n+1] += tcomplex(real_part, -imag_part);
                else if ( opts.GetDynAv() == ALL )
                    FID[n+1] += tcomplex(real_part, -imag_part);

            }
            it_token+=7;
        }

        // truncate FID if long, prob need a warning here...
        if ( FID.size() > 16384 )
            FID.resize(16384);

        // vector is cleared because mega-press data gets loaded twice
        // when using the gui
        m_fid.ClearVector();

        // only one col now since data has been collapsed
	    m_fid.SetCols(1);

        m_fid.AppendFromVector(FID);
    }
    */
    //else if ( opts.GetDynAv() == NONE )
    //{
    std::istringstream ins;

    // allocate space for signal	
    cvm::cvector FID(nSamples);

    for ( int m = 0; m < m_fid.GetVoxelCount(); m++ )
    {
        // read in points
        for(std::size_t n = 0; n < nSamples; n++) 
        {
            double real_part, imag_part;	
            ins.clear();
            if ( it_token > end_token )
                throw Exception("error file ended unexpectedly");

            ins.str(it_token->first);

            ins >> real_part;

            if( ins.fail() ) 
            {
                throw Exception("error reading FID tokens: '%s' should be a number.", it_token->first.c_str());
            }

            it_token++;

            ins.clear();

            if ( it_token > end_token )
                throw Exception("error file ended unexpectedly");

            ins.str(it_token->first);
            ins >> imag_part;

            if( ins.fail() ) 
            {
                throw Exception("error reading FID tokens: '%s' should be a number.", it_token->first.c_str());
            }

            it_token+=3;

            FID[n+1] = tcomplex(real_part, -imag_part);
        }

        // truncate FID if long, prob need a warning here...
        if ( FID.size() > 16384 )
            FID.resize(16384);

        m_fid.AppendFromVector(FID);

        it_token+=7;
    }
    //}
    /*else
        throw Exception("Dynamic averaging scheme incompatible with this data format");*/
}

