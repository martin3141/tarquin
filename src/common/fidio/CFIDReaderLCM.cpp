#include "CFIDReaderLCM.hpp"
#include "CBoswell.hpp"
#include "CFID.hpp"

#include <string>
#include <sstream>

tarquin::CFIDReaderLCM::CFIDReaderLCM(CFID& fid, CBoswell& log) : 
		CFIDReader(fid, log)
{
}

void tarquin::CFIDReaderLCM::Load(std::string strFilename, const Options& opts, CBoswell& log)
{
	// use the boost tokenizer that keeps things like strings in quotes intact
	// middle argument is token separator 
	boost::escaped_list_separator<char> sep("\\", "= ", "'");

	Lex(strFilename, sep);

	EatTokens();
}

bool IsNumber(const std::string& str)
{
	std::istringstream strm(str);

	tarquin::treal x;

	strm >> x;

	return not strm.fail();
}

void tarquin::CFIDReaderLCM::EatTokens()
{
	// the FID will be multiplied by this number if we find it
	treal tramp = 1.0;

	for( TokenList::iterator it = m_tokens.begin(); it != m_tokens.end(); ++it ) 
	{
		std::string strKey = it->first;

		// skip block markers
		if( strKey == "$SEQPAR" ) 
			continue;

		if( strKey == "$NMID") 
			continue;

		std::string strValue = (it+1)->first;

		// magic conversion
		std::istringstream strmValue;
		strmValue.clear();
		strmValue.str(strValue);

		if( strKey == "echot" ) 
		{
			treal tau;
			strmValue >> tau;

			m_fid.SetEchoTime(tau/1000);				
		}

		else if( strKey == "tramp" ) 
		{
			strmValue >> tramp;
		}

		// there may be many of these, so check if the next token is a number
		else if( strKey == "$END" ) 
		{
			if( true == IsNumber(strValue) ) 
			{
				// grab all (real,imag) time domain samples
				cvm::cvector y;
				double real_part, imag_part;	

				std::istringstream ins;

				for( TokenList::iterator in = it+1; in != m_tokens.end(); ) 
				{
					ins.clear();
					ins.str(in->first);
					ins >> real_part;

					if( ins.fail() ) 
						throw Exception("error reading FID tokens: '%s' should be a number.", in->first.c_str());

					++in;

					ins.clear();
					ins.str(in->first);
					ins >> imag_part;

					if( ins.fail() ) 
						throw Exception("error reading FID tokens: '%s' should be a number.", in->first.c_str());

					++in;

					y.push_back(tcomplex(real_part, imag_part));
				}

				y /= tramp;

				m_fid.AppendFromVector(y);
				return;
			}
		}

		// if we got here then we did indeed find a (key,value) pair
		++it;
	}

	throw Exception("got to end of file without finding FID data");
}

