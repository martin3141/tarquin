#include "CFIDReaderBruker.hpp"
#include "CBoswell.hpp"
#include "CFID.hpp"
#include "common.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <limits>
#include <string>
#include <sstream>
#include <bitset>
#include <boost/filesystem.hpp>

tarquin::CFIDReaderBruker::CFIDReaderBruker(CFID& fid, CBoswell& log) : 
	CFIDReader(fid, log)
{
}

void tarquin::CFIDReaderBruker::Load(std::string strFilename, const Options& opts, CBoswell& log)
{
	// strip off the filename and replace with acqus
	boost::filesystem::path fid_path;
	fid_path = strFilename;

    boost::filesystem::path fid_dir = fid_path.parent_path();
	
    boost::filesystem::path acqus_path = fid_dir;
    acqus_path /= "acqus";

	std::string str_acqus = acqus_path.generic_string();

	// open the procpar text file for parsing
	std::string line;
	std::ifstream acqus ( str_acqus.c_str() );

	// number of data points
	int N;
	int decim;
	int dspfvs;
	int bytorda;
	int grpdly;
	int digmod;
	if ( acqus.is_open() )
	{
		while (! acqus.eof() )
		{
			getline (acqus,line);
	        size_t len = line.size();
			if ( line.substr(0,8) == "##$SW_h=" )
			{
				float fs;
				from_string<float>(fs, line.substr(8,len), std::dec);
				m_fid.SetSamplingFrequency(fs);
				//std::cout << fs << std::endl;
			} 
			else if ( line.substr(0,8) == "##$SFO1=" )
			{
				float ft;
				from_string<float>(ft, line.substr(8,len), std::dec);
				m_fid.SetTransmitterFrequency(ft*1e6);
				//std::cout << ft << std::endl;
			}
			else if ( line.substr(0,6) == "##$TD=" )
			{
				from_string<int>(N, line.substr(6,len), std::dec);
				//std::cout << N << std::endl;
			}
			else if ( line.substr(0,9) == "##$DECIM=" )
			{
				from_string<int>(decim, line.substr(9,len), std::dec);
				//std::cout << decim << std::endl;
			}
			else if ( line.substr(0,10) == "##$DSPFVS=" )
			{
				from_string<int>(dspfvs, line.substr(10,len), std::dec);
				//std::cout << dspfvs << std::endl;
			}
			else if ( line.substr(0,11) == "##$BYTORDA=" )
			{
				from_string<int>(bytorda, line.substr(11,len), std::dec);
				//std::cout << bytorda << std::endl;
			}
			else if ( line.substr(0,10) == "##$GRPDLY=" )
			{
				from_string<int>(grpdly, line.substr(10,len), std::dec);
				//std::cout << grpdly << std::endl;
			}
			else if ( line.substr(0,10) == "##$DIGMOD=" )
			{
				from_string<int>(digmod, line.substr(10,len), std::dec);
				//std::cout << digmod << std::endl;
			}
		}
		acqus.close();
	}
	else 
	{
		throw Exception("Unable to open acqus file"); 
	}

	// number of digital fiter points to consider 
	int skip_pts;
	if ( digmod == 0 )
		skip_pts = 0; 
	else if ( 20 <= dspfvs && dspfvs <= 23 )
		skip_pts = grpdly;
	else
		skip_pts = (int) lookup_skip_points(decim, dspfvs);


	//std::cout << skip_pts << std::endl;

	// read file as int 32
	std::ifstream file(strFilename.c_str(), std::ios::binary);

	if( true == file.fail() )
		throw Exception("failed when reading FID segment of file");

	std::vector<tcomplex> samples;

	while( false == file.eof() ) 
	{
		int real_part;
		file.read((char *) &real_part, sizeof real_part);
		int imag_part;
		file.read((char *) &imag_part, sizeof imag_part);

		// eof is not true until we have actually tried to read past it
		if( true == file.eof() )
			break;

		// byteswap where appropriate
		if ( bytorda == 1 )
		{
			real_part = LongSwap(real_part);
			imag_part = LongSwap(imag_part);	
		}

		tcomplex z(real_part, -imag_part);

		// store this sample
		samples.push_back(z);
	}

	/*std::cout << samples.size() << std::endl;
	  std::cout << N << std::endl;
	  std::cout << samples[0] << std::endl;
	  std::cout << samples[1] << std::endl;
	  std::cout << samples[2] << std::endl;
	  */

	std::vector<tcomplex> guff;

	for (int n = 0; n < skip_pts; n++)
	{
		// populate the dsp filter guff vector
		guff.push_back(samples[skip_pts-n-1]);
		// apend some zeros to the samples
		samples.push_back(tcomplex(0,0));
	}

	samples.erase(samples.begin(), samples.begin() + skip_pts);

	for (int n = 0; n < skip_pts; n++)
		samples[n] = samples[n] - guff[n];

	// copy to fixed-size vector
	if( samples.size() != 0 ) 
	{
		cvm::cvector FID(samples.size());

		int n = 1;
		for( std::vector<tcomplex>::iterator it = samples.begin(); it != samples.end(); ++it, ++n ) 
			FID[n] = *it;
        
        // truncate FID if long, prob need a warning here...
        if ( FID.size() > 16384 )
            FID.resize(16384);

		m_fid.AppendFromVector(FID);
	}
}

/*
// find the average value at the end of the fid (last 5%)
int pts = np/2 - np/40;
float real_av;
float cplx_av;
for ( int n = np/40; n <= np/2; n++) {
real_av = real_av + m_fid.m_cvmFID(n).real()/pts;
cplx_av = cplx_av + m_fid.m_cvmFID(n).imag()/pts;
}

// subtract
std::complex<double> av(real_av, cplx_av);
for ( int n = 1; n <= np/2; n++) 
m_fid.m_cvmFID(n) = m_fid.m_cvmFID(n) - av;

*/

// silly Bruker lookup function
float tarquin::CFIDReaderBruker::lookup_skip_points( int decim, int dspfvs )
{
	switch ( dspfvs )
	{
		case 10 :
			switch ( decim )
			{
				case 2 : return 44.7500f;
				case 3 : return 33.5000f;
				case 4 : return 66.6250f;
				case 6 : return 59.0833f;
				case 8 : return 68.5625f;
				case 12 : return 60.3750f;
				case 16 : return 69.5313f;
				case 24 : return 61.0208f;
				case 32 : return 70.2578f;
				case 48 : return 61.3438f;
				case 64 : return 70.2578f;
				case 96 : return 61.5052f;
				case 128 : return 70.3789f;
				case 192 : return 61.5859f;
				case 256 : return 70.4395f;
				case 384 : return 61.6263f;
				case 512 : return 70.4697f;
				case 768 : return 61.6465f;
				case 1024 : return 70.4849f;
				case 1536 : return 61.6566f;
				case 2048 : return 70.4924f;
				default : return 0;
			}
		case 11 :
			switch ( decim )
			{
				case 2 : return 46.0000f;
				case 3 : return 36.5000f;
				case 4 : return 48.0000f;
				case 6 : return 50.1667f;
				case 8 : return 53.2500f;
				case 12 : return 69.5000f;
				case 16 : return 72.2500f;
				case 24 : return 70.1667f;
				case 32 : return 72.7500f;
				case 48 : return 70.5000f;
				case 64 : return 73.0000f;
				case 96 : return 70.6667f;
				case 128 : return 72.5000f;
				case 192 : return 71.3333f;
				case 256 : return 72.2500f;
				case 384 : return 71.6667f;
				case 512 : return 72.1250f;
				case 768 : return 71.8333f;
				case 1024 : return 72.0625f;
				case 1536 : return 71.9167f;
				case 2048 : return 72.0313f;
				default : return 0;
			}
		case 12 :
			switch ( decim )
			{
				case 2 : return 46.311f;
				case 3 : return 36.530f;
				case 4 : return 47.870f;
				case 6 : return 50.229f;
				case 8 : return 53.289f;
				case 12 : return 69.551f;
				case 16 : return 71.600f;
				case 24 : return 70.184f;
				case 32 : return 72.138f;
				case 48 : return 70.528f;
				case 64 : return 72.348f;
				case 96 : return 70.700f;
				case 128 : return 72.524f;
				default : return 0;
			}
		default : return 0;
	}
}


template <class T>
bool tarquin::CFIDReaderBruker::from_string(T& t, const std::string& s, std::ios_base& (*f)(std::ios_base&))
{
	std::istringstream iss(s);
	return !(iss >> f >> t).fail();
}

template<typename T>
std::string tarquin::CFIDReaderBruker::xxx_to_bin(const T& value)
{
	const std::bitset<std::numeric_limits<T>::digits + 1> bs(value);
	const std::string s(bs.to_string());
	const std::string::size_type pos(s.find_first_not_of('0'));
	return pos == std::string::npos ? "0" : s.substr(pos);
}


// this function was stolen from http://www.gamedev.net/reference/articles/article2091.asp
float tarquin::CFIDReaderBruker::FloatSwap( float f )
{
	union
	{
		float f;
		unsigned char b[4];
	} dat1, dat2;
	dat1.f = f;
	dat2.b[0] = dat1.b[3];
	dat2.b[1] = dat1.b[2];
	dat2.b[2] = dat1.b[1];
	dat2.b[3] = dat1.b[0];
	return dat2.f;
}

int tarquin::CFIDReaderBruker::LongSwap (int i)
{
	unsigned char b1, b2, b3, b4;

	b1 = i & 255;
	b2 = ( i >> 8 ) & 255;
	b3 = ( i>>16 ) & 255;
	b4 = ( i>>24 ) & 255;

	return ((int)b1 << 24) + ((int)b2 << 16) + ((int)b3 << 8) + b4;
}


short tarquin::CFIDReaderBruker::ShortSwap( short s )
{
	unsigned char b1, b2;

	b1 = s & 255;
	b2 = (s >> 8) & 255;

	return (b1 << 8) + b2;
}

