#include "CFIDReaderVarian.hpp"
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


tarquin::CFIDReaderVarian::CFIDReaderVarian(CFID& fid, CBoswell& log) : 
	CFIDReader(fid, log)
{
}


void tarquin::CFIDReaderVarian::Load(std::string strFilename, const Options& opts, CBoswell& log)
{
	// open the fid at the start
	std::ifstream fid_file(strFilename.c_str(), std::ios::in | std::ios::binary );
	if(!fid_file)
		throw Exception("couldn't open file '%s'", strFilename.c_str());

	// read data file header
	int nblocks;
	fid_file.read((char *) &nblocks, sizeof nblocks);
	nblocks = LongSwap(nblocks);
	int ntraces;
	fid_file.read((char *) &ntraces, sizeof ntraces);
	ntraces = LongSwap(ntraces);
	int np;
	fid_file.read((char *) &np, sizeof np);
	np = LongSwap(np);
	int ebytes;
	fid_file.read((char *) &ebytes, sizeof ebytes);
	ebytes = LongSwap(ebytes);
	int tbytes;
	fid_file.read((char *) &tbytes, sizeof tbytes);
	tbytes = LongSwap(tbytes);
	int bbytes;
	fid_file.read((char *) &bbytes, sizeof bbytes);
	bbytes = LongSwap(bbytes);
	short vers_id;
	fid_file.read((char *) &vers_id, sizeof vers_id);
	vers_id = ShortSwap(vers_id);
	short status;
	fid_file.read((char *) &status, sizeof status);
	status = ShortSwap(status);
	int nbheaders;
	fid_file.read((char *) &nbheaders, sizeof nbheaders);
	nbheaders = LongSwap(nbheaders);

	// read the bits from the status variable
	std::string str_bits = xxx_to_bin(status);
	int bit_num = str_bits.length();
	//char s_hyper = str_bits[bit_num-6];
	//char s_complex = str_bits[bit_num-5];
	char s_float = str_bits[bit_num-4];
	char s_32 = str_bits[bit_num-3];
	//char s_spec = str_bits[bit_num-2];
	//char sdata = str_bits[bit_num-1];

	/*std::cout << sdata << std::endl;
	  std::cout << s_spec << std::endl;
	  std::cout << s_32 << std::endl;
	  std::cout << s_float << std::endl;
	  std::cout << s_complex << std::endl;
	  std::cout << s_hyper << std::endl;*/

	// read the block header 
	short scale;
	fid_file.read((char *) &scale, sizeof scale);
	scale = ShortSwap(scale);
	short bstatus;
	fid_file.read((char *) &bstatus, sizeof bstatus);
	bstatus = ShortSwap(bstatus);
	short index;
	fid_file.read((char *) &index, sizeof index);
	index = ShortSwap(index);
	short mode;
	fid_file.read((char *) &mode, sizeof mode);
	mode = ShortSwap(mode);
	int ctcount;
	fid_file.read((char *) &ctcount, sizeof ctcount);
	ctcount = LongSwap(ctcount);
	float lpval;
	fid_file.read((char *) &lpval, sizeof lpval);
	lpval = FloatSwap(lpval);
	float rpval;
	fid_file.read((char *) &rpval, sizeof rpval);
	rpval = FloatSwap(rpval);
	float lvl;
	fid_file.read((char *) &lvl, sizeof lvl);
	lvl = FloatSwap(lvl);
	float tlt;
	fid_file.read((char *) &tlt, sizeof tlt);
	tlt = FloatSwap(tlt);

	// create the FID
	cvm::cvector FID(np/2);

	// now read in the first np points of the fid
	// format depends on the s_float and s_32 variables

	if ( s_32 == '1' ) 
	{
		std::vector<int> n(np, 0);
		fid_file.read((char *)&n[0], np*4);
		// change endian and populate FID
		for(int i = 0; i < np; i = i + 2)
			FID((i+2)/2) = tcomplex(LongSwap(n[i]), LongSwap(n[i+1]));
	}  
	else if ( s_float == '1' ) 
	{
		std::vector<float> n(np, 0);
		fid_file.read((char *)&n[0], np*4);
		// change endian and populate FID
		for(int i = 0; i < np; i = i + 2)
			FID((i+2)/2) = tcomplex(FloatSwap(n[i]), FloatSwap(n[i+1]));
	}  
	else 
	{
		std::vector<short> n(np, 0);
		fid_file.read((char *)&n[0], np*2);
		// change endian and populate FID
		for(int i = 0; i < np; i = i + 2)
			FID((i+2)/2) = tcomplex(ShortSwap(n[i]), ShortSwap(n[i+1]));
	}  

	// close the file
	fid_file.close();

	
	// read the sampling frequency and feild strength from the procpar file if available
	// strip off the last 3 charectors (assuming they are fid) and replace with procpar
	size_t len = strFilename.size();
	std::string str_procpar = strFilename.substr(0, len-3);
	str_procpar.append("procpar");

	// open the procpar text file for parsing
	std::string line;
	std::ifstream procpar ( str_procpar.c_str() );

    // number of scans, found for scaling the FID
	float ns;

	if ( procpar.is_open() )
	{
		while (! procpar.eof() )
		{
			getline (procpar,line);
			if ( line.substr(0,3) == "sw " )
			{
				// read the next line and trim the first two chars
				size_t len = line.size();
				getline (procpar,line);
				float fs;
				from_string<float>(fs, line.substr(2,len-2), std::dec);
				m_fid.SetSamplingFrequency(fs);
			} 
			else if ( line.substr(0,5) == "sfrq " )
			{
				// read the next line and trim the first two chars
				size_t len = line.size();
				getline (procpar,line);
				float ft;
				from_string<float>(ft, line.substr(2,len-2), std::dec);
				m_fid.SetTransmitterFrequency(ft*1e6);
			}
            else if ( line.substr(0,3) == "nt " )
			{
				// read the next line and trim the first two chars
				size_t len = line.size();
				getline (procpar,line);
				from_string<float>(ns, line.substr(2,len-2), std::dec);
			}
		}
		procpar.close();
	}
	else
		throw Exception("unable to open procpar file");


	//plot(m_fid.m_cvmFID);
	//std::cout << np << std::endl;
	//std::cout << m_fid.m_cvmFID.size() << std::endl;

    // find the average value at the end of the fid (last 5%)
    /*
	int pts = np/2 - np/40;
	float real_av = 0.0f;
	float cplx_av = 0.0f;
	for ( int n = np/2 - np/40; n <= np/2; n++) 
	{
		real_av = real_av + float(FID(n).real()/pts);
		cplx_av = cplx_av + float(FID(n).imag()/pts);
	}
    */

    // scale by number of scans 
	std::complex<double> nt(ns, 0);
    //std::cout << std::endl << nt << " scans for " << strFilename << std::endl;

	for ( int n = 1; n <= np/2; n++) 
		FID(n) = FID(n) / nt;

    // truncate FID if long, prob need a warning here...
    if ( FID.size() > 16384 )
        FID.resize(16384);

	m_fid.AppendFromVector(FID);

}

template <class T>
bool tarquin::CFIDReaderVarian::from_string(T& t, const std::string& s, std::ios_base& (*f)(std::ios_base&))
{
	std::istringstream iss(s);
	return !(iss >> f >> t).fail();
}

template<typename T>
std::string tarquin::CFIDReaderVarian::xxx_to_bin(const T& value)
{
	const std::bitset<std::numeric_limits<T>::digits + 1> bs(value);
	const std::string s(bs.to_string());
	const std::string::size_type pos(s.find_first_not_of('0'));
	return pos == std::string::npos ? "0" : s.substr(pos);
}


// this function was stolen from http://www.gamedev.net/reference/articles/article2091.asp
float tarquin::CFIDReaderVarian::FloatSwap( float f )
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

int tarquin::CFIDReaderVarian::LongSwap (int i)
{
	unsigned char b1, b2, b3, b4;

	b1 = i & 255;
	b2 = ( i >> 8 ) & 255;
	b3 = ( i>>16 ) & 255;
	b4 = ( i>>24 ) & 255;

	return ((int)b1 << 24) + ((int)b2 << 16) + ((int)b3 << 8) + b4;
}


short tarquin::CFIDReaderVarian::ShortSwap( short s )
{
	unsigned char b1, b2;

	b1 = s & 255;
	b2 = (s >> 8) & 255;

	return (b1 << 8) + b2;
}

