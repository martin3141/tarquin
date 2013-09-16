#include "CFIDReaderPhilips.hpp"
#include "CBoswell.hpp"
#include "CFID.hpp"

#include <string>
#include <sstream>


tarquin::CFIDReaderPhilips::CFIDReaderPhilips(CFID& fid, CBoswell& log)
: CFIDReader(fid, log)
{
    m_FFT = false;
}

// this works so that the users can select either the sdat or the spar file
void tarquin::CFIDReaderPhilips::Load(std::string strFilename, const Options& opts, CBoswell& log)
{
	std::string strExt  = GetFilenameExtension(strFilename);
	std::string strBase = GetFilenameBase(strFilename);

	// name of parameter file
	std::string strFileParam = strBase + ".spar";

	// can we open this, or do we need to change case of extension?
	if( false == FileExists(strFileParam) ) 
	{
		strFileParam = strBase + ".SPAR";

		if( false == FileExists(strFileParam) )
			throw Exception("SPAR file '%s' does not exist", strFileParam.c_str());
	}

	// parse the parameters and populate the useful data fields
	ReadParamFile(strFileParam);

	// name of data file
	std::string strFileData = strBase + ".sdat";

	// can we open this, or do we need to change case of extension?
	if( false == FileExists(strFileData) ) 
	{
		strFileData = strBase + ".SDAT";

		if( false == FileExists(strFileData) )
			throw Exception("SDAT file '%s' does not exist", strFileData.c_str());
	}

	// get the FID data
	ReadFIDFile(strFileData, opts);
}

// parse an SPAR file - turn into a list of tokens, then pull out fields
// TODO: can't we just use the default Lexer?
void tarquin::CFIDReaderPhilips::ReadParamFile(std::string strFileParam)
{
	std::ifstream file;

	// attempt to open
	file.open(strFileParam.c_str());

	if( true == file.fail() )
		throw Exception("failed to open Phillips parameter file '%s'", strFileParam.c_str());

	// turn file into a list of tokens, except for comments
	size_t nLine = 1;
	for( std::string strLine; getline(file, strLine); ++nLine ) 
	{
		// if empty, skip line
		if( strLine.length() == 0)
			continue;

		// if commented, skip line
		if( strLine.at(0) == '!' )
			continue;

		// if \r, skip line
		if( strLine.at(0) == '\r' )
			continue;

        //std::cout << strLine << std::endl;
        size_t found;
        found = strLine.find(":",0);
        
        // skip if a colon isn't found
        if ( found == std::string::npos)
			continue;

		std::string strKey   = strLine.substr(0, found);
		std::string strValue = strLine.substr(strKey.size()+1, strLine.size() - strKey.size()-1 );

		LooseWhiteSpace(strKey);

		// add to collection of tokens
		m_tokens.push_back(std::make_pair(strKey, nLine));
		m_tokens.push_back(std::make_pair(strValue, nLine));
	}

	// turn these tokens into parameter values
	EatTokens();
}


/*
 * The Vax conversion routines require that this copyright notice is included.
 *
 Copyright(c) 1982 Association of Universities for Research in Astronomy Inc.

 The FOCAS software is publicly available, but is NOT in the public domain.
 The difference is that copyrights granting rights for unrestricted use and
 redistribution have been placed on all of the software to identify its authors.
 You are allowed and encouraged to take this software and use it as you wish,
 subject to the restrictions outlined below.

 Permission to use, copy, modify, and distribute this software and its
 documentation is hereby granted without fee, provided that the above
 copyright notice appear in all copies and that both that copyright notice
 and this permission notice appear in supporting documentation, and that
 references to the Association of Universities for Research in Astronomy
 Inc. (AURA), the National Optical Astronomy Observatories (NOAO), or the
 Faint Object Classification and Analysis System (FOCAS) not be used in
 advertising or publicity pertaining to distribution of the software without
 specific, written prior permission from NOAO.  NOAO makes no
 representations about the suitability of this software for any purpose.  It
 is provided "as is" without express or implied warranty.

 NOAO DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE, INCLUDING ALL
 IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS, IN NO EVENT SHALL NOAO
 BE LIABLE FOR ANY SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION
 OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN 
 CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 */
inline short int vax2s(short int le_short)
{
	union	{
		short	bs;
		char	bc[2];
	} bsbuf;
	char	chold;

	bsbuf.bs = le_short;
	chold = bsbuf.bc[0];
	bsbuf.bc[0] = bsbuf.bc[1];
	bsbuf.bc[1] = chold;

	return (bsbuf.bs);
}

// Vax float to IEEE 32-bit float
inline float vax2f(unsigned char le_flt[4])
{
	union {
		unsigned char bytes[4];
		float  ff;
		short fs[2];
		struct {				/* VAX F-floating */
			unsigned int sign:1;
			unsigned int exponent:8;
			unsigned int mantissa:23;
		} v;
	} vaxbuf;
	float	ieeebuf;

	vaxbuf.bytes[0] = le_flt[3];
	vaxbuf.bytes[1] = le_flt[2];
	vaxbuf.bytes[2] = le_flt[1];
	vaxbuf.bytes[3] = le_flt[0];

	vaxbuf.fs[0] = vax2s (vaxbuf.fs[0]);
	vaxbuf.fs[1] = vax2s (vaxbuf.fs[1]);

	// this doesn't appear to be doing the right thing - so I took it out - GMR
	//if (vaxbuf.v.exponent < 3) 		/* prevent underflow */ {
	//			ieeebuf = 0.0;
	//			std::cout << std::endl << "underflow detected - result is zero";
	//		}
	//		else
	ieeebuf = vaxbuf.ff / 4.0f;

	return (ieeebuf);
}


void tarquin::CFIDReaderPhilips::EatTokens()
{
    // useful for mega-press
    int rows = 1;

	for(TokenList::iterator it = m_tokens.begin(); it != m_tokens.end(); it+=2 ) 
	{
		std::string strKey   = it->first;
		std::string strValue = (it+1)->first;

		// magic conversion
		std::istringstream strmValue;
		strmValue.clear();
		strmValue.str(strValue);

		// crude, but I have yet to find a more elegant way to do this
		if( 0 == strKey.compare("sample_frequency") ) 
		{
			treal fs;
			strmValue >> fs;

			m_fid.SetSamplingFrequency(fs);				
		}
		else if( 0 == strKey.compare("synthesizer_frequency") ) 
		{
			treal ft;
			strmValue >> ft;

			m_fid.SetTransmitterFrequency(ft);
		}
        else if( 0 == strKey.compare("rows") ) 
		{
			strmValue >> rows;
		}
		else if( 0 == strKey.compare("averages") ) 
		{
			int nAverages;
			strmValue >> nAverages;

			m_fid.SetAverages(nAverages);
		}
        else if( 0 == strKey.compare("echo_time") ) 
		{
            if (not(m_fid.IsKnownEchoTime()))
            {
                treal nEcho;
                strmValue >> nEcho;
                m_fid.SetEchoTime(nEcho/1000);
            }
		}
        else if( 0 == strKey.compare("dim1_ext") || 0 == strKey.compare("SUN_dim1_ext") ) 
        {
            std::string ext;
            strmValue >> ext;
            if ( 0 == ext.compare("[ppm]") )
                m_FFT = true;
        }
        else if( 0 == strKey.compare("dim1_pnts") || 0 == strKey.compare("SUN_dim1_pnts") ) 
        {
            int pts;
            strmValue >> pts;
            m_fid.SetNumberOfPoints(pts);
        }
        else if( 0 == strKey.compare("dim2_pnts") || 0 == strKey.compare("SUN_dim2_pnts") ) 
        {
            int dim_rows;
            strmValue >> dim_rows;
			m_fid.SetRows(dim_rows);
        }
        else if( 0 == strKey.compare("dim3_pnts") || 0 == strKey.compare("SUN_dim3_pnts") ) 
        {
            int cols;
            strmValue >> cols;
            m_fid.SetCols(cols);
        }
    }

    // looks like a dynamic scan
    if ( m_fid.GetRows() == 1 && m_fid.GetCols() == 1 && rows != 1 )
    {
        m_fid.SetCols(rows);
        m_fid.SetDyn(true);
    }

}

// read the complex FID from the file
// maybe there is a better way to do this - is there a field in the param file that
// dictates the number of samples?
void tarquin::CFIDReaderPhilips::ReadFIDFile(std::string strFilename, const Options& opts)
{
    std::ifstream file(strFilename.c_str(), std::ios::binary);

    if( true == file.fail() )
        throw Exception("failed to read Phillips data file: '%s'", strFilename.c_str());

    //if ( opts.GetDynAv() == ALL || opts.GetDynAv() == ODD || opts.GetDynAv() == EVEN || opts.GetDynAv() == SUBTRACT )
    if ( false )
    {
        int n = 0;
        std::vector<tcomplex> samples;
        //while( n < m_fid.GetNumberOfPoints() ) 
        while( false == file.eof() ) 
        {
            unsigned char real_part_vax[4] = {0};
            unsigned char imag_part_vax[4] = {0};
            file.read((char*)real_part_vax, 4);
            file.read((char*)imag_part_vax, 4);

            // eof is not true until we have actually tried to read past it
            if( true == file.eof() )
                break;

            // convert from Vax floats to IEEE floats (see "common.hpp" for credits)
            float real_part_ieee = vax2f(real_part_vax);
            float imag_part_ieee = vax2f(imag_part_vax);	

            tcomplex z(real_part_ieee, -imag_part_ieee);

            // store this sample

            div_t divresult;
            divresult = div(n, m_fid.GetNumberOfPoints());
            if ( ( divresult.quot % 2 == 0 ) && ( opts.GetDynAv() == ODD ) )
            {
                if ( ( n % m_fid.GetNumberOfPoints() ) < int(samples.size()) )
                {
                    samples[n % m_fid.GetNumberOfPoints()] += z;
                }
                else
                {
                    samples.push_back(z);
                }
            }
            else if ( ( divresult.quot % 2 == 1 ) && ( opts.GetDynAv() == EVEN ) )
            {
                if ( ( n % m_fid.GetNumberOfPoints() ) < int(samples.size()) )
                {
                    samples[n % m_fid.GetNumberOfPoints()] += z;
                }
                else
                {
                    samples.push_back(z);
                }
            }
            else if ( opts.GetDynAv() != EVEN && opts.GetDynAv() != ODD )
            {
                // should the even scans be negated?
                if ( divresult.quot % 2 == 1 && ( opts.GetDynAv() == SUBTRACT ) )
                    z = -z;

                if ( ( n % m_fid.GetNumberOfPoints() ) < int(samples.size()) )
                {
                    samples[n % m_fid.GetNumberOfPoints()] += z;
                }
                else
                {
                    samples.push_back(z);
                }
            }

            ++n;
        }
        
        int dyn_scans;
        if ( opts.GetDynAv() != EVEN && opts.GetDynAv() != ODD )
            dyn_scans = n/m_fid.GetNumberOfPoints();
        else
            dyn_scans = n/m_fid.GetNumberOfPoints()/2;

        // copy to fixed-size vector
        if( samples.size() != 0 ) 
        {
            cvm::cvector FID(samples.size());

            int n = 1;
            for( std::vector<tcomplex>::iterator it = samples.begin(); it != samples.end(); ++it, ++n ) 
                FID[n] = *it;
            
            // vector is cleared because mega-press data gets loaded twice
            // when using the gui
            m_fid.ClearVector();

            if ( m_FFT )
            {
                FID = fftshift(FID);
                FID = ifft(FID);
            }

            m_fid.AppendFromVector(FID/dyn_scans);
        }
        m_fid.SetCols(1);
    }
    //else if ( opts.GetDynAv() == NONE )
    else if ( true )
    {
        for ( int m = 0; m < m_fid.GetVoxelCount(); m++ )
        {
            int n = 0;
            std::vector<tcomplex> samples;
            //while( false == file.eof() ) 
            while( n < m_fid.GetNumberOfPoints() ) 
            {
                unsigned char real_part_vax[4] = {0};
                unsigned char imag_part_vax[4] = {0};
                file.read((char*)real_part_vax, 4);
                file.read((char*)imag_part_vax, 4);

                // eof is not true until we have actually tried to read past it
                if( true == file.eof() )
                    break;

                // convert from Vax floats to IEEE floats (see "common.hpp" for credits)
                float real_part_ieee = vax2f(real_part_vax);
                float imag_part_ieee = vax2f(imag_part_vax);	

                tcomplex z(real_part_ieee, -imag_part_ieee);

                // store this sample
                samples.push_back(z);
                ++n;
            }

            // copy to fixed-size vector
            if( samples.size() != 0 ) 
            {
                cvm::cvector FID(samples.size());

                int n = 1;
                for( std::vector<tcomplex>::iterator it = samples.begin(); it != samples.end(); ++it, ++n ) 
                    FID[n] = *it;

                if ( m_FFT )
                {
                    FID = fftshift(FID);
                    FID = ifft(FID);
                }

                m_fid.AppendFromVector(FID);
            }
        }
    }
    else
        throw Exception("Dynamic averaging scheme incompatible with this data format");
}

