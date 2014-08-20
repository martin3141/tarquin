#include "CFIDReaderGE.hpp"
#include "CBoswell.hpp"
#include "CFID.hpp"
#include <boost/filesystem.hpp>

#include <string>
#include <sstream>
#include <fstream>

/* Notes on the GE format...
GE P files contain time domain data for SVS and k-space data for CSI.  Each average is stored in the file for SVS, not sure about CSI.  Water reference data is also stored in the file when avaialble for SVS, not sure about CSI.

Basic file structure:

----------bof
Header data
----------data offset
MRS data
----------eof

Every p-files has an rdb_header_rev number.  Every rdb_header_rev number has an associated data offset for early versions < 11, for later versions the offset can be read from the header. The data from each coil is stored seperatly, before the acquired data frames there should be one prescan frame of just noise.  This may be handy for correctly scaling the data. The MRS data is stored as follows:

----------data offset + [(nCoil=0)*scans - 1]*nFieldSize*2*pt_size
nCoil=0 prescan data, one frame
----------data offse + (nCoil=0)*scans*nFieldSize*2*pt_size
water unsupressed data one frame per average
----------data offset + [(nCoil=0)*scans + w_scans]*nFieldSize*2*pt_size 
water supressed data one frame per average
----------data offset + [(nCoil=1)*scans - 1]*nFieldSize*2*pt_size 
nCoil=1 prescan data, one frame
----------data offset + (nCoil=1)*scans*FieldSize*2*pt_size 
water unsupressed data one frame per average
----------data offset + [(nCoil=1)*scans + w_scans]*nFieldSize*2*pt_size 
water supressed data one frame per average
----------data offset + [(nCoil=2)*scans - 1]*nFieldSize*2*pt_size
nCoil=2 prescan data, one frame
and so on...
----------eof

data offset     p-file data offset
nCoil           coil number from 0 to coils-1 
w_scans         number of water reference scans = rows*cols*slices*w_averages
ws_scans        number of water supressed scans = rows*cols*slices*w_averages
scans           w_scans + ws_scans + 1 (+ 1 is needed for the prescan frame)
nFieldSize      complex points per fid
pt_size         sizeof(int)

Question - does single coil data still contain a prescan?

MW Apr 2011
*/

tarquin::CFIDReaderGE::CFIDReaderGE(CFID& fid, CBoswell& log) : 
	CFIDReader(fid, log)
{
}

void tarquin::CFIDReaderGE::Load(std::string strFilename, const Options& opts, CBoswell& log)
{
	DiscoverOptions(strFilename, log);
	LoadFromOptions(strFilename, opts, true);
}

void tarquin::CFIDReaderGE::LoadW(std::string strFilename, const Options& opts, CBoswell& log)
{
	DiscoverOptions(strFilename, log);
	LoadFromOptions(strFilename, opts, false);
}

void tarquin::CFIDReaderGE::DiscoverOptions(std::string strFilename, CBoswell& log) 
{
	//std::cout << std::endl << "GE p-file info" << std::endl;
	std::size_t nFileSize = GetFileSize(strFilename);

	if( 0 == nFileSize )
		throw Exception("the file '%s' is 0 bytes in length", strFilename.c_str());

	// attempt to open the file
	std::ifstream file(strFilename.c_str(), std::ios::binary);

	if( true == file.fail() )
		throw Exception("failed to open '%s'", strFilename.c_str());

	// find the rdb version number of the p-file
	file.seekg(0, std::ios_base::beg);
	float rdb_header_rev;
	file.read((char*)&rdb_header_rev, 4);

	int p_file_off = 0;
    bool swap_end = false;
	// find the p-file offset 
    if ( ( FloatSwap(rdb_header_rev) > 5 ) && ( FloatSwap(rdb_header_rev) < 6 ))
    {
        swap_end = true;
        p_file_off = 39940;
    }
    else if ( (int) FloatSwap(rdb_header_rev) == 7 )
    {
        swap_end = true;
        p_file_off = 39984;
    }
    else if ( (int) FloatSwap(rdb_header_rev) == 8 )
    {
        swap_end = true;
        p_file_off = 60464;
    }
	else if ( (int) rdb_header_rev == 9 )
    {
		p_file_off = 61464;
    }
	else if ( (int) rdb_header_rev == 11 )
    {
		p_file_off = 66072;
    }
	else if ( (int) rdb_header_rev > 11 )
	{
		file.seekg(1468, std::ios_base::beg);
		file.read((char*)&p_file_off, 4);
	}
	else
		throw Exception("Problem with the rdb header number", strFilename.c_str());

    if ( swap_end )
        rdb_header_rev = FloatSwap(rdb_header_rev);

	//std::cout << std::endl << "Rdb header rev      : " << rdb_header_rev << std::endl;

	log.LogMessage(LOG_INFO, "Rdb reader rev : %f",rdb_header_rev);
	log.LogMessage(LOG_INFO, "P-file offset  : %i",p_file_off);

	//std::cout << "P-file offset : " << p_file_off << std::endl;

	// find the offset to the data	
	short frame_size;
	file.seekg(80, std::ios_base::beg);
	file.read((char*)&frame_size, 2);
    if ( swap_end )
        frame_size = ShortSwap(frame_size);

	// convert to bits ?
    int frame_size_int;
	frame_size_int = frame_size * 2 * 4;
	log.LogMessage(LOG_INFO, "Frame size : %i",frame_size_int);
    int ver_offset = 0;
    if ( rdb_header_rev > 14 ) // could this be required for lower versions?
        ver_offset = frame_size_int;

	int data_offset = p_file_off + ver_offset;
	log.LogMessage(LOG_INFO, "Offset to data : %i",data_offset);

    // find the number of echoes
	short necho;
	file.seekg(70, std::ios_base::beg);
	file.read((char*)&necho, 2);
    if ( swap_end )
        necho = ShortSwap(necho);
	
    log.LogMessage(LOG_INFO, "Nechoes : %i",necho);

	// find the number of data frames	
	short nframes;
	file.seekg(74, std::ios_base::beg);
	file.read((char*)&nframes, 2);
    if ( swap_end )
        nframes = ShortSwap(nframes);

	// find number of water reference frames
	float num_ref_frames;
	file.seekg(292, std::ios_base::beg);
	file.read((char*)&num_ref_frames, 4);
    if ( swap_end )
        num_ref_frames = FloatSwap(num_ref_frames);

    
	if( num_ref_frames > 32768 || num_ref_frames < 0 )
        throw Exception("Unrealistic number of reference frames, check data type.");

	log.LogMessage(LOG_INFO, "Water ref frames : %f",num_ref_frames);
	//std::cout << "Water ref frames : " << num_ref_frames << std::endl;
	float num_sig_frames = nframes - num_ref_frames;

	if( num_sig_frames > 32768 || num_sig_frames < 0 )
        throw Exception("Unrealistic number of reference frames, check data type.");
	
    log.LogMessage(LOG_INFO, "Water sup frames : %f",num_sig_frames);
	//std::cout << "Water supressed frames : " << num_sig_frames << std::endl;

	// find number of coils
	std::vector<short> coil_array;
	short coil_tmp;
	file.seekg(200, std::ios_base::beg);
	for ( size_t n = 0; n < 8; n++ )
	{
		file.read((char*)&coil_tmp, 2);
        if ( swap_end )
            coil_tmp = ShortSwap(coil_tmp);
		coil_array.push_back(coil_tmp);
	}

	int coils = 0;
	for ( size_t n = 0; n < 8; n = n + 2 )
		if ( ( coil_array[n] != 0 ) || ( coil_array[n+1] != 0 ) )
			coils += (int) coil_array[n+1] - (int) coil_array[n] + 1;

	if ( coils == 0 )
		coils = 1;

    log.LogMessage(LOG_INFO, "Number of coils : %i", coils);
    //std::cout << coils << std::endl;

	// find the sampling frequency 
	float fs = 0;
	file.seekg(368, std::ios_base::beg);
	file.read((char*)&fs, 4);
    if ( swap_end )
        fs = FloatSwap(fs);

	//std::cout << "Sampling freq : " << fs << std::endl;

	short csi_dims = 0;
	file.seekg(372, std::ios_base::beg);
	file.read((char*)&csi_dims, 2);
    if ( swap_end )
        csi_dims = ShortSwap(csi_dims);

	short xcsi = 0;
	file.read((char*)&xcsi, 2);
    if ( swap_end )
        xcsi = ShortSwap(xcsi);

    log.LogMessage(LOG_INFO, "x dim : %i",xcsi);
    //std::cout << "xcsi : " << xcsi << std::endl;
    short ycsi = 0;
	file.read((char*)&ycsi, 2);
    if ( swap_end )
        ycsi = ShortSwap(ycsi);

    log.LogMessage(LOG_INFO, "y dim : %i",ycsi);
    //std::cout << "ycsi : " << ycsi << std::endl;
    short zcsi = 0;
	file.read((char*)&zcsi, 2);
    if ( swap_end )
        zcsi = ShortSwap(zcsi);

    log.LogMessage(LOG_INFO, "z dim : %i",zcsi);
    //std::cout << "zcsi : " << zcsi << std::endl;

	// find the field strength
	int ft = 0;
	file.seekg(424, std::ios_base::beg);
	file.read((char*)&ft, 4);
    if ( swap_end )
        ft = LongSwap(ft);
	//std::cout << "Trans freq : " << ft << std::endl;

    int TE = 0;

    //int image_offset = 0;
    // this info is from SIVIC svkGEPFileReader.cc
    if ( rdb_header_rev > 6.95 && rdb_header_rev < 8.0 )
    {
        file.seekg(39148, std::ios_base::beg);
        file.read((char*)&TE, 4);
    }
    else if ( rdb_header_rev == 9.0 )
    {
        file.seekg(60136, std::ios_base::beg);
        file.read((char*)&TE, 4);
    }
    else if ( rdb_header_rev == 11.0 )
    {
        file.seekg(65032, std::ios_base::beg);
        file.read((char*)&TE, 4);
    }
    else if ( rdb_header_rev > 11.0 )
    {
        // comes from noise_weights_te.m by Fred J. Frigo
        /*
        file.seekg(1504, std::ios_base::beg);
        file.read((char*)&image_offset, 4);
        file.seekg(image_offset + 1064, std::ios_base::beg);
        file.read((char*)&TE, 4);
        */

        // comes from my own hacking about
        file.seekg(4*303, std::ios_base::beg);
        file.read((char*)&TE, 4);
    }

    if ( swap_end )
        TE = LongSwap(TE);

    /*
    else if ( rdb_header_rev == 14.0 )
        file.seekg(144580, std::ios_base::beg);
    else if ( rdb_header_rev == 15.0 )
        file.seekg(144580, std::ios_base::beg);
    else if ( rdb_header_rev == 20.0 )
        file.seekg(148404, std::ios_base::beg);
    */

    // find out the offset to rdb_hdr_image
    /*
    int image_offset = 0;
    int TE = 0;
	file.seekg(1504, std::ios_base::beg);
    file.read((char*)&image_offset, 4);
	file.seekg(image_offset + 1064, std::ios_base::beg);
    file.read((char*)&TE, 4);
    if ( swap_end )
        TE = LongSwap(TE);
        */

    
	//std::cout << "TE : " << TE << std::endl;

	// convert to seconds
	treal te_float = TE * 1e-6;
	//std::cout << "Echo time : " << te_float << std::endl;

	// put values in the appropriate structures	
	m_fid.SetEchoTime(te_float);
	m_fid.SetSamplingFrequency(fs);
	m_fid.SetTransmitterFrequency(ft / 10);
    m_fid.SetRows(xcsi); // TODO check these
    m_fid.SetCols(ycsi); // TODO check these
    m_fid.SetSlices(zcsi);
	m_options.nOffset      = data_offset;
	m_options.nCoils       = coils;
	m_options.nFieldSize   = frame_size_int / 8;
	m_options.nWaterFrames = static_cast<int>(num_ref_frames);
	m_options.nWSFrames    = static_cast<int>(num_sig_frames);
	m_options.nHeaderRev   = rdb_header_rev;
	m_options.nEchoes      = necho;
}

void tarquin::CFIDReaderGE::LoadSHF(std::string strFilename, const Options& opts, bool WS, CBoswell& log)
{
	// use the boost tokenizer that keeps things like strings in quotes intact
	// middle argument is token separator 
	boost::escaped_list_separator<char> sep("\\", " \t\r", "\"");

	Lex(strFilename, sep);

	EatTokens(log); 

	// this just assumes that the data file is the same name as the SHF
	// file but without the name

	std::string strFilenameFull = GetFilenameBase(strFilename);
    // changed by MW 06/10/09 to correct assumption for phased array data
	std::string strFilenameShort = m_fid.GetFilename(); 

    int pos = strFilenameFull.find(strFilenameShort);
    
    //std::cout << strFilenameFull.substr(0,pos-1) + strFilenameShort << std::endl; 

	LoadFromOptions(strFilenameFull.substr(0,pos) + strFilenameShort, opts, WS);
}

void tarquin::CFIDReaderGE::EatTokens(CBoswell& log)
{
	// we want the number of points in the time domain
	std::string strDomain = "None";

	for( TokenList::iterator it = m_tokens.begin(); it+1 != m_tokens.end(); ++it ) 
	{
		std::string strKey = it->first;


		// value may actually be the next key, unless we we recognise the key
		std::istringstream strmValue;
		strmValue.clear();
		strmValue.str((it+1)->first);

		// options about reading the file
		if( strKey == "offset_to_data" ) 
		{
			strmValue >> m_options.nOffset;
	        log.LogMessage(LOG_INFO, "Offset to data : %i",m_options.nOffset);
		}
		else if( strKey == "signa_rdb_hdr_rev" ) 
        {
            strmValue >> m_options.nHeaderRev;
	        log.LogMessage(LOG_INFO, "Rdb reader rev : %f",m_options.nHeaderRev);
        }
		else if( strKey == "domain" ) 
		{
			strDomain = strmValue.str();
		}
		else if( strKey == "num_sig_frames" ) 
		{
			strmValue >> m_options.nWSFrames;
            log.LogMessage(LOG_INFO, "Water sup frames : %i",m_options.nWaterFrames);
		}
		else if( strKey == "num_ref_frames" ) 
		{
			strmValue >> m_options.nWaterFrames;
	        log.LogMessage(LOG_INFO, "Water ref frames : %i",m_options.nWaterFrames);
		}
		else if( strKey == "num_pts" ) 
		{
			if( strDomain == "Time" ) 
			{
				strmValue >> m_options.nFieldSize;
	            log.LogMessage(LOG_INFO, "Number of pts    : %i",m_options.nFieldSize);
			}
		}
		else if( strKey == "num_coils" ) 
		{
            //std::cout << strKey << std::endl;
            //std::cout << strmValue << std::endl;
			strmValue >> m_options.nCoils;
	        log.LogMessage(LOG_INFO, "Number of coils : %i",m_options.nCoils);
		}

        else if( strKey == "data_format" ) 
		{
		    strmValue.str((it+2)->first);
			strmValue >> m_options.nFormat;
		}

		// options about parameters, etc.
		else if( strKey == "center_freq" ) 
		{
			treal ft;
			strmValue >> ft;
			m_fid.SetTransmitterFrequency(ft * 1e6);
		}
		else if( strKey == "width" ) 
		{
			if( strDomain == "Time" ) 
			{
				treal fs;
				strmValue >> fs;
				m_fid.SetSamplingFrequency(fs);
                log.LogMessage(LOG_INFO, "Sampling frequency: %f",fs);
			}
		}
		else if( strKey == "TE" ) 
		{
			treal tau;
			strmValue >> tau;
			// Apparently, time is always quoted in ms.
			m_fid.SetEchoTime(tau * 1e-3);
            log.LogMessage(LOG_INFO, "Echo time    : %f",tau*1e-3);
		}
		else if( strKey == "data_file_name" ) 
		{
            std::string f_name;
			strmValue >> f_name;
            m_fid.SetFilename(f_name);
		}
		/*else if( strKey == "ppm_offset" ) 
		{
			treal ref;
			strmValue >> ref;
			m_fid.SetPPMRef(ref);
		}*/


	}
    
    //m_fid.SetRows(16); // TODO check these
    //m_fid.SetCols(16); // TODO check these
    //m_fid.SetSlices(1);

}

void tarquin::CFIDReaderGE::LoadFromOptionsCSI(std::string strFilename) 
{
    // attempt to open the file
	std::ifstream file(strFilename.c_str(), std::ios::binary);
	if( true == file.fail() )
		throw Exception("failed to open file: '%s'", strFilename.c_str());
    
	std::size_t nOffset      = m_options.nOffset;
	std::size_t nWaterFrames = m_options.nWaterFrames;
	std::size_t nWSFrames    = m_options.nWSFrames;
	std::size_t nFieldSize   = m_options.nFieldSize;
	std::size_t nCoils       = m_options.nCoils;
	std::size_t nFrames      = nWaterFrames + nWSFrames + 1;
	std::string format       = m_options.nFormat;
	float header_rev         = m_options.nHeaderRev;
    //int voxels               = m_fid.GetRows() * m_fid.GetCols() * m_fid.GetSlices();
    //int ws_avgs              = nWSFrames / voxels;
    //int w_avgs               = nWaterFrames / voxels;

	// flag that water ref data is included in the file	
	if ( nWaterFrames > 0 )
		m_fid.SetCWF(true);
	else
		m_fid.SetCWF(false);

    //std::cout << "nOffset      : " << nOffset << std::endl;
    //std::cout << "nWaterFrames : " << nWaterFrames << std::endl;
    //std::cout << "nWSFrames    : " << nWSFrames << std::endl;
    //std::cout << "nFieldSize   : " << nFieldSize << std::endl;
    //std::cout << "nCoils       : " << nCoils << std::endl;
    //std::cout << "nFrames      : " << nFrames << std::endl;
    //std::cout << "Rows         : " << m_fid.GetRows() << std::endl;
    //std::cout << "Cols         : " << m_fid.GetCols() << std::endl;
    //std::cout << "Slices       : " << m_fid.GetSlices() << std::endl;
    //std::cout << "WS avgs      : " << ws_avgs << std::endl;
    //std::cout << "W avgs       : " << w_avgs << std::endl;

    int pt_size = sizeof(int);

    // read in the first data point for each voxel for each coil
    // for phasing purposes (can this just be for each coil instead?)
    // assumes one scan only 2D CSI with no water ref data

    cmat_stdvec coil_vec;  // one row x col mat for each coil
    // for each coil
    for ( size_t coil = 0; coil < nCoils; coil++ )
    {
        cvm::cmatrix kspace(m_fid.GetRows(), m_fid.GetCols());
        std::size_t coil_loc = nOffset + coil*nFrames*nFieldSize*2*pt_size;
        cvm::cmatrix pt(m_fid.GetRows(), m_fid.GetCols());
        // for each column
        for ( int col = 0; col < m_fid.GetCols(); col++ )
        {
            // for each row
            for ( int row = 0; row < m_fid.GetRows(); row++ )
            {
                std::size_t vox_loc; 
                vox_loc = coil_loc + ( col + row * m_fid.GetRows() ) * nFieldSize*2*pt_size;
                file.seekg(vox_loc, std::ios_base::beg);
                // read the first data point
                int sampleReal;
                int sampleImag;
                file.read((char*)&sampleReal, 4);
                file.read((char*)&sampleImag, 4);
                // byteswap where necessary TODO move this outside the loop?
                if ( header_rev < 9 )
                {
                    sampleReal = LongSwap(sampleReal);
                    sampleImag = LongSwap(sampleImag);
                }
                pt(row + 1,col + 1) = tcomplex(sampleReal, sampleImag);
            }
        }
        coil_vec.push_back(pt);
    }
    //std::cout << coil_vec[0](1,1) << std::endl;

    // Transform the first point in k-space to get the first point in
    // the time domain for each coil. 
    // Use this for phasing the data to combine the coils.
    cmat_stdvec phase_mat_vec;  // one row x col mat for each coil
    for ( std::size_t n = 0; n < coil_vec.size(); n++ )
    {
        cvm::cmatrix tmp_mat(coil_vec[n].msize(), coil_vec[n].nsize());
        tmp_mat = fft(coil_vec[n]);
        // transpose the data
        tmp_mat = ~tmp_mat;
        tmp_mat.assign_imag(-tmp_mat.imag());
        // perform the second fft
        tmp_mat = fft(tmp_mat);
        // transpose back
        tmp_mat = ~tmp_mat;
        tmp_mat.assign_imag(-tmp_mat.imag());
        // do fftshift
        tmp_mat = fftshift(tmp_mat);
        // find the phase correction term
        cvm::cmatrix phase_mat(coil_vec[n].msize(), coil_vec[n].nsize());
        for ( int p = 1; p < coil_vec[n].msize() + 1; p++ )
            for ( int q = 1; q < coil_vec[n].nsize() + 1; q++ )
                phase_mat(p,q) = exp(-tcomplex(0,1)*arg(tcomplex(tmp_mat(p,q))));

        phase_mat_vec.push_back(phase_mat);
    }
    
    /* useful for testing
    std::cout << phase_mat_vec[0] << std::endl;
    std::cout << data_mat_vec.size() << std::endl;
    std::cout << data_mat_vec[0].msize() << std::endl;
    std::cout << data_mat_vec[0].nsize() << std::endl;
    cvm::rmatrix plot_mat = data_mat_vec[0].real();
    cvm::rvector plot_vec(plot_mat(8));
    plot(plot_vec);
    */
    
    // the final data structure
    cmat_stdvec data_mat_vec;
    // fill with zeros
    cvm::cmatrix data(m_fid.GetRows(), m_fid.GetCols());
    data.vanish();
    for ( size_t n = 0; n < nFieldSize; n++ )
        data_mat_vec.push_back(data);

    // for each coil
    for ( size_t coil = 0; coil < nCoils; coil++ )
    {
        cvm::cmatrix kspace(m_fid.GetRows(), m_fid.GetCols());
        std::size_t coil_loc = nOffset + coil*nFrames*nFieldSize*2*pt_size;
        // for each data point
        for( std::size_t n = 0; n < nFieldSize; n++ ) 
        {
            cvm::cmatrix tmp_mat(m_fid.GetRows(), m_fid.GetCols());
            // for each col
            for ( int col = 0; col < m_fid.GetCols(); col++ )
            {
                // for each row
                for ( int row = 0; row < m_fid.GetRows(); row++ )
                {
                    int data_pt = coil_loc + (( col + row*m_fid.GetRows() )*nFieldSize + n)*2*pt_size;
                    file.seekg(data_pt, std::ios_base::beg);
                    int sampleReal;
                    int sampleImag;
                    file.read((char*)&sampleReal, 4);
                    file.read((char*)&sampleImag, 4);

                    // byteswap where necessary
                    if ( header_rev < 9 )
                    {
                        sampleReal = LongSwap(sampleReal);
                        sampleImag = LongSwap(sampleImag);
                    }
                    tmp_mat(row+1, col+1) = tcomplex(sampleReal, sampleImag);
                }
            }
            // do the fft
            tmp_mat = fft(tmp_mat);
            // transpose the data
            tmp_mat = ~tmp_mat;
            tmp_mat.assign_imag(-tmp_mat.imag());
            // perform the second fft
            tmp_mat = fft(tmp_mat);
            // transpose back
            tmp_mat = ~tmp_mat;
            tmp_mat.assign_imag(-tmp_mat.imag());
            // do fftshift
            tmp_mat = fftshift(tmp_mat);
            // correct phase and add to data mat
            for ( int x = 1; x < data_mat_vec[n].msize()+1; x++ )
                for ( int y = 1; y < data_mat_vec[n].nsize()+1; y++ )
                    data_mat_vec[n](x,y) += phase_mat_vec[coil](x,y)*tmp_mat(x,y);
        }
    }

    // add the data to fid object
    for ( int col = 0; col < m_fid.GetCols(); col++ )
    {
        for ( int row = 0; row < m_fid.GetRows(); row++ )
        {
            cvm::cvector y(nFieldSize);
            for ( size_t n = 0; n < nFieldSize; n++ )
            {
                y(n+1) = data_mat_vec[n](row+1,col+1);
            }
            m_fid.AppendFromVector(y);
        }
    }
}

void tarquin::CFIDReaderGE::LoadFromOptionsSVS(std::string strFilename, const Options& opts, bool WS) 
{
    // attempt to open the file
	std::ifstream file(strFilename.c_str(), std::ios::binary);
    
	if( true == file.fail() )
		throw Exception("failed to open file: '%s'", strFilename.c_str());
    
	std::size_t nOffset      = m_options.nOffset;
	std::size_t nWaterFrames = m_options.nWaterFrames;
	std::size_t nWSFrames    = m_options.nWSFrames;
	std::size_t nFieldSize   = m_options.nFieldSize;
	std::size_t nCoils       = m_options.nCoils;
	std::size_t nFrames      = nWaterFrames + nWSFrames + 1; // 1+ to skip the noise frame
	std::string format       = m_options.nFormat;
	float header_rev         = m_options.nHeaderRev;
	std::size_t nEchoes      = m_options.nEchoes;
    //int voxels               = m_fid.GetRows() * m_fid.GetCols() * m_fid.GetSlices();
    //int ws_avgs              = nWSFrames / voxels;
    //int w_avgs               = nWaterFrames / voxels;

	// flag that water ref data is included in the file	
	if ( nWaterFrames > 0 )
		m_fid.SetCWF(true);
	else
		m_fid.SetCWF(false);
    
    //std::cout << "nOffset      : " << nOffset << std::endl;
    //std::cout << "nWaterFrames : " << nWaterFrames << std::endl;
    //std::cout << "nWSFrames    : " << nWSFrames << std::endl;
    //std::cout << "nFieldSize   : " << nFieldSize << std::endl;
    //std::cout << "nCoils       : " << nCoils << std::endl;
    //std::cout << "nFrames      : " << nFrames << std::endl;
    //std::cout << "Rows         : " << m_fid.GetRows() << std::endl;
    //std::cout << "Cols         : " << m_fid.GetCols() << std::endl;
    //std::cout << "Slices       : " << m_fid.GetSlices() << std::endl;
    //std::cout << "WS avgs      : " << ws_avgs << std::endl;
    //std::cout << "W avgs       : " << w_avgs << std::endl;

	// the FID will be this size
	//y.resize(nFieldSize);

	    
    // vector of first order phase correction values for each coil
    std::vector<tcomplex> phases;
    std::vector<double> amps;
    
    // for Ralphs MEGA-PRESS sequence
    if ( nEchoes == 2 )
        nCoils = nCoils * 2;

    m_fid.ClearVector();
    
	// read data from each frame 
	for( std::size_t nFrame = 1; nFrame <= nFrames-1; nFrame++ ) 
	{
		// advance to the point in the file where this coil's data begins
		//std::size_t nOffsetCoil = nOffset + nCoil*nFrames*nFieldSize*2*sizeof(int);
		//file.seekg(nOffsetCoil, std::ios_base::beg);
            
        cvm::cvector y(nFieldSize);
        y.set(tcomplex(0.0, 0.0));

        cvm::cvector yw(nFieldSize);
        yw.set(tcomplex(0.0, 0.0));

		// for each coil
	    for( std::size_t nCoil = 0; nCoil < nCoils; nCoil++ ) 
		{
            std::size_t nOffsetframe = nOffset + nCoil*nFrames*nFieldSize*2*sizeof(int) + (nFrame-1)*nFieldSize*2*sizeof(int);
		    file.seekg(nOffsetframe, std::ios_base::beg);

			// read in the FID
			for( std::size_t n = 0; n < nFieldSize; n++ ) 
			{
				int sampleReal;
                int sampleImag;
                file.read((char*)&sampleReal, 4);
                file.read((char*)&sampleImag, 4);
                
                // byteswap where necessary
                if ( header_rev < 9 )
                {
                    sampleReal = LongSwap(sampleReal);
                    sampleImag = LongSwap(sampleImag);
                }

                // if this is the first point of the first wus fid find the phase
                // for correction purposes ***fixme should be the average wus
                if ( ( n == 0 ) && ( nFrame <= nWaterFrames ) && ( phases.size() < nCoil + 1 ) && ( nWaterFrames > 0 ) )
                {
                    phases.push_back(exp(-tcomplex(0,1)*arg(tcomplex(sampleReal, sampleImag))));
                    amps.push_back(abs(tcomplex(sampleReal, sampleImag)));

                    // no corrections
                    //phases.push_back(exp(-tcomplex(0,1)));
                    //amps.push_back(1.0);
                    
                    //std::cout << std::endl << phases[nCoil];
                    //std::cout << std::endl << amps[nCoil];
                }
                
                // if no water data available use the WUS data
                if ( ( n == 0 ) && ( nFrame <= nFrames-1 ) && ( phases.size() < nCoil + 1 ) && ( nWaterFrames == 0 ) )
                {
                    phases.push_back(exp(-tcomplex(0,1)*arg(tcomplex(sampleReal, sampleImag))));
                    amps.push_back(abs(tcomplex(sampleReal, sampleImag)));
                    //std::cout << phases[nCoil] << std::endl;
                    //std::cout << amps[nCoil] << std::endl;
                }

                // read in this FID
                // and apply first order phase correction
                if ( nEchoes != 2 )
                {
                    if( nFrame <= nWaterFrames )
                        yw[n+1] += amps[nCoil]*phases[nCoil]*tcomplex(sampleReal, sampleImag);
                    else if( nFrame <= nFrames-1 )
                        y[n+1] += amps[nCoil]*phases[nCoil]*tcomplex(sampleReal, sampleImag);
                }
                else
                {
                    if( nFrame <= nWaterFrames )
                        yw[n+1] += amps[nCoil % nCoils]*phases[nCoil % nCoils]*tcomplex(sampleReal, sampleImag);
                    else if( nFrame <= nFrames-1 )
                        y[n+1] += amps[nCoil % nCoils]*phases[nCoil % nCoils]*tcomplex(sampleReal, sampleImag);
                }
            }
        }

        // normalise by number of coils and sum of amplitudes
        treal amp_sum = 0;
        for ( size_t n = 0; n < amps.size(); n++ )
            amp_sum = amp_sum + amps[n];

        //std::cout << amp_sum << std::endl;

        y /= treal(nCoils)*amp_sum;
        if ( nWaterFrames > 0 )
            yw /= treal(nCoils)*amp_sum;

        if ( WS && nFrame > nWaterFrames )
            m_fid.AppendFromVector(y);

        if ( !WS && nFrame <= nWaterFrames )
            m_fid.AppendFromVector(yw);
    }

    m_fid.SetCols(nWSFrames);

    if ( nWSFrames > 1 )
        m_fid.SetDyn(true);

	// normalise by number of water suppressed and coils
    /*
	y /= treal(nWSFrames*nCoils);
    if ( nWaterFrames > 0 )
        yw /= treal(nWaterFrames*nCoils);
        */
    
    // clear any old data, needed for MEGA-PRESS where data is loaded twice
    //m_fid.ClearVector();
    /*
    if ( WS )
        m_fid.AppendFromVector(y);
    else
        m_fid.AppendFromVector(yw);
        */
}

void tarquin::CFIDReaderGE::LoadFromOptions(std::string strFilename, const Options& opts, bool WS) 
{
    int voxels = m_fid.GetRows() * m_fid.GetCols() * m_fid.GetSlices();
    if ( !boost::filesystem::exists(strFilename) )
        throw Exception("Failed to open %s, has the GE filename been renamed?",strFilename.c_str());

    if ( voxels > 1 )
        LoadFromOptionsCSI(strFilename);
    else
        LoadFromOptionsSVS(strFilename, opts, WS);

}

int tarquin::CFIDReaderGE::LongSwap (int i)
{
	unsigned char b1, b2, b3, b4;

	b1 = i & 255;
	b2 = ( i >> 8 ) & 255;
	b3 = ( i>>16 ) & 255;
	b4 = ( i>>24 ) & 255;

	return ((int)b1 << 24) + ((int)b2 << 16) + ((int)b3 << 8) + b4;
}

short tarquin::CFIDReaderGE::ShortSwap( short s )
{
	unsigned char b1, b2;

	b1 = s & 255;
	b2 = (s >> 8) & 255;

	return (b1 << 8) + b2;
}

float tarquin::CFIDReaderGE::FloatSwap( float f )
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

