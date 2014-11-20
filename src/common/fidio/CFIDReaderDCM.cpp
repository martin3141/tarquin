#include "CFIDReaderDCM.hpp"
#include "CFID.hpp"
#include "CBoswell.hpp"
#include "CDICOMFile.hpp"

#include <iostream>

tarquin::CFIDReaderDCM::CFIDReaderDCM(CFID& fid, CBoswell& log) : 
		CFIDReader(fid, log)
{
}

void tarquin::CFIDReaderDCM::Load(std::string strFilename, const Options& opts, CBoswell& log)
{
    // NOTE tags must be read in the order that they appear in the DICOM
    // file.


    size_t file_size = GetFileSize(strFilename);
    
	if( !file_size )
		throw Exception("file size is zero");

	// attempt to open the file in DICOM format
	CDICOMFile file;

	file.Open(strFilename);

    // find sampling frequency tag
    long bytes = file.MoveToTag("0018", "9052"); 
    if( -1 == bytes )
    {
		throw Exception("Could not find the sampling frequency tag.");
    }
	double fs = 0;
	file.GetFileStream().read((char*)&fs, 8);
    m_fid.SetSamplingFrequency(fs);

    //std::cout << "fs : " << fs << std::endl;
    
    // field strength tag
    bytes = file.MoveToTag("0018", "9098");
    if( -1 == bytes )
    {
		throw Exception("Could not find the field strength tag.");
    }
	double ft = 0;
	file.GetFileStream().read((char*)&ft, 8);
    m_fid.SetTransmitterFrequency(ft*1e6);
    //std::cout << "ft : " << ft << std::endl;

    // Rows
    long bytes_rows = file.MoveToTag("0028", "0010"); 
    if( -1 != bytes_rows )
    {
        uint16_t rows = 0;
        file.GetFileStream().read((char*)&rows, 2);
        m_fid.SetRows(rows);
    }
    else
    {
        long bytes_acqrows = file.MoveToTag("0018", "9095"); 
        if( -1 != bytes_acqrows )
        {
            uint32_t acq_rows = 0;
            file.GetFileStream().read((char*)&acq_rows, 4);
            m_fid.SetRows(acq_rows);
        }
        else
		    throw Exception("Could not find the number of rows");
    }
    
    // Cols
    long bytes_cols = file.MoveToTag("0028", "0011"); 
    if( -1 != bytes_cols )
    {
        uint16_t cols = 0;
        file.GetFileStream().read((char*)&cols, 2);
        m_fid.SetCols(cols);
    }
    else
    {
        long bytes_acqcols = file.MoveToTag("0018", "9234"); 
        if( -1 != bytes_acqcols )
        {
            uint32_t acq_cols = 0;
            file.GetFileStream().read((char*)&acq_cols, 4);
            m_fid.SetCols(acq_cols);
        }
        else
		    throw Exception("Could not find the number of cols");
    }

    // Old code 
    // SpectroscopyAcquisitionPhaseRows 
    /*
    long bytes_acqrows = file.MoveToTag("0018", "9095"); 
    if( -1 != bytes_acqrows )
    {
	    uint32_t acq_rows = 0;
	    file.GetFileStream().read((char*)&acq_rows, 4);
        m_fid.SetRows(acq_rows);
    }
    else
    {
        bytes_acqrows = file.MoveToTag("0028", "0010"); 
        if( -1 == bytes_acqrows )
        {
            throw Exception("Could not find the field strength tag.");
        }
        uint16_t rows = 0;
        file.GetFileStream().read((char*)&rows, 2);
        m_fid.SetRows(rows);
    }

    // SpectroscopyAcquisitionPhaseColumns 
    long bytes_acqcols = file.MoveToTag("0018", "9234"); 
    if ( -1 != bytes_acqcols ) 
    {
        uint32_t acq_cols = 0;
	    file.GetFileStream().read((char*)&acq_cols, 4);
        m_fid.SetCols(acq_cols);
    }
    else
    {
        bytes_acqcols = file.MoveToTag("0028", "0011"); 
        if( -1 == bytes_acqcols )
        {
            throw Exception("Could not find the field strength tag.");
        }
        uint16_t cols = 0;
        file.GetFileStream().read((char*)&cols, 2);
        m_fid.SetCols(cols);
    }
    */
    
    //std::cout << "Cols : " <<  m_fid.GetCols() << std::endl;
    //std::cout << "Rows : " <<  m_fid.GetRows() << std::endl;

    // Find the number of data frames
    long frame_bytes = file.MoveToTag("0028", "0008"); 
    if( -1 == frame_bytes )
    {
		throw Exception("Could not find the number of frames");
    }
    std::vector<char> frames_char(frame_bytes, 0);
	file.GetFileStream().read(&frames_char[0], frame_bytes);
	std::string frames_str(frames_char.begin(), frames_char.end());
    std::istringstream iss(frames_str, std::istringstream::in);
    int frames = 0;
	iss >> frames;

    //std::cout << "Frames : " <<  frames << std::endl;

    // find number of data points
    bytes = file.MoveToTag("0028", "9002"); 
    if( -1 == bytes )
    {
		throw Exception("Could not find the number of data points");
    }
    //file.MoveToTag("0018", "9127"); 
	int N = 0;
	file.GetFileStream().read((char*)&N, 4);

    //std::cout << "N : " << N << std::endl;
    
    // find echo time tag
    bytes = file.MoveToTag("0018", "9082"); 
    if( -1 == bytes )
    {
		throw Exception("Could not find the echo time");
    }
	double echo;
	file.GetFileStream().read((char*)&echo, 8);

    m_fid.SetEchoTime(echo/1000);
    //std::cout << "TE : " << echo/1000 << std::endl;

    // find the acquisition type
    //long str_bytes = file.MoveToTag("0018","9200");
    /*
    file.MoveToTag("0018","9200");
    char type[6];
	file.GetFileStream().read((char*)&type, 6);
    //std::cout << "Acq type bytes : " << str_bytes << std::endl;
    std::cout << "Acq type : " << type << std::endl;
    */
    
    // ImageOrientationPatient
    long ori_bytes = file.MoveToTag("0020","0037");
    if( -1 == ori_bytes )
    {
		throw Exception("Could not find ImageOrientationPatient");
    }
    std::vector<char> ori(ori_bytes, 0);
	file.GetFileStream().read(&ori[0], ori_bytes);
    
	std::string ori_str(ori.begin(), ori.end());
    //std::cout << "Orientation : " << ori_str.substr(0,ori_bytes) << std::endl;
    std::string ori_str_cut = ori_str.substr(0,ori_bytes);
    std::vector<double> ori_vec;
    str2rvec(ori_str_cut, ori_vec);

    std::vector<double> row_ori;
    std::vector<double> col_ori;
    for ( size_t n = 0; n < 6; n++ )
        if ( n < 3 )
            row_ori.push_back(ori_vec[n]);
        else
            col_ori.push_back(ori_vec[n]);

    m_fid.SetRowDirn(row_ori);
    m_fid.SetColDirn(col_ori);
    
    // ImagePositionPatient
    long pos_bytes = file.MoveToTag("0020","0032");
    std::vector<char> pos(pos_bytes, 0);
	file.GetFileStream().read(&pos[0], pos_bytes);
    std::string pos_str(pos.begin(), pos.end());

    //std::cout << "Position : " << pos_str.substr(0,pos_bytes) << std::endl;
    std::string pos_str_cut = pos_str.substr(0,pos_bytes);
    std::vector<double> pos_vec;
    str2rvec(pos_str_cut, pos_vec);
    m_fid.SetPos(pos_vec);

    // SliceThickness
    long thick_bytes = file.MoveToTag("0018","0050");
    
    if( -1 == thick_bytes )
    {
		throw Exception("Could not find the slice thickness");
    }

    std::vector<char> thick(thick_bytes, 0);
	file.GetFileStream().read(&thick[0], thick_bytes);
    std::string thick_str(thick.begin(), thick.end());
    //std::cout << "Thickness : " << thick_str.substr(0,thick_bytes) << std::endl;
    // convert slicethickness to doubles
    std::istringstream istr_slice(thick_str.substr(0,thick_bytes));
    double vox_slice_dim = 0.0;
    istr_slice >> vox_slice_dim;

    // PixelSpacing
    long dims_bytes = file.MoveToTag("0028","0030");
    std::vector<char> dims(dims_bytes, 0);
	file.GetFileStream().read(&dims[0], dims_bytes);
    std::string dims_str(dims.begin(), dims.end());
    //std::cout << "Pixel Spacing : " << dims_str.substr(0,dims_bytes) << std::endl;
    // convert pixel spacing to doubles
    std::string str = dims_str.substr(0,dims_bytes);
    std::vector<double> row_col;
    str2rvec(str, row_col);
    
    // save voxel dims to fid object
    std::vector<double> voxel_dim;
    voxel_dim.resize(3);
    voxel_dim[0] = row_col[0];
    voxel_dim[1] = row_col[1];
    voxel_dim[2] = vox_slice_dim;
    m_fid.SetVoxelDim(voxel_dim);

    // find the manufacturer
    long manu_bytes = file.MoveToTag("0008","0070");
    std::vector<char> manu(manu_bytes, 0);
	file.GetFileStream().read(&manu[0], manu_bytes);
	std::string manu_str(manu.begin(), manu.end());

    //std::cout << "Manufacturer : " << manu_str << "M" << std::endl;

    // localisation paras (don't seem to work well for Philips data)
    long local_bytes = file.MoveToTag("0018", "9104"); 
    double slice_thick = 0;
    if ( -1 != local_bytes )
        file.GetFileStream().read((char*)&slice_thick, 8);

    //std::cout << "Slice thick : " << slice_thick << std::endl;

    //std::cout << "Bytes : " << local_bytes << std::endl;
    if ( ( -1 != local_bytes ) && ( manu_str.compare("SIEMENS") == 0 || manu_str.compare("SIEMENS ") == 0 ) )
    {
        file.MoveToTag("0018", "9104", false); 
        double row_thick;
        file.GetFileStream().read((char*)&row_thick, 8);
        //std::cout << "Row thick : " << row_thick << std::endl;

        file.MoveToTag("0018", "9104", false); 
        double col_thick;
        file.GetFileStream().read((char*)&col_thick, 8);
        //std::cout << "Col thick : " << col_thick << std::endl;

        std::vector<double> voi_dim;
        voi_dim.resize(3);
        voi_dim[0] = row_thick;
        voi_dim[1] = col_thick;
        voi_dim[2] = slice_thick;

        m_fid.SetVoiDim(voi_dim);
        
        m_fid.SetSlices(frames);
    }
    
    // this is required, probally because the reader gets past the eof if one of the tags are missing
	file.Close();
	file.Open(strFilename);
    
    // Are we looking at dynamic data? need to check TemporalPositionIndex tags (0020,9128)
    // if > 1 then this looks like dynamic data
    bool dyn_study = false;
    uint32_t dyn_num = 0;
	long nBytesDyn = file.MoveToTag("0020", "9128"); 
    if ( nBytesDyn != -1 )
    {
        uint32_t first = 0;
        file.GetFileStream().read((char*)&first, 4);
        //std::cout << first << std::endl;
        nBytesDyn = file.MoveToTag("0020", "9128", false); 
        if ( nBytesDyn != -1 )
        {
            file.GetFileStream().read((char*)&dyn_num, 4);
            //std::cout << first << std::endl;
            // looks like a dyanmic Philips study
            if ( dyn_num > 1 )
                dyn_study = true;

            // keep munching tags until we get to the end
            while ( nBytesDyn != -1 )
            {
                nBytesDyn = file.MoveToTag("0020", "9128", false); 
                if ( nBytesDyn != -1 )
                    file.GetFileStream().read((char*)&dyn_num, 4);
            }
            //std::cout << std::endl << dyn_num << " dynamic scans available." << std::endl;
        }
    }

    if ( frames > 1 && m_fid.GetRows()*m_fid.GetCols() == 1 && dyn_study )
    {
        m_fid.SetCols(dyn_num);
        m_fid.SetDyn(true);
        //m_fid.SetCols(frames);
    }
    
    // this is required, probally because the reader gets past the eof if one of the tags are missing
    file.Close();
	file.Open(strFilename);
    
	long nBytesInFID = file.MoveToTag("5600", "0020"); 
	if( -1 == nBytesInFID ) 
	{
		throw Exception("file did not contain the magic tag (5600, 0020) in WS data.");
	}
    
    //std::cout << std::endl << "FID data is " << nBytesInFID << " bytes" << std::endl;
    //std::cout << std::endl << "Data should be " << m_fid.GetRows()*m_fid.GetCols()*m_fid.GetSlices()*N*8  << " bytes" << std::endl;

    if ( ( nBytesInFID % (m_fid.GetRows()*m_fid.GetCols()*m_fid.GetSlices()*N*8)) != 0 )
    {
        std::ostringstream error_stream;
        error_stream << "file does not contain a consistent number of data points" << std::endl;
        error_stream << "FID data is " << nBytesInFID << " bytes" << std::endl;
        error_stream << "Data should be a multiple of " << m_fid.GetSlices()*m_fid.GetRows()*m_fid.GetCols()*N*8  << " bytes" << std::endl;
        error_stream << "N    = " << N << std::endl;
        error_stream << "Rows = " << m_fid.GetRows() << std::endl;
        error_stream << "Cols = " << m_fid.GetCols() << std::endl;
        error_stream << "Cols = " << m_fid.GetSlices() << std::endl;
        error_stream << nBytesInFID/(8.0*N) << " voxels in file" << std::endl;
		throw Exception(error_stream.str().c_str());
    }

    long avgs = nBytesInFID/(m_fid.GetRows()*m_fid.GetCols()*m_fid.GetSlices()*N*8);
    //std::cout << std::endl << "This file contains " << avgs << " averages" << std::endl;

    // if fid is double the expected and SVS then there is a
    // good chance the WUS is in there too
    //if ( avgs == 2 && m_fid.GetRows()*m_fid.GetCols() == 1 )
    if ( avgs == 2 )
    {
        //std::cout << std::endl << "This file contains W data also" << std::endl;
        m_fid.SetCWF(true);
    }

	ReadFIDData(file.GetFileStream(), m_fid.GetSlices()*m_fid.GetRows()*m_fid.GetCols()*N*8, 0, 1);
	
    // swap data organisation for Siemens data
    //std::cout << manu_str << std::endl;

    if ( manu_str.compare("SIEMENS") == 0 || manu_str.compare("SIEMENS ") == 0 )
    {
        m_fid.swap_row_col();
    }
    else if ( manu_str.compare("Philips Medical Systems") == 0 || manu_str.compare("Philips Medical Systems ") == 0 )
    {
        
        // Philips don't really get their definitions of IPP and IOP correct 
        // so you need to delve into the private togs if the VOI is non
        // central or rotated WRT the FOV...grim

        bytes = file.MoveToTag("2005", "1055"); // rotation angle
        float angle = 0;
        if( -1 != bytes ) // otherwise just assume no angle
            file.GetFileStream().read((char*)&angle, 4);

        bytes = file.MoveToTag("2005", "1071"); // row angle 
        float row_ang = 0;
        if( -1 != bytes )
            file.GetFileStream().read((char*)&row_ang, 4);

        bytes = file.MoveToTag("2005", "1072"); // col angle
        float col_ang = 0;
        if( -1 != bytes )
            file.GetFileStream().read((char*)&col_ang, 4);

        bytes = file.MoveToTag("2005", "1073"); // slice angle
        float slice_ang = 0;
        if( -1 != bytes )
            file.GetFileStream().read((char*)&slice_ang, 4);

        bytes = file.MoveToTag("2005", "107A"); // row offset
        float row_off = pos_vec[0];
        if( -1 != bytes ) // otherwise just assume no offset 
            file.GetFileStream().read((char*)&row_off, 4);

        bytes = file.MoveToTag("2005", "1078"); // col offset
        float col_off = pos_vec[1];
        if( -1 != bytes ) // otherwise just assume no offset 
            file.GetFileStream().read((char*)&col_off, 4);

        bytes = file.MoveToTag("2005", "1079"); // slice offset
        float slice_off = pos_vec[2];
        if( -1 != bytes ) // otherwise just assume no offset 
            file.GetFileStream().read((char*)&slice_off, 4);

        //std::cout << "Slice off : " << slice_off << std::endl;
        //m_fid.SetPos(pos_vec);

        cvm::rvector cvm_row_ori(3);
        cvm_row_ori(1) = row_ori[0];
        cvm_row_ori(2) = row_ori[1];
        cvm_row_ori(3) = row_ori[2];

        cvm::rvector cvm_col_ori(3);
        cvm_col_ori(1) = col_ori[0];
        cvm_col_ori(2) = col_ori[1];
        cvm_col_ori(3) = col_ori[2];

        cvm::rvector cvm_true_row(3);
        cvm_true_row(1) = 1;
        cvm_true_row(2) = 0;
        cvm_true_row(3) = 0;
        
        cvm::rvector cvm_true_col(3);
        cvm_true_col(1) = 0;
        cvm_true_col(2) = 1;
        cvm_true_col(3) = 0;

        cvm::rvector cvm_true_slice(3);
        cvm_true_slice(1) = 0;
        cvm_true_slice(2) = 0;
        cvm_true_slice(3) = 1;

        /*double rotate = acos(cvm_row_ori*cvm_true_row)*180.0/M_PI;

        std::cout << std::endl;
        std::cout << "Row         : " << cvm_row_ori << std::endl;
        std::cout << "Col         : " << cvm_col_ori << std::endl;
        std::cout << "Row angle   : " << row_ang << std::endl;
        std::cout << "Col angle   : " << col_ang << std::endl;
        std::cout << "Slice angle : " << slice_ang << std::endl;
        std::cout << "Angle       : " << angle << std::endl;
        std::cout << "Rotation    : " << rotate << std::endl;

        angle = 0;
        rotate = 0;
        */
        
        /*
        std::cout << std::endl;
        std::cout << "Row angle   : " << row_ang << std::endl;
        std::cout << "Col angle   : " << col_ang << std::endl;
        std::cout << "Slice angle : " << slice_ang << std::endl;
        std::cout << "Row offset   : " << row_off << std::endl;
        std::cout << "Col offset   : " << col_off << std::endl;
        std::cout << "Slice offset : " << slice_off << std::endl;
        */


        cvm::rvector new_cvm_row_ori(3);
        rotate_vec(cvm_true_row, cvm_true_slice, (col_ang)*M_PI/180.0, new_cvm_row_ori);
        rotate_vec(new_cvm_row_ori, cvm_true_col, (row_ang)*M_PI/180.0, new_cvm_row_ori);
        rotate_vec(new_cvm_row_ori, cvm_true_row, (slice_ang)*M_PI/180.0, new_cvm_row_ori);
        //std::cout << new_cvm_row_ori << std::endl;

        cvm::rvector new_cvm_col_ori(3);
        rotate_vec(cvm_true_col, cvm_true_slice, (col_ang)*M_PI/180.0, new_cvm_col_ori);
        rotate_vec(new_cvm_col_ori, cvm_true_col, (row_ang)*M_PI/180.0, new_cvm_col_ori);
        rotate_vec(new_cvm_col_ori, cvm_true_row, (slice_ang)*M_PI/180.0, new_cvm_col_ori);
        //std::cout << new_cvm_col_ori << std::endl;

        //cvm::rvector cvm_slice_ori(3);
        //m_fid.cross(cvm_row_ori, cvm_col_ori, cvm_slice_ori);

        // rotate row and col by angle
        //cvm::rvector new_cvm_row_ori(3);
        //rotate_vec(cvm_row_ori, cvm_slice_ori, (-angle)*M_PI/180.0, new_cvm_row_ori);
        cvm_row_ori = new_cvm_row_ori; 
        row_ori[0] = cvm_row_ori(1);
        row_ori[1] = cvm_row_ori(2);
        row_ori[2] = cvm_row_ori(3);
        m_fid.SetRowDirn(row_ori);

        //cvm::rvector new_cvm_col_ori(3);
        //rotate_vec(cvm_col_ori, cvm_slice_ori, (-angle)*M_PI/180.0, new_cvm_col_ori);
        cvm_col_ori = new_cvm_col_ori; 
        col_ori[0] = cvm_col_ori(1);
        col_ori[1] = cvm_col_ori(2);
        col_ori[2] = cvm_col_ori(3);
        m_fid.SetColDirn(col_ori);


        cvm::rvector cvm_pos_vec(3);
        /*cvm_pos_vec(1) = pos_vec[0];
        cvm_pos_vec(2) = pos_vec[1];
        cvm_pos_vec(3) = pos_vec[2];*/
        cvm_pos_vec(1) = row_off;
        cvm_pos_vec(2) = col_off;
        cvm_pos_vec(3) = slice_off;

        const std::vector<double>& voxel_dim = m_fid.GetVoxelDim();

        cvm_pos_vec += -cvm_col_ori*voxel_dim[0]*(0.5*(m_fid.GetRows()-1)) -cvm_row_ori*voxel_dim[1]*(0.5*(m_fid.GetCols()-1));

        pos_vec[0] = cvm_pos_vec(1);
        pos_vec[1] = cvm_pos_vec(2);
        pos_vec[2] = cvm_pos_vec(3);

        m_fid.SetPos(pos_vec);
    }
}

void tarquin::CFIDReaderDCM::LoadW(std::string strFilename, const Options& opts)
{
    size_t file_size = GetFileSize(strFilename);

	if( !file_size )
		throw Exception("file size is zero");

	// attempt to open the file in DICOM format
	CDICOMFile file;

	file.Open(strFilename);
    
    // Rows
    long bytes_rows = file.MoveToTag("0028", "0010"); 
    if( -1 != bytes_rows )
    {
        uint16_t rows = 0;
        file.GetFileStream().read((char*)&rows, 2);
        m_fid.SetRows(rows);
    }
    else
    {
        long bytes_acqrows = file.MoveToTag("0018", "9095"); 
        if( -1 != bytes_acqrows )
        {
            uint32_t acq_rows = 0;
            file.GetFileStream().read((char*)&acq_rows, 4);
            m_fid.SetRows(acq_rows);
        }
        else
		    throw Exception("Could not find the number of rows");
    }
    
    // Cols
    long bytes_cols = file.MoveToTag("0028", "0011"); 
    if( -1 != bytes_cols )
    {
        uint16_t cols = 0;
        file.GetFileStream().read((char*)&cols, 2);
        m_fid.SetCols(cols);
    }
    else
    {
        long bytes_acqcols = file.MoveToTag("0018", "9234"); 
        if( -1 != bytes_acqcols )
        {
            uint32_t acq_cols = 0;
            file.GetFileStream().read((char*)&acq_cols, 4);
            m_fid.SetCols(acq_cols);
        }
        else
		    throw Exception("Could not find the number of cols");
    }

    // Old code
    // SpectroscopyAcquisitionPhaseRows 
    /*
    long bytes_acqrows = file.MoveToTag("0018", "9095"); 
    if( -1 != bytes_acqrows )
    {
	    uint32_t acq_rows = 0;
	    file.GetFileStream().read((char*)&acq_rows, 4);
        m_fid.SetRows(acq_rows);
    }
    else
    {
        file.MoveToTag("0028", "0010"); 
        uint16_t rows = 0;
        file.GetFileStream().read((char*)&rows, 2);
        m_fid.SetRows(rows);
    }

    // SpectroscopyAcquisitionPhaseColumns 
    long bytes_acqcols = file.MoveToTag("0018", "9234"); 
    if ( -1 != bytes_acqcols ) 
    {
        uint32_t acq_cols = 0;
	    file.GetFileStream().read((char*)&acq_cols, 4);
        m_fid.SetCols(acq_cols);
    }
    else
    {
        file.MoveToTag("0028", "0011"); 
        uint16_t cols = 0;
        file.GetFileStream().read((char*)&cols, 2);
        m_fid.SetCols(cols);
    }*/
    
    // find number of data points
    file.MoveToTag("0018", "9127"); 
	int N;
	file.GetFileStream().read((char*)&N, 4);
	    
    // Find the number of data frames
    long frame_bytes = file.MoveToTag("0028", "0008"); 
    if( -1 == frame_bytes )
    {
		throw Exception("Could not find the number of frames");
    }
    std::vector<char> frames_char(frame_bytes, 0);
	file.GetFileStream().read(&frames_char[0], frame_bytes);
	std::string frames_str(frames_char.begin(), frames_char.end());
    std::istringstream iss(frames_str, std::istringstream::in);
    int frames = 0;
	iss >> frames;

    // Are we looking at dynamic data? need to check TemporalPositionIndex tags (0020,9128)
    // if > 1 then this looks like dynamic data
    bool dyn_study = false;
    uint32_t dyn_num = 0;
	long nBytesDyn = file.MoveToTag("0020", "9128"); 
    if ( nBytesDyn != -1 )
    {
        uint32_t first = 0;
        file.GetFileStream().read((char*)&first, 4);
        //std::cout << first << std::endl;
        nBytesDyn = file.MoveToTag("0020", "9128", false); 
        if ( nBytesDyn != -1 )
        {
            file.GetFileStream().read((char*)&dyn_num, 4);
            //std::cout << first << std::endl;
            // looks like a dyanmic Philips study
            if ( dyn_num > 1 )
                dyn_study = true;

            // keep munching tags until we get to the end
            while ( nBytesDyn != -1 )
            {
                nBytesDyn = file.MoveToTag("0020", "9128", false); 
                if ( nBytesDyn != -1 )
                    file.GetFileStream().read((char*)&dyn_num, 4);
            }
            //std::cout << std::endl << dyn_num << " dynamic scans available." << std::endl;
        }
    }

    /*std::cout << std::endl << m_fid.GetCols() << std::endl;
    std::cout << m_fid.GetRows() << std::endl;
    std::cout << frames << std::endl;
    std::cout << dyn_num << std::endl;
    std::cout << dyn_study << std::endl;*/

    if ( frames > 1 && m_fid.GetRows()*m_fid.GetCols() == 1 && dyn_study )
    {
        //m_fid.SetCols(frames);
        m_fid.SetCols(dyn_num);
        m_fid.SetDyn(true);
        //std::cout << m_fid.GetCols() << std::endl;
        //std::cout << m_fid.GetRows() << std::endl;
    }
    
    // next two lines are a bodge, not sure why we need them, something to do with
    // reading past eof in above while loop...

	file.Close();
	file.Open(strFilename);

    long nBytesInFID = file.MoveToTag("5600", "0020"); 
	if( -1 == nBytesInFID ) 
	{
		throw Exception("file did not contain the magic tag (5600, 0020) in WUS data.");
	}

    ReadFIDData(file.GetFileStream(), m_fid.GetRows()*m_fid.GetCols()*N*16, m_fid.GetRows()*m_fid.GetCols()*N*8, 1);
}


void tarquin::CFIDReaderDCM::ReadFIDData(std::ifstream& file, std::size_t nLength, std::size_t byte_offset, std::size_t averages)
{
    //std::cout << std::endl << "nLength : " << nLength << std::endl;
    //std::cout << "byte offset : " << byte_offset << std::endl;

    // skip to the part of the file we want
    file.seekg(byte_offset, std::ios::cur);

	// floats are in IEEE 32 bit, little endian
	// (Siemens scanners use glorified PCs)
	float real_part_ieee;
	float imag_part_ieee;
    
	long nSamples = (nLength - byte_offset) / 8 / m_fid.GetRows() / m_fid.GetCols() / m_fid.GetSlices();
    std::cout << "samples : " << nSamples << std::endl;
    for ( int avs = 0; avs < averages; avs++ )
    {
        for ( int fid = 0; fid < ( m_fid.GetRows() * m_fid.GetCols() * m_fid.GetSlices() ); fid++ ) 
        {
            //std::cout << fid << std::endl;
            int n = 0;
            cvm::cvector FID(nSamples);
            for( std::size_t nBytesRead = 0; nBytesRead < (nLength - byte_offset) / m_fid.GetRows() / m_fid.GetCols() / m_fid.GetSlices(); nBytesRead+= 8 ) 
            {
                file.read((char*)&real_part_ieee, 4);
                file.read((char*)&imag_part_ieee, 4);

                if( true == file.eof() ) 
                    throw Exception("ran out file before reading in target number of bytes");

                // only for big endian boxes (Sun?)
#ifdef __BIG_ENDIAN__
                BYTESWAP(real_part_ieee);
                BYTESWAP(imag_part_ieee);
#endif
                // store this sample
                FID[n+1] =  tcomplex(real_part_ieee, -imag_part_ieee);
                n++;
            }
            //std::cout << FID.size() << std::endl;
            //plot(FID);
            if ( avs == 0 )
            {
                m_fid.AppendFromVector(FID);
            }
            else
            {
                cvm::cvector& current_fid = m_fid.GetVectorFID(fid);
                current_fid += FID;
            }
        }
    }
    //std::cout << m_fid.GetVectorFID().size() << " fids saved" << std::endl;
}

template <class T>
bool tarquin::CFIDReaderDCM::from_string(T& t, const std::string& s, std::ios_base& (*f)(std::ios_base&))
{
	std::istringstream iss(s);
	return !(iss >> f >> t).fail();
}

void tarquin::CFIDReaderDCM::rotate_vec(const cvm::rvector &vec_in, const cvm::rvector &ax, double theta, cvm::rvector& vec_out)
{
    double ct = cos(theta);
    double st = sin(theta);
    cvm::rmatrix rotate_mat(3,3);
    rotate_mat(1,1) = ct + pow(ax(1),2)*(1-ct);
    rotate_mat(1,2) = ax(1)*ax(2)*(1-ct)-ax(3)*st;
    rotate_mat(1,3) = ax(1)*ax(3)*(1-ct)+ax(2)*st;
    rotate_mat(2,1) = ax(2)*ax(1)*(1-ct)+ax(3)*st;
    rotate_mat(2,2) = ct+pow(ax(2),2)*(1-ct);
    rotate_mat(2,3) = ax(2)*ax(3)*(1-ct)-ax(1)*st;
    rotate_mat(3,1) = ax(3)*ax(1)*(1-ct)-ax(2)*st;
    rotate_mat(3,2) = ax(3)*ax(2)*(1-ct)+ax(1)*st;
    rotate_mat(3,3) = ct+pow(ax(3),2)*(1-ct);

    vec_out.resize(3);
    vec_out = rotate_mat*vec_in;

    // normalise at the end
    vec_out = vec_out/vec_out.norm2();
}

