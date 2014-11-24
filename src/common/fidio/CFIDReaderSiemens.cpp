#include "CFIDReaderSiemens.hpp"
#include "CFID.hpp"
#include "CBoswell.hpp"
#include "CDICOMFile.hpp"

#include <iostream>

tarquin::CFIDReaderSiemens::CFIDReaderSiemens(CFID& fid, CBoswell& log) : 
		CFIDReader(fid, log)
{
}

void tarquin::CFIDReaderSiemens::Load(std::string strFilename, const Options& opts, CBoswell& log)
{
	size_t file_size = GetFileSize(strFilename);

	if( !file_size )
		throw Exception("file size is zero");

	// read in the whole file
	std::vector<unsigned char> buffer(file_size);
	{
		std::ifstream fin( strFilename.c_str(), std::ios_base::binary);

		if( !fin.is_open() )
			throw Exception("failed to open file");

		fin.read(reinterpret_cast<char*>(&buffer[0]), file_size);
	}

    // some variables needed for CSI data
    double rotation = 0;
    double col_dim = 0;
    double row_dim = 0;
    double slice_dim = 0;
    double norm_tra = 0;
    double norm_cor = 0;
    double norm_sag = 0;
    double pos_tra = 0;
    double pos_cor = 0;
    double pos_sag = 0;
    double voi_thick = 0;
    double voi_col = 0;
    double voi_row = 0;

	// go through each 'line' of the file, a line is determined by the point
	// where we find a '\n'
    bool ascconv_start = false;	
	for( size_t n = 0; n < buffer.size(); ++n )
	{
		std::string line;

		// read until a newline
		for( ; n < buffer.size(); ++n )
			if( buffer[n] == 0x0A )
				break;
			else 
				line += buffer[n];

		// blam! we have found a 'line'

		if ( line.substr(0,9) == "ulVersion" )
        {
            size_t found = line.find("=");
            //std::cout << line << std::endl;
            //std::cout << found << std::endl;
            // the following is a sanity check, sometimes ulVersion appears in the binary header
            // for some reason
            if ( found < 50 )
                ascconv_start = true;
        }
        
        /*if ( line.size() > 21 )
            if ( line.substr(line.size() - 21 ,line.size()) == "### ASCCONV BEGIN ###" )
                ascconv_start = true;
                */
        
        if ( line.substr(0,19) == "### ASCCONV END ###"  && ascconv_start )
            break;

        if ( ascconv_start )
        {
            //std::cout << line << std::endl;
            if ( line.substr(0,22) == "sRXSPEC.alDwellTime[0]" )
            {
                // find the = char 
                size_t found = line.find("=");
                float fs = 0;
                from_string<float>(fs, line.substr(found+2,line.size()), std::dec);
                m_fid.SetSamplingFrequency(1e9/(fs*2));
            } 
            else if ( line.substr(0,35) == "sTXSPEC.asNucleusInfo[0].lFrequency" )
            {
                // find the = char 
                size_t found = line.find("=");
                float ft = 0;
                from_string<float>(ft, line.substr(found+2,line.size()), std::dec);
                m_fid.SetTransmitterFrequency(ft);
            }
            else if ( line.substr(0,7) == "alTE[0]" )
            {
                // find the = char 
                size_t found = line.find("=");
                float te;
                from_string<float>(te, line.substr(found+2,line.size()), std::dec);
                m_fid.SetEchoTime(te/1e6);
            }
            else if ( line.substr(0,31) == "sSpecPara.lFinalMatrixSizePhase" )
            {
                // find the = char 
                size_t found = line.find("=");
                int cols;
                from_string<int>(cols, line.substr(found+2,line.size()), std::dec);
                m_fid.SetCols(cols);
            }
            else if ( line.substr(0,30) == "sSpecPara.lFinalMatrixSizeRead" )
            {
                // find the = char 
                size_t found = line.find("=");
                int rows;
                from_string<int>(rows, line.substr(found+2,line.size()), std::dec);
                m_fid.SetRows(rows);
            }
            else if ( line.substr(0,31) == "sSpecPara.lFinalMatrixSizeSlice" )
            {
                // find the = char 
                size_t found = line.find("=");
                int slices;
                from_string<int>(slices, line.substr(found+2,line.size()), std::dec);
                m_fid.SetSlices(slices);
            }

            else if ( line.substr(0,29) == "sSpecPara.sVoI.sPosition.dSag" )
            {
                // find the = char 
                size_t found = line.find("=");
                from_string<double>(pos_sag, line.substr(found+2,line.size()), std::dec);
                //std::cout << "Sag : " << pos_sag << std::flush << std::endl;
            }

            else if ( line.substr(0,29) == "sSpecPara.sVoI.sPosition.dCor" )
            {
                // find the = char 
                size_t found = line.find("=");
                from_string<double>(pos_cor, line.substr(found+2,line.size()), std::dec);
                //std::cout << "Cor : " << pos_cor << std::flush <<  std::endl;
            }

            else if ( line.substr(0,29) == "sSpecPara.sVoI.sPosition.dTra" )
            {
                // find the = char 
                size_t found = line.find("=");
                from_string<double>(pos_tra, line.substr(found+2,line.size()), std::dec);
                //std::cout << "Tra : " << pos_tra << std::flush << std::endl;
            }

            else if ( line.substr(0,27) == "sSpecPara.sVoI.sNormal.dSag" )
            {
                // find the = char 
                size_t found = line.find("=");
                from_string<double>(norm_sag, line.substr(found+2,line.size()), std::dec);
                //std::cout << "Sag : " << norm_sag << std::endl;
            }

            else if ( line.substr(0,27) == "sSpecPara.sVoI.sNormal.dCor" )
            {
                // find the = char 
                size_t found = line.find("=");
                from_string<double>(norm_cor, line.substr(found+2,line.size()), std::dec);
                //std::cout << "Cor : " << norm_cor << std::flush << std::endl;
            }

            else if ( line.substr(0,27) == "sSpecPara.sVoI.sNormal.dTra" )
            {
                // find the = char 
                size_t found = line.find("=");
                from_string<double>(norm_tra, line.substr(found+2,line.size()), std::dec);
                //std::cout << "Tra : " << norm_tra << std::endl;
            }

            else if ( line.substr(0,33) == "sSliceArray.asSlice[0].dThickness" )
            {
                // find the = char 
                size_t found = line.find("=");
                from_string<double>(slice_dim, line.substr(found+2,line.size()), std::dec);
                //std::cout << "Thick : " << slice_dim << std::endl;
            }

            else if ( line.substr(0,32) == "sSliceArray.asSlice[0].dPhaseFOV" )
            {
                // find the = char 
                size_t found = line.find("=");
                from_string<double>(row_dim, line.substr(found+2,line.size()), std::dec);
                //std::cout << "Row : " << row_dim << std::endl;
            }

            else if ( line.substr(0,34) == "sSliceArray.asSlice[0].dReadoutFOV" )
            {
                // find the = char 
                size_t found = line.find("=");
                from_string<double>(col_dim, line.substr(found+2,line.size()), std::dec);
                //std::cout << "Col : " << col_dim << std::endl;
            }

            else if ( line.substr(0,34) == "sSliceArray.asSlice[0].dInPlaneRot" )
            {
                // find the = char 
                size_t found = line.find("=");
                from_string<double>(rotation, line.substr(found+2,line.size()), std::dec);
                //std::cout << "Rotation : " << rotation << std::endl;
            }

            else if ( line.substr(0,25) == "sSpecPara.sVoI.dThickness" )
            {
                // find the = char 
                size_t found = line.find("=");
                from_string<double>(voi_thick, line.substr(found+2,line.size()), std::dec);
                //std::cout << "VOI thickness : " << voi_thick << std::endl;
            }

            else if ( line.substr(0,24) == "sSpecPara.sVoI.dPhaseFOV" )
            {
                // find the = char 
                size_t found = line.find("=");
                from_string<double>(voi_row, line.substr(found+2,line.size()), std::dec);
                //std::cout << "VOI row thickness : " << voi_col << std::endl;
            }

            else if ( line.substr(0,26) == "sSpecPara.sVoI.dReadoutFOV" )
            {
                // find the = char 
                size_t found = line.find("=");
                from_string<double>(voi_col, line.substr(found+2,line.size()), std::dec);
                //std::cout << "VOI col thickness : " << voi_row << std::endl;
            }
        }
        
    }

    //if ( m_fid.GetVoxelCount() > 1 )
        {
            std::vector<double> voxel_dim;
            voxel_dim.resize(3);
            if ( m_fid.GetVoxelCount() > 1 )
            {
                voxel_dim[0] = row_dim/m_fid.GetRows();
                voxel_dim[1] = col_dim/m_fid.GetCols();
                voxel_dim[2] = slice_dim;
            }
            else
            {
                voxel_dim[0] = voi_row;
                voxel_dim[1] = voi_col;
                voxel_dim[2] = voi_thick;
            }
            m_fid.SetVoxelDim(voxel_dim);

            std::vector<double> voi_dim;
            voi_dim.resize(3);
            voi_dim[0] = voi_row;
            voi_dim[1] = voi_col;
            voi_dim[2] = voi_thick;
            m_fid.SetVoiDim(voi_dim);

            cvm::rvector ima_norm(3);
            ima_norm(1) = norm_sag;
            ima_norm(2) = norm_cor;
            ima_norm(3) = norm_tra;

            // convert the plane normal cosine and in-plane rotation to
            // DICOM style row and column direction cosines 
            cvm::rvector new_x;
            const double x[] = {1., 0., 0.};
            cvm::rvector x_dirn(3);
            x_dirn.assign(x);
            cvm::rvector x_new;
            rotate_vec( x_dirn, ima_norm, -rotation, x_new);
            cvm::rvector new_col;
            m_fid.cross(ima_norm, x_new, new_col);
            cvm::rvector new_row;
            m_fid.cross(new_col, ima_norm, new_row);
            
            std::vector<double> new_row_vec;
            new_row_vec.resize(3);
            new_row_vec[0] = new_row(1);
            new_row_vec[1] = new_row(2);
            new_row_vec[2] = new_row(3);

            std::vector<double> new_col_vec;
            new_col_vec.resize(3);
            new_col_vec[0] = new_col(1);
            new_col_vec[1] = new_col(2);
            new_col_vec[2] = new_col(3);

            m_fid.SetRowDirn(new_row_vec);
            m_fid.SetColDirn(new_col_vec);

            //std::cout << std::endl;
            //std::cout << "Row dirn : " << new_row << std::endl;
            //std::cout << "Col dirn : " << new_col << std::endl;

            // convert CSI grid ceter position to 
            // DICOM style row[0], col[0] position
            cvm::rvector ima_pos(3);
            ima_pos(1) = pos_sag;
            ima_pos(2) = pos_cor;
            ima_pos(3) = pos_tra;
            
            cvm::rvector dcm_pos(3);
            dcm_pos = ima_pos - new_row * ( m_fid.GetRows() / 2.0 - 0.5 ) * row_dim/m_fid.GetRows()
                    - new_col * ( m_fid.GetCols() / 2.0 - 0.5 ) * col_dim/m_fid.GetCols();

            std::vector<double> dcm_pos_vec;
            dcm_pos_vec.resize(3);
            dcm_pos_vec[0] = dcm_pos(1);
            dcm_pos_vec[1] = dcm_pos(2);
            dcm_pos_vec[2] = dcm_pos(3);
            m_fid.SetPos(dcm_pos_vec);
        }


	// attempt to open the file in DICOM format
	CDICOMFile file;

	file.Open(strFilename);

	// move to point where FID starts
	long nBytesInFID = file.MoveToTag("7FE1", "1010"); 

	if( -1 == nBytesInFID ) 
	{
		throw Exception("file did not contain the magic tag (7FE1, 1010)");
	}

	ReadFIDData(file.GetFileStream(), nBytesInFID);

    m_fid.swap_row_col();
}

void tarquin::CFIDReaderSiemens::ReadFIDData(std::ifstream& file, std::size_t nLength)
{
	// floats are in IEEE 32 bit, little endian
	// (Siemens scanners use glorified PCs)
	float real_part_ieee;
	float imag_part_ieee;

	int n_fids = m_fid.GetRows() * m_fid.GetCols() * m_fid.GetSlices();

	long nSamples = nLength / 8 / n_fids ;

	for ( int fids = 0; fids < n_fids; fids++ )
	{
		cvm::cvector FID(nSamples);
		int n = 0;
		for( std::size_t nBytesRead = 0; nBytesRead < nLength / n_fids; nBytesRead+= 8 ) 
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
			FID[n+1] =  tcomplex(real_part_ieee, imag_part_ieee);
			n++;
		}
		m_fid.AppendFromVector(FID);
	}
}

template <class T>
bool tarquin::CFIDReaderSiemens::from_string(T& t, const std::string& s, std::ios_base& (*f)(std::ios_base&))
{
	std::istringstream iss(s);
	return !(iss >> f >> t).fail();
}

void tarquin::CFIDReaderSiemens::rotate_vec(const cvm::rvector &vec_in, const cvm::rvector &ax, double theta, cvm::rvector& vec_out)
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


