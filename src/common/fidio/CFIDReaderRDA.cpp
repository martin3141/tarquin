#include "CFIDReaderRDA.hpp"
#include "CBoswell.hpp"
#include "CFID.hpp"
#include <istream>
#include <boost/tokenizer.hpp>

using namespace boost;


tarquin::CFIDReaderRDA::CFIDReaderRDA(CFID& fid, CBoswell& log) : 
	CFIDReader(fid, log)
{
}

void tarquin::CFIDReaderRDA::Load(std::string strFilename, const Options& opts, CBoswell& log)
{
	// tokenize the file and get the position of the binary data
	std::streampos binpos = LexRDA(strFilename);

	if( -1 == binpos ) 
		throw Exception("file does not look like correct format");

	// look for the "VectorSize" token
	std::istringstream ins;
	long nSamples = -1;
	for( TokenList::iterator it = m_tokens.begin(); it != m_tokens.end(); ++it ) 
	{
		// clear up for this conversion
		ins.clear();

		if( it->first == "VectorSize" ) 
		{
			it++;
			ins.str(it->first);
			ins >> nSamples;
		}

		if( it->first == "DwellTime" ) 
		{
			it++;
			ins.str(it->first);
			treal dwell_time;
			ins >> dwell_time;

			// it appears dwell time is in microseconds (obviously...)
			treal fs = 1000000 / dwell_time;

			m_fid.SetSamplingFrequency(fs);
		} 

		if( it->first == "MRFrequency" ) 
		{
			it++;
			ins.str(it->first);
			treal ft;
			ins >> ft;

			// it appears transmitter frequency is in MHz
			m_fid.SetTransmitterFrequency(1e6 * ft);
		}

		if( it->first == "TE" ) 
		{
			it++;
			ins.str(it->first);
			treal tau;
			ins >> tau;

			// it appears echo time is in ms
			m_fid.SetEchoTime(tau/1000.0);
		}

		if( it->first == "NumberOfAverages" ) 
		{
			it++;
			ins.str(it->first);
			integer nAverages;
			ins >> nAverages;

			m_fid.SetAverages(nAverages);
		}

		if( it->first == "CSIMatrixSize[0]" ) 
        {
            it++;
			ins.str(it->first);
			integer rows;
			ins >> rows;

            m_fid.SetRows(rows);
        }

		if( it->first == "CSIMatrixSize[1]" ) 
        {
            it++;
			ins.str(it->first);
			integer cols;
			ins >> cols;
			
            m_fid.SetCols(cols);
        }
        
        if( it->first == "CSIMatrixSize[2]" ) 
        {
            it++;
			ins.str(it->first);
			integer slices;
			ins >> slices;
			
            m_fid.SetSlices(slices);
        }
        
        if( it->first == "PositionVector[0]" ) 
        {
            it++;
			ins.str(it->first);
            std::vector<double> vec;
			double tmp;
			ins >> tmp;
            vec.push_back(tmp);
            it++;
            it++;
			ins.str(it->first);
			ins >> tmp;
            vec.push_back(tmp);
            it++;
            it++;
			ins.str(it->first);
			ins >> tmp;
            vec.push_back(tmp);
            m_fid.SetPos(vec);
        }
    
        if( it->first == "RowVector[0]" ) 
        {
            it++;
			ins.str(it->first);
            std::vector<double> vec;
			double tmp;
			ins >> tmp;
            vec.push_back(tmp);
            it++;
            it++;
			ins.str(it->first);
			ins >> tmp;
            vec.push_back(tmp);
            it++;
            it++;
			ins.str(it->first);
			ins >> tmp;
            vec.push_back(tmp);
            m_fid.SetRowDirn(vec);
        }

        if( it->first == "ColumnVector[0]" ) 
        {
            it++;
			ins.str(it->first);
            std::vector<double> vec;
			double tmp;
			ins >> tmp;
            vec.push_back(tmp);
            it++;
            it++;
			ins.str(it->first);
			ins >> tmp;
            vec.push_back(tmp);
            it++;
            it++;
			ins.str(it->first);
			ins >> tmp;
            vec.push_back(tmp);
            m_fid.SetColDirn(vec);
        }

        if( it->first == "VOIThickness" ) 
        {
            it++;
			ins.str(it->first);
            std::vector<double> vec;
			double tmp;
			ins >> tmp;
            vec.push_back(tmp);
            it++;
            it++;
			ins.str(it->first);
			ins >> tmp;
            std::vector<double>::iterator insert_it;
            insert_it = vec.begin();
            vec.insert(insert_it, tmp);
            it++;
            it++;
			ins.str(it->first);
			ins >> tmp;
            insert_it = vec.begin()+1;
            vec.insert(insert_it, tmp);
            m_fid.SetVoiDim(vec);
        }
        
        if( it->first == "PixelSpacingRow" ) 
        {
            it++;
			ins.str(it->first);
            std::vector<double> vec;
			double tmp;
			ins >> tmp;
            vec.push_back(tmp);
            it++;
            it++;
			ins.str(it->first);
			ins >> tmp;
            vec.push_back(tmp);
            it++;
            it++;
			ins.str(it->first);
			ins >> tmp;
            vec.push_back(tmp);
            m_fid.SetVoxelDim(vec);
        }
	}

    // correct the position vector to be DICOM-like
    // by subtracting half a voxel row and col from position
    if ( m_fid.IsKnownVoxelDim() && m_fid.IsKnownPos() )
    {
        std::vector<double> row_dirn = m_fid.GetRowDirn();
        std::vector<double> col_dirn = m_fid.GetColDirn();
        std::vector<double> vox_dim = m_fid.GetVoxelDim();
        std::vector<double> pos = m_fid.GetPos();

        pos[0] = pos[0] + vox_dim[0] / 2.0 * row_dirn[0];
        pos[1] = pos[1] + vox_dim[0] / 2.0 * row_dirn[1];
        pos[2] = pos[2] + vox_dim[0] / 2.0 * row_dirn[2];

        pos[0] = pos[0] + vox_dim[1] / 2.0 * col_dirn[0];
        pos[1] = pos[1] + vox_dim[1] / 2.0 * col_dirn[1];
        pos[2] = pos[2] + vox_dim[1] / 2.0 * col_dirn[2];
        
        m_fid.SetPos(pos);
    }

	if( -1 == nSamples ) 
	{
		throw Exception("did not find the key field 'VectorSize'");
	}

	// now reopen file in binary mode and advance to the right point
	std::ifstream file(strFilename.c_str(), std::ios::binary);

	if( true == file.fail() )
		throw Exception("failed to open '%s' in binary mode", strFilename.c_str());
	
	file.seekg(binpos, std::ios_base::beg);

	ReadFIDData(file, nSamples);

    m_fid.swap_row_col();
}

void tarquin::CFIDReaderRDA::ReadFIDData(std::ifstream& file, std::size_t nSamples)
{
	double real_part;
	double imag_part;

	int n_fids = m_fid.GetRows() * m_fid.GetCols() * m_fid.GetSlices();

    for ( int fid = 0; fid < n_fids; fid++)
    {
        cvm::cvector FID(nSamples);
        int n = 0;
        for( std::size_t nSamplesRead = 0; nSamplesRead < nSamples; nSamplesRead++  ) 
        {
            file.read((char*)&real_part, 8);
            file.read((char*)&imag_part, 8);

            if( true == file.eof() ) 
                throw Exception("ran out file before reading in target number of samples");

            // store this sample
            FID[n+1] =  tcomplex(real_part, imag_part);
            n++;
        }
        m_fid.AppendFromVector(FID);
    }
}

std::streampos tarquin::CFIDReaderRDA::LexRDA(std::string strFilename)
{
	// attempt to open file
	std::ifstream fin(strFilename.c_str(), std::ios_base::binary);

	if( false == fin.is_open() )
		throw Exception("failed to open '%s' for lexing", strFilename.c_str());

	size_t nLine = 0;

	bool bASCIIRead = false;

	// read line by line 
	for( std::string strLine; getline(fin, strLine); ) 
	{
		nLine++;

		// skip over empty lines
		if( 0 == strLine.length() )
			continue;

		// skip over header
		if( strLine.find(">>> Begin of header <<<",0) == 0 )
			continue;

		// stop at end of ASCII bit
		if( strLine.find(">>> End of header <<<",0) == 0 ) 
		{
			bASCIIRead = true;
			break;
		}

		// use the boost tokenizer that keeps things like strings in quotes intact
		// middle argument is token separator 
		escaped_list_separator<char> sep("\\", ": ", "\"");

		// tokenize this line
		tokenizer< escaped_list_separator<char> > tokens_line(strLine, sep);

		// variables used when loading 
		tokenizer< escaped_list_separator<char> >::iterator it_token = tokens_line.begin();

		// check for bad file format, i.e. failure to tokenize
		if( it_token == tokens_line.end() ) 
			throw Exception("failed to tokenize line: %d of '%s'", nLine, strFilename.c_str());

		// add line tokens to top level collection
		for( ; it_token != tokens_line.end(); it_token++ ) 
		{
			if( it_token->length() > 0 ) 
				m_tokens.push_back(make_pair(*it_token, nLine));
		}
	}

	// if we have read the ASCII bit OK, then return the point in the file where
	// the binary bit begins
	if( true == bASCIIRead ) 
		return fin.tellg();
	else
		return -1;
}


