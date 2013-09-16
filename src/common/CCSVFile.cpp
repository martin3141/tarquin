#include "CCSVFile.hpp"

using namespace std;

namespace tarquin {

	inline bool StringToDouble(std::string& str, double& d) 
	{
		istringstream myStream(str);
		if( myStream >> d )
			return true;
		else
			return false;
	}

	bool CCSVFile::load(std::string strFilename) {

		// attempt to open the file
		ifstream infile(strFilename.c_str());

		if( !infile ) {
			cerr << "failed to open " << strFilename;
			return false;
		}

        m_file_name = strFilename;

		// read file into vector of lines
		string strLine;
		vector<string> vecLines;

		while (! infile.eof() ) {
			getline(infile,strLine);
			vecLines.push_back(strLine);
		}
		infile.close(); 


		int rows = vecLines.size()-1;

		// old way ... count number of columns by looking at second line
		// count number of columns by looking at first line
		int cols = 1;	
		for(size_t n = 0; n < vecLines[0].size(); n++)
			if (vecLines[0][n] == ',')
				cols++;

		size_t p;
		size_t pter;
		string tempstr;
		vector<string> strvec;

		for(int r = 0; (r < rows); r++) {
			pter = 0;
			for(int c = 0; (c < (cols+1)); c++) {
				p = 0;

				// Greg replaced the two lines below (the order of evaluation is different on M$ compiler)
				// and the character is examined before checking that it is legal to do so
				//while ((vecLines[r][p+pter] != ',') && ((p+pter) < vecLines[r].size()))
				//	p++;
				for( p = 0; p+pter < vecLines[r].size(); p++ )
					if( vecLines[r][p+pter] == ',' )
						break;

				tempstr = "";
				for (size_t q = 0; (q < p); q++)
					tempstr += vecLines[r][pter+q];

				strvec.push_back(tempstr);
				pter = pter + p + 1;
			}
			m_strMat.push_back(strvec);
			strvec.clear();
		}	
		vector<double> doubvec;

		// convert string matrix to double matrix
		double result;
		for(int r = 0; (r < rows); r++) {
			for(int c = 0; (c < cols); c++) {
				if (StringToDouble(m_strMat[r][c],result)) {
					doubvec.push_back(result);
					m_strMat[r][c] = "";
				}
				else
					doubvec.push_back(std::numeric_limits<double>::quiet_NaN());
			}
			m_doubMat.push_back(doubvec);
			doubvec.clear();
		}

		return true;
	}

	// output to stream; conventionally this is a friend not a member
	ostream& operator<<(ostream& os, const CCSVFile& rhs)
	{
		// set width and precision for output stream
		os.setf(ios::fixed | ios::right);
		os.precision(2);
		os.width(4);
		os.fill(' ');

		os << "Double matrix:";
		for( vector<vector<double> >::const_iterator itRow = rhs.m_doubMat.begin(); itRow != rhs.m_doubMat.end(); itRow++ ) {
			os << endl;

			for( vector<double>::const_iterator itCol = itRow->begin(); itCol != itRow->end(); itCol++ ) 
				os << *itCol << " ";
		}

		os << endl << "String matrix:";
		for( vector<vector<string> >::const_iterator itRow = rhs.m_strMat.begin(); itRow != rhs.m_strMat.end(); itRow++ ) {
			os << endl;

			for( vector<string>::const_iterator itCol = itRow->begin(); itCol != itRow->end(); itCol++ ) 
				os << *itCol << " ";
		}

		return os;
	}
}
