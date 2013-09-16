#ifndef __CCSVFILE__
#define __CCSVFILE__

#include "../common/common.hpp"
#include <iostream>

namespace tarquin {

	/*!
	 * A class for representing basis-vector descriptions, as contained in the comma-separated-value files.
	 *
	 * In the original code, this was called "Signal".
	 */

	class CCSVFile {

		public:
			CCSVFile()
			{

			}

			~CCSVFile()
			{

			}


			//! Load a CSV file into the class.
			bool load(std::string strFilename);

			//! Dump to a stream.
      friend std::ostream& operator<<(std::ostream& os, const CCSVFile& rhs);

		private:

			//! Matrix of doubles, first dimension is line in file?
			std::vector< std::vector<double> > m_doubMat;

			//! Matrix of strings, first dimension is line?
			std::vector< std::vector<std::string> > m_strMat;
			
            std::string m_file_name;

			//
			// Accessors below here.
			//
		public:

			//! Get a reference to the "matrix" of doubles.
			inline std::vector< std::vector<double> >& getDoubleMatrix()
			{
				return m_doubMat;
			}

			//! Get a reference to the "matrix" of strings.
			inline std::vector< std::vector<std::string> >& getStringMatrix()
			{
				return m_strMat;
			}
    

            inline std::string getFileName()
            {
                return m_file_name;
            }

	};

}

#endif
