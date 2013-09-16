#ifndef __CFIDReaderRDA__
#define __CFIDReaderRDA__

#include "CFIDReader.hpp"

namespace tarquin 
{

class CFIDReaderRDA : public CFIDReader 
{

	public:

		CFIDReaderRDA(CFID& fid, CBoswell& log);

		void Load(std::string strFilename, const Options& opts, CBoswell& log);

	private:

		std::streampos LexRDA(std::string strFilename);

		void ReadFIDData(std::ifstream& file, std::size_t nSamples);
};

}

#endif
