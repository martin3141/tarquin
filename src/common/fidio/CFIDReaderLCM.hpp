#ifndef __CFIDReaderLCM__
#define __CFIDReaderLCM__

#include "CFIDReader.hpp"

namespace tarquin 
{

class CFIDReaderLCM : public CFIDReader 
{

	public:

		CFIDReaderLCM(CFID& fid, CBoswell& log);

		void Load(std::string strFilename, const Options& opts, CBoswell& log);

	private:

		void EatTokens();
};

}

#endif
