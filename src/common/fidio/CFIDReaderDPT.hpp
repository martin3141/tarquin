#ifndef __CFIDReaderDPT__
#define __CFIDReaderDPT__

#include "CFIDReader.hpp"

namespace tarquin {

    class CFIDReaderDPT : public CFIDReader {

	public:

	    CFIDReaderDPT(CFID& fid, CBoswell& log);

	    void Load(std::string strFilename, const Options& opts, CBoswell& log);

	private:

	    void EatTokens();

	    void EatTokensFID(TokenList::iterator it_token, std::size_t nSamples);

    };
}

#endif
