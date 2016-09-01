#ifndef __CFIDReaderJMRUITXT__
#define __CFIDReaderJMRUITXT__

#include "CFIDReader.hpp"
#include "CFID.hpp"

namespace tarquin {

    class CFIDReaderJMRUITXT : public CFIDReader {

	public:

	    CFIDReaderJMRUITXT(CFID& fid, CBoswell& log);

	    void Load(std::string strFilename, const Options& opts, CBoswell& log);

	private:

	    void EatTokens(const Options& opts);

	    void EatTokensFID(TokenList::iterator it_token, std::size_t nSamples, const Options& opts, TokenList::iterator end_token, int skip);

    };
}

#endif
