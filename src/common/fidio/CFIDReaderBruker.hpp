#ifndef __CFIDReaderBruker__
#define __CFIDReaderBruker__

#include "CFIDReader.hpp"

namespace tarquin {

    class CFIDReaderBruker : public CFIDReader {

	public:

	    CFIDReaderBruker(CFID& fid, CBoswell& log);

	    void Load(std::string strFilename, const Options& opts, CBoswell& log);

	private:
        
        float lookup_skip_points( int decim, int dspfvs );

        float FloatSwap( float f );
        
        int LongSwap (int i);

        short ShortSwap( short s );

        template<typename T>
            std::string xxx_to_bin(const T& value);

        template <class T>
            bool from_string(T& t, const std::string& s, std::ios_base& (*f)(std::ios_base&));

    };
}

#endif
