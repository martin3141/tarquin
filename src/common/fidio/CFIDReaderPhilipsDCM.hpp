#ifndef __CFIDREADERPHILPSDCM__
#define __CFIDREADERPHILPSDCM__

#include "CFIDReader.hpp"
#include <fstream>

namespace tarquin {

    class CFIDReaderPhilipsDCM : public CFIDReader {

	public:

	    CFIDReaderPhilipsDCM(CFID& fid, CBoswell& log);

	    void Load(std::string strFilename, const Options& opts, CBoswell& log);

	private:

	    void ReadFIDData(std::ifstream& file, std::size_t nLength);
        
        template <class T>
            bool from_string(T& t, const std::string& s, std::ios_base& (*f)(std::ios_base&));



    };

}

#endif
