#ifndef __CFIDREADERSIEMENS__
#define __CFIDREADERSIEMENS__

#include "CFIDReader.hpp"
#include <fstream>

namespace tarquin {

    class CFIDReaderSiemens : public CFIDReader {

	public:

	    CFIDReaderSiemens(CFID& fid, CBoswell& log);

	    void Load(std::string strFilename, const Options& opts, CBoswell& log);

	private:

	    void ReadFIDData(std::ifstream& file, std::size_t nLength);
        
        template <class T>
            bool from_string(T& t, const std::string& s, std::ios_base& (*f)(std::ios_base&));
        
        void rotate_vec(const cvm::rvector &vec_in, const cvm::rvector &ax, double theta, cvm::rvector& vec_out);

    };

}

#endif
