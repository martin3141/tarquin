#ifndef __CFIDREADERDCM__
#define __CFIDREADERDCM__

#include "CFIDReader.hpp"
#include <fstream>

namespace tarquin {

    class CFIDReaderDCM : public CFIDReader {

	public:

	    CFIDReaderDCM(CFID& fid, CBoswell& log);

	    void Load(std::string strFilename, const Options& opts, CBoswell& log);
	    
        void LoadW(std::string strFilename, const Options& opts);

	private:

	    void ReadFIDData(std::ifstream& file, std::size_t nLength, std::size_t byte_offset, std::size_t averages);
        
        template <class T>
            bool from_string(T& t, const std::string& s, std::ios_base& (*f)(std::ios_base&));

        void rotate_vec(const cvm::rvector &vec_in, const cvm::rvector &ax, double theta, cvm::rvector& vec_out);



    };

}

#endif
