#ifndef __CFIDReaderPhilipsDCM__
#define __CFIDReaderPhilipsDCM__

#include "CFIDReader.hpp"

namespace tarquin 
{

class CFIDReaderPhilips : public CFIDReader 
{

	public:

		CFIDReaderPhilips(CFID& fid, CBoswell& log);

		void Load(std::string strFilename, const Options& opts, CBoswell& log);

	private:

		void ReadParamFile(std::string strFileParam);

		void EatTokens();

		void ReadFIDFile(std::string strFilename, const Options& opts);

		void EatTokensFID(TokenList::iterator it_token, std::size_t nSamples);
    
        void rotate_vec(const cvm::rvector &vec_in, const cvm::rvector &ax, double theta, cvm::rvector& vec_out);

        // do we need to fft the data?
        bool m_FFT;
};

}

#endif
