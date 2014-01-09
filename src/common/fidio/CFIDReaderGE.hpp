#ifndef __CFIDReaderGE__
#define __CFIDReaderGE__

#include "CFIDReader.hpp"
#include "GEOptions.hpp"
#include "Options.hpp"

namespace tarquin 
{

class CFIDReaderGE : public CFIDReader 
{

	private:

		//! Offsets into the GE file, etc.
		GEOptions m_options;

	public:

		CFIDReaderGE(CFID& fid, CBoswell& log);

		/*!
		 * Use this function if the SHF files are not available, it may fail since
		 * the parameters are only guestimates.
		 */
		void Load(std::string strFilename, const Options& opts, CBoswell& log);

		void LoadW(std::string strFilename, const Options& opts, CBoswell& log);

		/*!
		 * Use this function if the SHF files are provided, give it the name of the
		 * first SHF file.
		 */
		void LoadSHF(std::string strFilename, const Options& opts, bool WS, CBoswell& log);
		
		/*!
		 * If these options are set from the command line, initialise the structure
		 */
		void SetOptions(const GEOptions& opts)
		{
			m_options = opts;
		}

	private:

		void DiscoverOptions(std::string strFilename, CBoswell& log); 

		//! Turns the token stream into the options structure (if loading SHF files).
		void EatTokens(CBoswell& log);

		//! This is called once the options structure has been initialised.
		void LoadFromOptions(std::string strFilename, const Options& opts, bool WS);

		void LoadFromOptionsCSI(std::string strFilename);

		void LoadFromOptionsSVS(std::string strFilename, const Options& opts, bool WS);
		
        int LongSwap (int i);

        short ShortSwap( short s );

        float FloatSwap( float f );
};

}

#endif
