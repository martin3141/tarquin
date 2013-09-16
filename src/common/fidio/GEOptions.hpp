#ifndef __GEOPTIONS__
#define __GEOPTIONS__

namespace tarquin {

    //! Options specific to GE files; these come from the SHF files or from guestimates.
    struct GEOptions {

	GEOptions()
	{
	    nOffset      = 162292;
	    nWaterFrames = 2;
	    nWSFrames    = 32;
	    nFieldSize   = 2048;
	    nCoils       = 8;
	    nFormat      = "Raw 15.0 PROBE Multicoil";
        nHeaderRev   = 12;
        nEchoes      = 0;
	}

	//! Offset to data in GE file.
	integer nOffset;

	//! Number of GE averages (WS).
	integer nWSFrames;

	//! Number of GE averages (W);
	integer nWaterFrames;

	//! Dimension of vector.
	integer nFieldSize;

	//! Number of coils.
	integer nCoils;
    
    //! Data format.
	std::string nFormat;

    //! Header version number
    float nHeaderRev;

    //! Number of echoes
	integer nEchoes;

	};
}

#endif
