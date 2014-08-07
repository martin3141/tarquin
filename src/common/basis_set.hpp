#ifndef BASIS_PARAMETER_HARDCODED_INCLUDED
#define BASIS_PARAMETER_HARDCODED_INCLUDED

#include <vector>
#include <string>
#include <limits>
#include <math.h>
#include <iostream>

namespace tarquin 
{
	enum basis_vector_e
	{
		BV_CRCH2,
		BV_CRCH2_RT,
		BV_ALA,
		BV_ASP,
		BV_CR,
		BV_CR_RT,
		BV_CHO_RT,
		BV_PCR,
		BV_GABA,
		BV_MEGAGABA,
		BV_GLC,
		BV_GLN,
		BV_GLU,
		BV_GLU_RT,
		BV_GPC,
		BV_GUA,
		BV_INS,
		BV_INS_RT,
		BV_LAC,
		BV_LAC_RT,
		BV_LIP09,
		BV_LIP13A,
		BV_LIP13B,
		BV_LIP20,
		BV_MM09,
		BV_MM12,
		BV_MM14,
		BV_MM17,
		BV_MM20,
		BV_MM30,
		BV_MM38,
		BV_NAA,
		BV_NAA_RT,
		BV_NAAG,
		BV_PCH,
		BV_SCYLLO,
		BV_TAU,
		BV_GLY,
		BV_CIT,
		BV_PETH,

        BV_GLTH,

		BV_TCHO,
		BV_TNAA,
		
        BV_GABA_A,
        BV_GABA_B,
        BV_NAA_MEGA,
        BV_GLX_A,
        BV_GLX_B,
        BV_GLX_C,
        BV_GLX_D,
        BV_MM09_MEGA,

		BV_31P_ATP,
		BV_31P_GPC,
		BV_31P_GPE,
		BV_31P_NADH,
		BV_31P_PCH,
		BV_31P_PCR,
		BV_31P_PE,
		BV_31P_PI
	};

	
	void getMetaboliteMatrix(
			basis_vector_e metab,
			std::vector<std::vector<double> >& doubmat,
			double tf
			);

	void getMetaboliteDescription(basis_vector_e metab, std::string& desc);

} // namespace tarquin

#endif // BASIS_PARAMETER_HARDCODED_INCLUDED
