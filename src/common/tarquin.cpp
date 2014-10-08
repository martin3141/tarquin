/*!
 * All the code uses the "levmar" library to do the optimisation, but we provide the objective function, Jacobian
 * and constraints. 
 *
 * The objective function makes use of (Greg's) implementation of the Non-Negative Least Squares algorithm. It
 * should probably be replaced by TNNLS.
 * 
 * The implementation here is by no means efficient as it could be, a better implementation must derive
 * from CNNLSSolver in order to be compatible with this code. TODO: Remove this dependence, it is silly.
 *
 * Notation: suffixed by 'p'  = prime (e.g. Yp = Y')
 *           suffixed by 'pp' = double prime (e.g. Ypp = Y'')
 */

#include <iostream>
#include <set>
#include <cmath>
#include <string.h>
#include "../levmar/levmar-2.5/levmar.h" // this conflicts with some other header on win32, so its path must be explicit 
#include "tarquin.hpp"
#include "CNNLSSolverLH.hpp"
#include "NNLSSolverBJ.hpp"
#include "FNNLSSolverBJ.hpp"
#include "td_conv_ws.hpp"


namespace tarquin 
{

void residual_objective_all(
		treal* pp, 
		treal* px, 
		integer nParams, 
		integer activeNdbl, 
		void* pParams
		);

void jacobian_func(
		treal* pp, 
		treal* pjac, 
		integer nParams, 
		integer activeNdbl, 
		void* pParams
		);

/*
 * This is for holding the parameters to the object function - i.e. all the variables that
 * must persist between calls to the objective function.
 */	
struct SResParams 
{

	SResParams(
			const cvm::cvector& y,
			const CBasis& basis,
			cvm::rvector& yActive, 
			const cvm::cmatrix& G, 
			const cvm::cmatrix& GDFT, 
			const cvm::rmatrix& A,
			const cvm::rvector& t, 
			cvm::rvector& ahat, 
			integer nStart,
			cvm::rvector& shift_freq_range,
			integer activeN,
			integer M,
			integer Q,
			std::vector<std::string>& metab_names,
			treal lambda,
			CBoswell& log,
            Options& opts ) : 
		m_y(y),
		m_basis(basis),
		m_yActive(yActive),
		m_G(G),
		m_GDFT(GDFT),
		m_A(A),
		m_t(t),
		m_ahat(ahat),
		m_nStart(nStart),
		m_shift_freq_range(shift_freq_range),
		m_metab_names(metab_names),
		m_lambda(lambda),
		m_log(log),
		m_opts(opts)
	{ 
		// modified group matrix [real; imag]
		m_Gp.resize(2*activeN, M);
		m_Gpp.resize(2*activeN, M);

		m_gc.resize(m_G.msize());
		m_gcp.resize(m_G.msize());
		m_GC.resize(m_G.msize());
		m_GCp.resize(m_G.msize());

		// modified basis matrix [real; imag]
		m_Sp.resize(2*activeN, Q);
		m_Spp.resize(2*activeN, Q);

		// intermediate matrix for Jacobian
		m_dX.resize(Q, 2*activeN);

		// intermediate matrix for Jacobian
		m_dSaug.resize(2*activeN, Q);

		// intermediate matrix for Jacobian
		m_dSaugF.resize(2*activeN, Q);

		// final Jacobian (transposed)
		integer P = M*2 + 1 + 1 + 1;
		m_final_jac.resize(2*activeN, P);
	}

	// the unmodified original FID
	const cvm::cvector& m_y;

	//! The basis used for anaylsis.
	const CBasis& m_basis;

	// chunk of signal we are fitting to [real; imag]
	cvm::rvector& m_yActive;

	// orginal group matrix
	const cvm::cmatrix& m_G;

	// Fourier transform of each column
	const cvm::cmatrix& m_GDFT;

	// phased column of this matrix
	cvm::cvector m_gc;
	cvm::cvector m_gcp;

	// DFT of phased column of group matrix
	cvm::cvector m_GC;
	cvm::cvector m_GCp;

	// summation matrix
	const cvm::rmatrix& m_A;

	// time
	const cvm::rvector& m_t;

	// amplitude estimates
	cvm::rvector& m_ahat;

	// starting sample
	integer m_nStart;

	// frequency range
	const cvm::rvector& m_shift_freq_range;

	// modified group matrix [real; imag]	
	cvm::rmatrix m_Gp;	

	// modified group matrix [real; imag] multiplied by frequency term	
	cvm::rmatrix m_Gpp;	

	// modified basis matrix [real; imag]
	cvm::rmatrix m_Sp;	

	// modified basis matrix [real; imag] multiplied by frequency term
	cvm::rmatrix m_Spp;

	// metabolite names
	std::vector<std::string> m_metab_names;	

    // soft cons value
    treal m_lambda;

	// intermediate matrix used for Jacobian
	cvm::rmatrix m_dX;

	// intermediate matrix used for Jacobian
	cvm::rmatrix m_dSaug;

	// intermediate matrix used for Jacobian
	cvm::rmatrix m_dSaugF;

	// final Jacobian (we keep for CRLB calculation)
	cvm::rmatrix m_final_jac;

	// logging object
	CBoswell& m_log;
    
    // options
    const Options& m_opts;
};

} // namespace tarquin

/*
 * This is the objective function. This is called by the optimiser repeatedly.
 * generate f(p), the levmar code will minimise ||x-f(p)||, where x is the observed signal, 
 * in this case yActive.
 */
void tarquin::residual_objective_all(
		treal* pp, 
		treal* px, 
		integer nParams, 
		integer activeNdbl, 
		void* pParams
		)
{
    // following added to stop compiler warnings
    activeNdbl = activeNdbl;

	SResParams& params = *reinterpret_cast<SResParams*>(pParams);

	// vector of parameters (shares memory)
	cvm::rvector vParams(pp, nParams);

	// some variables we use locally
	integer M         = params.m_G.nsize();
	integer nStart    = params.m_nStart;
	integer nActive   = params.m_yActive.size() / 2;
	cvm::cvector& gc  = params.m_gc;
	cvm::cvector& gcp = params.m_gcp;
	cvm::cvector& GC  = params.m_GC;
	cvm::cvector& GCp = params.m_GCp;
	const cvm::rvector& shift_freq_range = params.m_shift_freq_range;
	CBoswell& log     = params.m_log;
	const Options& opts     = params.m_opts;

	//
	// apply parameters to group matrix (all parameters except phi0, phi1)
	//
	
	// indices into the parameter array of these parameters
	std::size_t nIdxBeta  = 2*M+1;
	std::size_t nIdxPhi0  = 2*M+2; 
	std::size_t nIdxPhi1  = 2*M+3;

	// parameters common to all columns
	treal phi1 = vParams(nIdxPhi1);
	treal phi0 = vParams(nIdxPhi0); 

	// compute the exp(j\phi_1 \omega) term to be applied to the frequency domain of each column
	cvm::cvector linear_phase(GC.size());

	for( integer n = 0; n < GC.size(); ++n )
	{
		treal freq = shift_freq_range[n+1];
		//treal dw   = 2.0*M_PI*freq;
		linear_phase[n+1] = std::exp( tcomplex(0, phi1*2.0*M_PI*freq) );
	}

	// we only want to compute a plan once, so set these up here
	ifft_wrap ifft_GC_gc(GC, gc);
	ifft_wrap ifft_GCp_gcp(GCp, gcp);
	
	// for each column
	for( integer c = 1; c <= M; c++ ) 
	{
		// apply phi1 term to frequency domain
		for( integer n = 0; n < GC.size(); n++ ) 
		{
			treal dw = 2.0*M_PI*shift_freq_range[n+1];

			// the phi1 term is the same for all columns of the group matrix
			GC[n+1]  = params.m_GDFT(n+1, c) * linear_phase[n+1];
			GCp[n+1] = tcomplex(-GC[n+1].imag()*dw, GC[n+1].real()*dw);
		}

		//ifft(GC, gc);
		//ifft(GCp, gcp);
		ifft_GC_gc();
		ifft_GCp_gcp();

		// indices to parameter vector	
		std::size_t nIdxShift = c;
		std::size_t nIdxAlpha = c+M;
	
		// the parameters to apply
		treal omega = vParams(nIdxShift)*2.0*M_PI;
		treal alpha = vParams(nIdxAlpha);
		treal beta  = pow(vParams(nIdxBeta),1)*opts.GetBetaScale(); // MAGIC BETA SCALE

		// for each (sample of the active part of the signal == row of G)
		for( integer n = 0; n < nActive; n++ ) 
		{
			// time starts at the offset because of nStart
			treal time = params.m_t(n+nStart);

            treal xx    = time*omega + phi0;
            treal ct    = std::cos(xx);
            treal st    = std::sin(xx);
            treal decay = std::exp(-time*(alpha + time*beta));

			// common sub expression
			//tcomplex zstar = std::exp(-time*(alpha + time*beta) + tcomplex(0, time*omega + phi0));
            tcomplex zstar(decay*ct, decay*st);

			// the modified sample of the group matrix 
			tcomplex z  = gc(n+nStart)  * zstar;
			tcomplex zp = gcp(n+nStart) * zstar;

			// modify sample and put in right place in active matrix
			params.m_Gp(n+1, c)         = z.real();
			params.m_Gp(n+1+nActive, c) = z.imag();	

			params.m_Gpp(n+1, c)         = zp.real();
			params.m_Gpp(n+1+nActive, c) = zp.imag();	
		}
	}

	// construct basis matrix from group matrix
	params.m_Sp  = params.m_Gp  * params.m_A;
	params.m_Spp = params.m_Gpp * params.m_A;

	//
	// generate candidate fit 
	//

	//CNNLSSolverLH nnls;
    //NNLSSolverBJ nnls;
    FNNLSSolverBJ nnls;

	// estimate amplitudes


	cvm::rvector yActive_aug = params.m_yActive;
	cvm::rmatrix Sp_aug = params.m_Sp;
    
    bool soft_cons = true;

	if ( soft_cons )
	{
		// use inital estimate to produce an estimate with soft constraints
		// NAA and NAAG
		int naa_indx = -1;
		int naag_indx = -1;
		int found = 0;

		for ( size_t n = 0; n < params.m_metab_names.size(); n++ )
		{
			if ( params.m_metab_names[n] == "NAA" )
			{
				found++;
				naa_indx = n+1;
			}
			if ( params.m_metab_names[n] == "NAAG" )
			{
				found++;
				naag_indx = n+1;
			}
		}
		// if both signals were found, proceed with nnls augmentation
		if ( found == 2 )
		{
			Sp_aug.resize(Sp_aug.msize()+2, Sp_aug.nsize());
			yActive_aug.resize(yActive_aug.size()+2);
			treal TNAA = params.m_ahat(naa_indx) + params.m_ahat(naag_indx);
			//float NAA = params.m_ahat(naa_indx);

			//float Sp_norm = params.m_Sp.norm1();
			// NAAG
			treal lambda = 0.05;
			//lambda = 0.01;
			//lambda = 10;
			yActive_aug(yActive_aug.size()-1) = TNAA*(0.15)/(0.15+1)*lambda*params.m_Sp(naag_indx).norm2();
			Sp_aug(Sp_aug.msize()-1, naag_indx) = lambda*params.m_Sp(naag_indx).norm2();

			// NAA 
			yActive_aug(yActive_aug.size()) = TNAA/(0.15+1)*lambda*params.m_Sp(naa_indx).norm2();
			Sp_aug(Sp_aug.msize(), naa_indx) = lambda*params.m_Sp(naa_indx).norm2();
		}

		// Lipids 
		int lip09_indx  = -1;
		int lip13a_indx = -1;
		int lip13b_indx = -1;
		int lip20_indx  = -1;
		found = 0;

		for ( size_t n = 0; n < params.m_metab_names.size(); n++ )
		{
			if ( params.m_metab_names[n] == "Lip09" )
			{
				found++;
				lip09_indx = n+1;
			}

			if ( params.m_metab_names[n] == "Lip13a" )
			{
				found++;
				lip13a_indx = n+1;
			}

			if ( params.m_metab_names[n] == "Lip13b" )
			{
				found++;
				lip13b_indx = n+1;
			}

			if ( params.m_metab_names[n] == "Lip20" )
			{
				found++;
				lip20_indx = n+1;
			}

		}
		// if all signals were found, proceed with nnls augmentation
		if ( found == 4 )
		{
			Sp_aug.resize(Sp_aug.msize()+2, Sp_aug.nsize());
			yActive_aug.resize(yActive_aug.size()+2);
			treal Lip13 = params.m_ahat(lip13a_indx) + params.m_ahat(lip13b_indx);
			//float Sp_norm = params.m_Sp.norm1();
			// Lip09 
			treal lambda = 0.5;
			//lambda = 5;
			//lambda = 0.05;			
			//yActive_aug(yActive_aug.size()-1) = Lip13*0.12*lambda*params.m_Sp(lip09_indx).norm2();
			yActive_aug(yActive_aug.size()-1) = Lip13*0.267*lambda*params.m_Sp(lip09_indx).norm2();
			Sp_aug(Sp_aug.msize()-1, lip09_indx) = lambda*params.m_Sp(lip09_indx).norm2();

			// Lip20
			lambda = 0.3;
			yActive_aug(yActive_aug.size()) = Lip13*0.15*lambda*params.m_Sp(lip20_indx).norm2();
			Sp_aug(Sp_aug.msize(), lip20_indx) = lambda*params.m_Sp(lip20_indx).norm2();
		}

		// MMs
		int MM09_indx = -1;
		int MM12_indx = -1;
		int MM14_indx = -1;
		int MM17_indx = -1;
		int MM20_indx = -1;
		found = 0;

		for ( size_t n = 0; n < params.m_metab_names.size(); n++ )
		{
			if ( params.m_metab_names[n] == "MM09" )
			{
				found++;
				MM09_indx = n+1;
			}

			if ( params.m_metab_names[n] == "MM12" )
			{
				found++;
				MM12_indx = n+1;
			}

			if ( params.m_metab_names[n] == "MM14" )
			{
				found++;
				MM14_indx = n+1;
			}

			if ( params.m_metab_names[n] == "MM17" )
			{
				found++;
				MM17_indx = n+1;
			}

			if ( params.m_metab_names[n] == "MM20" )
			{
				found++;
				MM20_indx = n+1;
			}

		}
		// if all signals were found, proceed with nnls augmentation
		if ( found == 5 )
		{
			//Sp_aug.resize(Sp_aug.msize()+5, Sp_aug.nsize());
			//yActive_aug.resize(yActive_aug.size()+5);
			Sp_aug.resize(Sp_aug.msize()+4, Sp_aug.nsize());
			yActive_aug.resize(yActive_aug.size()+4);
			//float TMM = params.m_ahat(MM09_indx) + params.m_ahat(MM12_indx)
		    //		+ params.m_ahat(MM14_indx) + params.m_ahat(MM17_indx) +
		    //		params.m_ahat(MM20_indx);
			float MM09 = params.m_ahat(MM09_indx);
			//float Sp_norm = params.m_Sp.norm1();
			float lambda;

			// MM09
			//lambda = 0.05;
			//lambda = 0.02;
			//yActive_aug(yActive_aug.size()-4) = TMM*(1/3.925)*lambda*Sp_norm/params.m_yActive.size();
			//Sp_aug(Sp_aug.msize()-4, MM09_indx) = lambda*Sp_norm/params.m_yActive.size();

			// MM12
			//lambda = 0.05;
			lambda = 0.9f;
			//yActive_aug(yActive_aug.size()-3) = TMM*(0.3/3.925)*lambda*Sp_norm/params.m_yActive.size();
			yActive_aug(yActive_aug.size()-3) = MM09*0.3*lambda*params.m_Sp(MM12_indx).norm2();
;
			Sp_aug(Sp_aug.msize()-3, MM12_indx) = lambda*params.m_Sp(MM12_indx).norm2();

			// MM14
			//lambda = 0.01;
			lambda = 0.3f;
			//yActive_aug(yActive_aug.size()-2) = TMM*(0.75/3.925)*lambda*Sp_norm/params.m_yActive.size();
			yActive_aug(yActive_aug.size()-2) = MM09*0.75*lambda*params.m_Sp(MM14_indx).norm2();

			Sp_aug(Sp_aug.msize()-2, MM14_indx) = lambda*params.m_Sp(MM14_indx).norm2();


			// MM17
			//lambda = 0.02;
			lambda = 0.8f;
			//yActive_aug(yActive_aug.size()-1) = TMM*(0.375/3.925)*lambda*Sp_norm/params.m_yActive.size();
			yActive_aug(yActive_aug.size()-1) = MM09*0.375*lambda*params.m_Sp(MM17_indx).norm2();
			Sp_aug(Sp_aug.msize()-1, MM17_indx) = lambda*params.m_Sp(MM17_indx).norm2();

			// MM20
			lambda = 0.2f;
			//yActive_aug(yActive_aug.size()) = TMM*(1.5/3.925)*lambda*Sp_norm/params.m_yActive.size();
			yActive_aug(yActive_aug.size()) = MM09*1.5*lambda*params.m_Sp(MM20_indx).norm2();
			Sp_aug(Sp_aug.msize(), MM20_indx) = lambda*params.m_Sp(MM20_indx).norm2();

		}
        
        //treal lambda = 0.0006;
        //treal lambda = 0.00025;
       

		// if the matrix has been adjusted do NNLS with
		// soft contraints added
		//if ( yActive_aug.size() > params.m_yActive.size() )


		//std::cout << params.m_ahat << std::endl;
		//
	}
    
    //treal lambda = 0.05;
    treal lambda = params.m_lambda;

    //float Sp_norm = params.m_Sp.norm1();
    //std::cout << params.m_yActive.size();
    //float y_norm = params.m_yActive.norm2();
    //std::cout << "Y norm :" << y_norm << std::endl;
    Sp_aug.resize(Sp_aug.msize()+1, Sp_aug.nsize());
    yActive_aug.resize(yActive_aug.size()+1);
    for ( int n = 1; n < Sp_aug.nsize() + 1; n++)
    {
        //Sp_aug(Sp_aug.msize(), n) = lambda*params.m_Sp(n).norm2(); // old method, lambda approx 0.05
        //Sp_aug(Sp_aug.msize(), n) = lambda/params.m_Sp(n).norm1(); // 
        //Sp_aug(Sp_aug.msize(), n) = lambda/params.m_Sp(n).norm1()/params.m_Sp(n).norm1(); // approx 10000
        //Sp_aug(Sp_aug.msize(), n) = params.m_Sp(n).norminf()*lambda/params.m_Sp(n).norm2()/params.m_Sp(n).norm2(); // lambda approx 10
        //Sp_aug(Sp_aug.msize(), n) = lambda/params.m_Sp(n).norm2(); // lambda approx 10
        Sp_aug(Sp_aug.msize(), n) = params.m_Sp(n).norminf()*lambda; // lambda approx 0.2
        //Sp_aug(Sp_aug.msize(), n) = 1.0/params.m_Sp(n).norminf()*lambda/params.m_Sp(n).norm2(); // lambda approx 8
        //Sp_aug(Sp_aug.msize(), n) = params.m_Sp(n).norminf()*lambda*params.m_Sp(n).norm2(); // lambda approx ?
        //Sp_aug(Sp_aug.msize(), n) = lambda/params.m_Sp(n).norminf(); // lambda approx 0.5 not so good
        //Sp_aug(Sp_aug.msize(), n) = lambda; // 0.7 is good
        //plot(temp);
    }

    if ( soft_cons)
        nnls.solve(Sp_aug, yActive_aug, params.m_ahat, opts.GetFastFit());
    else
        nnls.solve(params.m_Sp, params.m_yActive, params.m_ahat, opts.GetFastFit());
    //
	//nnls.solve(params.m_Sp, params.m_yActive, params.m_ahat, true);
	
    // estimate signal - yhat (shares memory with px provided by levmar)
	cvm::rvector yhat( px, 2*nActive );

	// compute yhat (result is in px as well)
	yhat = params.m_Sp * params.m_ahat;
    
    // test code
    /*
    cvm::cvector yhat_full(yhat.size()/2);
	cvm::cvector y_full(yhat.size()/2);
    for ( int n = 1; n < yhat.size()/2 + 1; n++ )
    {
       yhat_full(n) = tcomplex(yhat(n), yhat(n+yhat.size()/2));
       y_full(n) = tcomplex(params.m_yActive(n), params.m_yActive(n+yhat.size()/2));
    }

    //plot(yhat_full);

    cvm::cvector linesh(yhat_full.size());
    for ( int n = 1; n < yhat_full.size() + 1; n++ ) 
    {
        linesh(n) = y_full(n)/(yhat_full(n)*0.999+y_full(n)*0.001);
    }
    
    cvm::cvector linesh_smo(yhat_full.size());
    td_conv_ws( linesh, linesh_smo, 100, 10);	
    //td_conv_ws( linesh, linesh_smo, yfid.GetSamplingFrequency()/20, 10);	

    for ( int n = 1; n < yhat_full.size() + 1; n++ ) 
        yhat_full(n) = yhat_full(n) * linesh_smo(n);

    for ( int n = 1; n < yhat.size()/2 + 1; n++ )
    {
       yhat(n) = yhat_full(n).real();
       yhat(n+yhat.size()/2) = yhat_full(n).imag();
    }

    cvm::cmatrix m_Sp_full(params.m_Sp.msize()/2,params.m_Sp.nsize());
    for ( int n = 1; n < params.m_Sp.nsize() + 1; n++ )
        for ( int m = 1; m < params.m_Sp.msize()/2 + 1; m++ )
           m_Sp_full(m,n) = tcomplex(params.m_Sp(m,n), params.m_Sp(m+yhat.size()/2,n));

    for ( int n = 1; n < params.m_Sp.nsize() + 1; n++ )
        for ( int m = 1; m < m_Sp_full.msize() + 1; m++ ) 
            m_Sp_full(m,n) = m_Sp_full(m,n) * linesh_smo(m);

    for ( int n = 1; n < params.m_Sp.nsize() + 1; n++ )
    {
        for ( int m = 1; m < m_Sp_full.msize()/2 + 1; m++ )
        {
            params.m_Sp(m,n) = m_Sp_full(m,n).real();
            params.m_Sp(m+params.m_Sp.msize()/2,n) = m_Sp_full(m,n).imag();
        }
    }

    cvm::cmatrix m_Spp_full(params.m_Spp.msize()/2,params.m_Spp.nsize());
    for ( int n = 1; n < params.m_Spp.nsize() + 1; n++ )
        for ( int m = 1; m < params.m_Spp.msize()/2 + 1; m++ )
           m_Spp_full(m,n) = tcomplex(params.m_Spp(m,n), params.m_Spp(m+yhat.size()/2,n));

    for ( int n = 1; n < params.m_Spp.nsize() + 1; n++ )
        for ( int m = 1; m < m_Spp_full.msize() + 1; m++ ) 
            m_Spp_full(m,n) = m_Spp_full(m,n) * linesh_smo(m);

    for ( int n = 1; n < params.m_Spp.nsize() + 1; n++ )
    {
        for ( int m = 1; m < m_Spp_full.msize()/2 + 1; m++ )
        {
            params.m_Spp(m,n) = m_Spp_full(m,n).real();
            params.m_Spp(m+params.m_Spp.msize()/2,n) = m_Spp_full(m,n).imag();
        }
    }
    */

	log.UpdateTask(".");

	cvm::rvector r(yhat.size());
	r = yhat - params.m_yActive;
	//std::cout << "\nresidual = " << r.norm2() << std::flush;
	//std::cout << "\nphi0 = " << phi0;
	//std::cout << "\nphi1 = " << phi1;
}

/*!
 * \brief Function to compute the Jacobian.
 */
void tarquin::jacobian_func(
		treal* pp, 
		treal* pjac, 
		integer nParams, 
		integer activeNdbl, 
		void* pParams
		)
{
    // following added to stop compiler warnings
    activeNdbl = activeNdbl;

	SResParams& params = *reinterpret_cast<SResParams*>(pParams);

	// vector of parameters (shares memory)
	cvm::rvector vParams(pp, nParams);

	// some variables we use locally
	integer M       = params.m_G.nsize();
	integer P       = vParams.size();
	integer Q       = params.m_Sp.nsize();
	integer nStart  = params.m_nStart;
	integer nActive = params.m_yActive.size() / 2;

	// the complex Jacobian would be have activeN rows and P columns 
	// the real Jacboian is [real; imag] so has 2*activeN rows and P columns
	// However, because FORTRAN arrays are the opposite away round from C++,
	// this matrix is actually the tranpose of what I defined on paper!
	cvm::rmatrix J(pjac, P, 2*nActive);

	// solution of non-negative least squares problem
	cvm::rvector& annls = params.m_ahat;

	// find inactive set (i.e. non-zero amplitudes of columns of S)
	std::set<integer> F;
	for( integer n = 1; n <= Q; n++ )
		if( params.m_ahat(n) > 0.0 ) 
			F.insert(n);

	// the basis set has been adjusted so that non of the basis vectors have a > 0 amplitude
	// at this point, we assume everything is bad and quite
	// possibly, we could set the Jacobian to be all zeros here, and then the exit may be more
	// elegant
	if( F.size() == 0 ) 
	{

		params.m_log.LogMessage(LOG_ERROR, "Warning. Data quality is very poor.");
		//params.m_log.LogMessage(LOG_ERROR, "Failed, data quality is very poor, not fitting.");
		//throw std::runtime_error("data quality too low to fit");
		return;
	}

	// construct SaugF = Saug(:,F) ( can be initialised in params to be full size, take a reference here )
	cvm::rmatrix SaugF(2*nActive, F.size());
	integer c = 1;
	for( std::set<integer>::iterator itf = F.begin(); itf != F.end(); itf++, c++ ) 
		SaugF(c) = params.m_Sp(*itf);


	// construct Y = SaugF'*SaugF 
	//cvm::rmatrix Y(F.size(), F.size());
	//Y = (~SaugF)*SaugF;

	// pseudo-inverse of Y
	//cvm::rmatrix W = Y.pinv(1e-13); 

	// this is stored externally
	cvm::rmatrix& dSaug = params.m_dSaug;

	// for each parameter
	for(integer q = 1; q <= P; q++ ) 
	{
		// 0 otherwise 
		dSaug.set(0);

		// frequency adjustment parameter
		if( (1 <= q) && (q <= M) ) 
		{
			// column of S to which this adjustment applies
			integer i = params.m_basis.GetBasisFromGroup(q-1);

			// column of G which this parameter effects
			integer k = q;

			for( integer n = 0; n < nActive; n++ ) 
			{
				treal time = params.m_t(n+nStart);

				tcomplex x = tcomplex(params.m_Gp(n+1, k), params.m_Gp(n+1+nActive, k));
				tcomplex z = tcomplex(0, 2.0 * M_PI * time) * x;

				dSaug(n+1, i)         = z.real();
				dSaug(n+1+nActive, i) = z.imag();
			} 
		}

		// alpha adjustment parameter
		else if( (M+1 <= q) && (q <= 2*M) ) 
		{
			// column of S to which this adjustment applies
			integer i = params.m_basis.GetBasisFromGroup(q-1-M);

			// column of G which this parameter effects
			integer k = q-M;

			for( integer n = 0; n < nActive; n++ ) 
			{
				treal time = params.m_t(n+nStart);

				tcomplex x = tcomplex(params.m_Gp(n+1, k), params.m_Gp(n+1+nActive, k));
				tcomplex z = -time * x;
				dSaug(n+1, i)         = z.real();
				dSaug(n+1+nActive, i) = z.imag();
			} 
		}
		// beta adjustment parameter
		else if( 2*M + 1 == q ) 
		{
			// applies to all columns (of basis matrix)
			for(integer i = 1; i <= Q; i++ ) 
				//for(integer i = 1; i < Q; i++ ) 
			{
				for( integer n = 0; n < nActive; n++ ) 
				{
					treal time = params.m_t(n+nStart);

					tcomplex x = tcomplex(params.m_Sp(n+1,i), params.m_Sp(n+1+nActive,i));
					tcomplex z = -(time * time) * x;

					dSaug(n+1, i)         = z.real();
					dSaug(n+1+nActive, i) = z.imag();
				} 
			}
		}
		// phi0 adjustment parameter
		else if( 2*M + 2 == q ) 
		{
			// applies to all columns (of basis matrix)
			for(integer i = 1; i <= Q; i++ ) 
				//for(integer i = 1; i < Q; i++ ) 
			{
				for( integer n = 0; n < nActive; n++ ) 
				{
					tcomplex x = tcomplex(params.m_Sp(n+1,i), params.m_Sp(n+1+nActive,i));
					tcomplex z = tcomplex(0,1) * x;

					dSaug(n+1, i)         = z.real();
					dSaug(n+1+nActive, i) = z.imag();
				} 
			}
		}
		// phi1 adjustment parameter
		else if( 2*M + 3 == q ) 
		{
			// applies to all columns (we applied it to y)
			for(integer i = 1; i <= Q; i++ ) 
				//for(integer i = 1; i < Q; i++ ) 
			{
				for( integer n = 0; n < nActive; n++ ) 
				{
					tcomplex z = tcomplex(params.m_Spp(n+1,i), params.m_Spp(n+1+nActive,i));

					dSaug(n+1, i)         = z.real();
					dSaug(n+1+nActive, i) = z.imag();
				} 
			}
		}

		// construct dSaugF (i.e. dSaug(:,F), storage is external)
		/*cvm::rmatrix dSaugF(params.m_dSaugF, 1, 1, 2*nActive, F.size());
		  c = 1;
		  for( std::set<integer>::iterator itf = F.begin(); itf != F.end(); itf++, c++ ) 
		  dSaugF(c) = dSaug(*itf);
		  */

		// the storage for this matrix is external
		//cvm::rmatrix dX(params.m_dX, 1, 1, F.size(), 2*nActive);

		// derivative of a with respect to theta_q
		cvm::rvector da(F.size());

		// derivative of pseudo-inverse of NNLS projection matrix (wow!)
		//dX = -(W*( (~dSaugF)*SaugF  + (~SaugF)*dSaugF)*W)*(~SaugF) + W*(~dSaugF);
		//da = dX*params.m_yActive;

		// the q'th column of the Jacobian as defined on paper == q'th row of levmar Jacobian
		J[q] = dSaug*annls; // + SaugF*da;
	}

	params.m_final_jac = ~J;
}

bool tarquin::RunTARQUIN(Workspace& work, CBoswell& log)
{
	try
	{
		// note that we are fitting to the processed FID, if available, else the raw one
		CFID& fid               = work.GetFID();
		CBasis& basis           = work.GetBasis();
		Options& options        = work.GetOptions();

		log.LogMessage(LOG_INFO, "Running optimiser");

		// perform a dummy initialisation to get the matrix sizes
		// for the group and summation matrices
		basis.Initialise(4.65f);
		if ( options.GetBasisComp() )
			basis.CompressGroupMatrix();
		options.Initialise(basis, fid, log);

		// the original group matrix
		const cvm::cmatrix& G = basis.GetGroupMatrix();

		// the summation matrix, S = A*G 
		const cvm::rmatrix A = basis.GetSummationMatrix().real();

		// number of samples
		integer N = G.msize();

		// number of groups vectors 
		integer M = G.nsize();

		// number of parameters (2 per group (delta f and delta alpha), commom: beta, phi0, phi1) 
		integer P = M*2 + 1 + 1 + 1;

		// number of basis vectors used for NNLS solution
		integer Q = basis.GetBasisMatrix().nsize();

		log.LogMessage(LOG_INFO, "%u group vectors in total", M);
		log.LogMessage(LOG_INFO, "%u basis vectors in total", Q);

		// lower constraint for parameters
		cvm::rvector vLowerBounds(P);

		// upper constraints for parameters
		cvm::rvector vUpperBounds(P);

		// parameters (to be determined)
		cvm::rvector vParams(P);

		// names of metabolites for identifying when
		// soft contraints can be used
		std::vector<std::string> metab_names = basis.GetSignalNames();

		// set the constraints corresponding to each column
		// format is:
		// [ frequency limits; alpha limits; beta limit; phi0; phi1 ]
		for( integer m = 1; m <= M; m++ ) {
			//std::cout << options.GetLowerLimitShift(m-1) << std::endl;
			vLowerBounds(m) = options.GetLowerLimitShift(m-1); 
			vUpperBounds(m) = options.GetUpperLimitShift(m-1); 
			vParams(m)      = options.GetTypicalShift(m-1);

			vLowerBounds(M + m) = options.GetLowerLimitAlpha(m-1); 
			vUpperBounds(M + m) = options.GetUpperLimitAlpha(m-1); 
			vParams(M + m)      = options.GetTypicalAlpha(m-1);
		}

		// limits for beta
		std::size_t nIdxBeta   = 2*M+1;
		vLowerBounds(nIdxBeta) = options.GetLowerLimitBeta(0); 
		vUpperBounds(nIdxBeta) = options.GetUpperLimitBeta(0);
		vParams(nIdxBeta)      = options.GetTypicalBeta(0); 

		// limits for phi0
		vLowerBounds(2*M + 2) = options.GetLowerLimitPhi0(); 
		vUpperBounds(2*M + 2) = options.GetUpperLimitPhi0(); 
		vParams(2*M + 2)      = options.GetTypicalPhi0(); 

		// limits for phi1
		vLowerBounds(2*M + 3) = options.GetLowerLimitPhi1(); 
		vUpperBounds(2*M + 3) = options.GetUpperLimitPhi1(); 
		vParams(2*M + 3)      = options.GetTypicalPhi1(); 

		// the sampling interval (time step)
		treal dt = 1.0 / fid.GetSamplingFrequency(); 

		// generate time signal
		cvm::rvector t(N);
		for(integer n = 0; n < N; n++)
			t(n+1) = n*dt;

		// the range over which we are fitting
		integer nStart = options.GetRangeStart();
		integer nEnd   = options.GetRangeEnd();

		log.LogMessage(LOG_INFO, "Final nStart = %i",options.GetRangeStart());
		log.LogMessage(LOG_INFO, "Final nEnd   = %i",options.GetRangeEnd());

		// the number of (complex) points we will use in the fit
		integer activeN = 1 + nEnd - nStart;

		// this is the frequency range over which to apply phi1
		cvm::rvector freq_range = fid.GetFreqScale();
		cvm::rvector shift_freq_range = fftshift(freq_range);

		const coord_vec& fit_list = options.GetFitList();
		int voxel_num = 0;

		log.LogMessage(LOG_INFO, "%u voxel(s) to be fitted", fit_list.size());
        
        // check for overlapping signals
        bool TNAA = false;
        int TNAA_ind = -1;
        int NAA_pos = -1;
        int NAAG_pos = -1;

        bool TCho = false;
        int TCho_ind = -1;
        int PCh_pos = -1;
        int GPC_pos = -1;

        bool TCr = false;
        int TCr_ind = -1;
        int Cr_pos = -1;
        int PCr_pos = -1;

        bool Glx = false;
        int Glx_ind = -1;
        int Glu_pos = -1;
        int Gln_pos = -1;

        bool TLM09 = false;
        int TLM09_ind = -1;
        int Lip09_pos = -1;
        int MM09_pos = -1;
        
        /*bool TL13 = false;
        int TL13_ind = -1;
        int Lip13a_pos = -1;
        int Lip13b_pos = -1;
        this would require a seperate CRLB calc as it overlapps with TLM13 */
        
        bool TLM13 = false;
        int TLM13_ind = -1;
        int Lip13a_pos = -1;
        int Lip13b_pos = -1;
        int MM12_pos = -1;
        int MM14_pos = -1;
        
        bool TLM20 = false;
        int TLM20_ind = -1;
        int Lip20_pos = -1;
        int MM20_pos = -1;

        bool TGABA = false;
        int TGABA_ind = -1;
        int GABAA_pos = -1;
        int GABAB_pos = -1;

        //bool TLM13 = false;

        std::vector<int> pos_list;
        int overlapping_signals = 0;

        int extra_cols = 0;

        bool NAA_found = false;
        if (std::find(metab_names.begin(), metab_names.end(), "NAA") != metab_names.end())
        {
            NAA_found = true;
            NAA_pos = 1 + std::find(metab_names.begin(), metab_names.end(), "NAA") - metab_names.begin();
        }
        bool NAAG_found = false;
        if (std::find(metab_names.begin(), metab_names.end(), "NAAG") != metab_names.end())
        {
            NAAG_found = true;
            NAAG_pos = 1 + std::find(metab_names.begin(), metab_names.end(), "NAAG") - metab_names.begin();
        }

        bool PCh_found = false;
        if (std::find(metab_names.begin(), metab_names.end(), "PCh") != metab_names.end())
        {
            PCh_found = true;
            PCh_pos = 1 + std::find(metab_names.begin(), metab_names.end(), "PCh") - metab_names.begin();
        }
        bool GPC_found = false;
        if (std::find(metab_names.begin(), metab_names.end(), "GPC") != metab_names.end())
        {
            GPC_found = true;
            GPC_pos = 1 + std::find(metab_names.begin(), metab_names.end(), "GPC") - metab_names.begin();
        }

        bool Cr_found = false;
        if (std::find(metab_names.begin(), metab_names.end(), "Cr") != metab_names.end())
        {
            Cr_found = true;
            Cr_pos = 1 + std::find(metab_names.begin(), metab_names.end(), "Cr") - metab_names.begin();
        }
        bool PCr_found = false;
        if (std::find(metab_names.begin(), metab_names.end(), "PCr") != metab_names.end())
        {
            PCr_found = true;
            PCr_pos = 1 + std::find(metab_names.begin(), metab_names.end(), "PCr") - metab_names.begin();
        }

        bool Glu_found = false;
        if (std::find(metab_names.begin(), metab_names.end(), "Glu") != metab_names.end())
        {
            Glu_found = true;
            Glu_pos = 1 + std::find(metab_names.begin(), metab_names.end(), "Glu") - metab_names.begin();
        }
        bool Gln_found = false;
        if (std::find(metab_names.begin(), metab_names.end(), "Gln") != metab_names.end())
        {
            Gln_found = true;
            Gln_pos = 1 + std::find(metab_names.begin(), metab_names.end(), "Gln") - metab_names.begin();
        }
        bool Lip09_found = false;
        if (std::find(metab_names.begin(), metab_names.end(), "Lip09") != metab_names.end())
        {
            Lip09_found = true;
            Lip09_pos = 1 + std::find(metab_names.begin(), metab_names.end(), "Lip09") - metab_names.begin();
        }
        bool MM09_found = false;
        if (std::find(metab_names.begin(), metab_names.end(), "MM09") != metab_names.end())
        {
            MM09_found = true;
            MM09_pos = 1 + std::find(metab_names.begin(), metab_names.end(), "MM09") - metab_names.begin();
        }
        bool Lip13a_found = false;
        if (std::find(metab_names.begin(), metab_names.end(), "Lip13a") != metab_names.end())
        {
            Lip13a_found = true;
            Lip13a_pos = 1 + std::find(metab_names.begin(), metab_names.end(), "Lip13a") - metab_names.begin();
        }
        bool Lip13b_found = false;
        if (std::find(metab_names.begin(), metab_names.end(), "Lip13b") != metab_names.end())
        {
            Lip13b_found = true;
            Lip13b_pos = 1 + std::find(metab_names.begin(), metab_names.end(), "Lip13b") - metab_names.begin();
        }
        bool MM12_found = false;
        if (std::find(metab_names.begin(), metab_names.end(), "MM12") != metab_names.end())
        {
            MM12_found = true;
            MM12_pos = 1 + std::find(metab_names.begin(), metab_names.end(), "MM12") - metab_names.begin();
        }
        bool MM14_found = false;
        if (std::find(metab_names.begin(), metab_names.end(), "MM14") != metab_names.end())
        {
            MM14_found = true;
            MM14_pos = 1 + std::find(metab_names.begin(), metab_names.end(), "MM14") - metab_names.begin();
        }
        bool Lip20_found = false;
        if (std::find(metab_names.begin(), metab_names.end(), "Lip20") != metab_names.end())
        {
            Lip20_found = true;
            Lip20_pos = 1 + std::find(metab_names.begin(), metab_names.end(), "Lip20") - metab_names.begin();
        }
        bool MM20_found = false;
        if (std::find(metab_names.begin(), metab_names.end(), "MM20") != metab_names.end())
        {
            MM20_found = true;
            MM20_pos = 1 + std::find(metab_names.begin(), metab_names.end(), "MM20") - metab_names.begin();
        }

        bool GABAA_found = false;
        if (std::find(metab_names.begin(), metab_names.end(), "GABA_A") != metab_names.end())
        {
            GABAA_found = true;
            GABAA_pos = 1 + std::find(metab_names.begin(), metab_names.end(), "GABA_A") - metab_names.begin();
        }
        bool GABAB_found = false;
        if (std::find(metab_names.begin(), metab_names.end(), "GABA_B") != metab_names.end())
        {
            GABAB_found = true;
            GABAB_pos = 1 + std::find(metab_names.begin(), metab_names.end(), "GABA_B") - metab_names.begin();
        }

        std::vector<std::string> &metab_names_comb = work.GetMetabNamesComb();

        int cnt = 1;
        if ( NAA_found && NAAG_found )
        {
            TNAA = true;
            overlapping_signals = overlapping_signals + 2;
            extra_cols = extra_cols + 1;
            pos_list.push_back(NAA_pos);
            pos_list.push_back(NAAG_pos);
            TNAA_ind = cnt;
            cnt++;
            metab_names_comb.push_back("TNAA");
        }
        if ( PCh_found && GPC_found )
        {
            TCho = true;
            overlapping_signals = overlapping_signals + 2;
            extra_cols = extra_cols + 1;
            pos_list.push_back(PCh_pos);
            pos_list.push_back(GPC_pos);
            TCho_ind = cnt;
            cnt++;
            metab_names_comb.push_back("TCho");
        }
        if ( Cr_found && PCr_found )
        {
            TCr = true;
            overlapping_signals = overlapping_signals + 2;
            extra_cols = extra_cols + 1;
            pos_list.push_back(Cr_pos);
            pos_list.push_back(PCr_pos);
            TCr_ind = cnt;
            cnt++;
            metab_names_comb.push_back("TCr");
        }
        if ( Glu_found && Gln_found )
        {
            Glx = true;
            overlapping_signals = overlapping_signals + 2;
            extra_cols = extra_cols + 1;
            pos_list.push_back(Glu_pos);
            pos_list.push_back(Gln_pos);
            Glx_ind = cnt;
            cnt++;
            metab_names_comb.push_back("Glx");
        }
        if ( Lip09_found && MM09_found )
        {
            TLM09 = true;
            overlapping_signals = overlapping_signals + 2;
            extra_cols = extra_cols + 1;
            pos_list.push_back(Lip09_pos);
            pos_list.push_back(MM09_pos);
            TLM09_ind = cnt;
            cnt++;
            metab_names_comb.push_back("TLM09");
        }
        //if ( Lip13a_found && Lip13b_found )
        //{
        //    TL13 = true;
        //    overlapping_signals = overlapping_signals + 2;
        //    extra_cols = extra_cols + 1;
        //    pos_list.push_back(Lip13a_pos);
        //    pos_list.push_back(Lip13b_pos);
        //    TL13_ind = cnt;
        //    cnt++;
        //    metab_names_comb.push_back("TL13");
        //}
        if ( Lip13a_found && Lip13b_found && MM12_found && MM14_found )
        {
            TLM13 = true;
            overlapping_signals = overlapping_signals + 4;
            extra_cols = extra_cols + 1;
            pos_list.push_back(Lip13a_pos);
            pos_list.push_back(Lip13b_pos);
            pos_list.push_back(MM12_pos);
            pos_list.push_back(MM14_pos);
            TLM13_ind = cnt;
            cnt++;
            metab_names_comb.push_back("TLM13");
        }
        if ( Lip20_found && MM20_found )
        {
            TLM20 = true;
            overlapping_signals = overlapping_signals + 2;
            extra_cols = extra_cols + 1;
            pos_list.push_back(Lip20_pos);
            pos_list.push_back(MM20_pos);
            TLM20_ind = cnt;
            cnt++;
            metab_names_comb.push_back("TLM20");
        }
        if ( GABAA_found && GABAB_found )
        {
            TGABA = true;
            overlapping_signals = overlapping_signals + 2;
            extra_cols = extra_cols + 1;
            pos_list.push_back(GABAA_pos);
            pos_list.push_back(GABAB_pos);
            TGABA_ind = cnt;
            cnt++;
            metab_names_comb.push_back("TGABA");
        }

        log.DebugMessage(DEBUG_LEVEL_1, "Overlapping signals found: %i", overlapping_signals);

		// iterate over the fit list	
		for( coord_vec::const_iterator fit_it = fit_list.begin(); fit_it != fit_list.end(); ++fit_it )	
		{

			// reset all the parameters for CSI
			for( integer m = 1; m <= M; m++ ) 
			{
				vParams(m)          = options.GetTypicalShift(m-1);
				vParams(M + m)      = options.GetTypicalAlpha(m-1);
			}
			vParams(nIdxBeta)       = options.GetTypicalBeta(0); 
			vParams(2*M + 2)        = options.GetTypicalPhi0(); 
			vParams(2*M + 3)        = options.GetTypicalPhi1(); 

			voxel_num++;
			log.LogMessage(LOG_INFO, "Fitting fid %u of %u", voxel_num, fit_list.size());


			log.LogMessage(LOG_INFO, "Adjusting basis set to %f ppm", fid.GetPPMRef(*fit_it));
			// shift the basis to match ref
			basis.Initialise(fid.GetPPMRef(*fit_it));
			if ( options.GetBasisComp() )
				basis.CompressGroupMatrix();

			// the fid we are fitting to
			cvm::cvector& y = fid.GetVectorFID(*fit_it);

			// turn FID into [real; imaginary]
			cvm::rvector yActive( 2*activeN );

			for(integer n = 0; n < activeN; n++) 
			{
				yActive(n+1)         = y(nStart + n).real();
				yActive(n+1+activeN) = y(nStart + n).imag();
			}


			// estimate the standard deviation of the noise in the time domain
			// Greg, I moved this to before fitting in case a big phi1
			// correction causes a distortion at the end of the fid

            int td_block_size = 100;
            double td_noise_min = std::numeric_limits<double>::infinity();
		    double td_noise_temp = 0;
            for ( int n = 1; n < (y.size()+1) - td_block_size; n = n + td_block_size ) 
            {
                td_noise_temp = stdev(y.real(),n,n+td_block_size-1);
                if ( td_noise_temp < td_noise_min )
                    td_noise_min = td_noise_temp;
            }


            // Save noise estimate
			std::vector<double>& TdNoise = work.GetTdNoise();
            TdNoise.push_back(td_noise_min);
			log.LogMessage(LOG_INFO, "Time domain noise  = %f", td_noise_min);

			//plot(y);

			// amplitude estimate (we need to size it since it is intialised here)
			rvec_stdvec& ahat_vec = work.GetAmplitudes();
			cvm::rvector ahat;
			ahat.resize(Q);

			// precompute the DFT of each column of the group matrix
			cvm::cmatrix GDFT(G.msize(), G.nsize());

			for( int c = 1; c <= G.nsize(); ++c )
			{
				cvm::cvector GDFTc = GDFT(c);

				fft(G(c), GDFTc);
			}

			treal lambda = options.GetLambda();
			// data structure of variables used by objective function
			SResParams params(y, basis, yActive, G, GDFT, A, t, ahat, nStart, shift_freq_range, 
					activeN, M, Q, metab_names, lambda, log, options);

			// before running the optimiser, we do an initial projection, to test whether the data
			// is worth running through the optimiser, i.e. we check that there is a solution to the NNLS 
			// problem such that some of the elements are non-zero (and positive)

			//CNNLSSolverLH nnls;
            //NNLSSolverBJ nnls;
            FNNLSSolverBJ nnls;

			for( integer c = 1; c <= M; c++ ) 
			{
				for( integer n = 0; n < activeN; n++ ) 
				{
					tcomplex z = G(n+nStart, c);

					// modify sample and put in right place in active matrix
					params.m_Gp(n+1, c)         = z.real();
					params.m_Gp(n+1+activeN, c) = z.imag();
				}
			}

			params.m_Sp = params.m_Gp * params.m_A;

			//std::cout << params.m_Sp << std::endl;
			//std::cout << params.m_yActive << std::endl;
			//std::cout << params.m_ahat << std::endl;

			//plot(params.m_Sp);
			//plot(params.m_yActive);
			//plot(params.m_ahat);

			// estimate amplitudes
			nnls.solve(params.m_Sp, params.m_yActive, params.m_ahat, false);

			//std::cout << params.m_ahat << std::endl;

			// find number of elements which are greater than zero
			std::set<integer> F;
			for( integer n = 1; n <= Q; n++ )
				if( params.m_ahat(n) > 0.0 ) 
					F.insert(n);

			// the data was so bad that no solution can be found
			if( 0 == F.size() ) 
			{
				log.LogMessage(LOG_ERROR, "Warning. Data quality is very poor.");
				//log.LogMessage(LOG_ERROR, "FAILED. Data quality is very poor, not fitting.");
				//return false;
			}

			// now an initial sweep of the beta value so that we get a sensible starting position
			//log.BeginTask("Estimating initial Gaussian damping term...");
			const std::size_t numBetaSteps = 200;

			treal betaStepSize = ( pow(options.GetUpperLimitBeta(0) - options.GetLowerLimitBeta(0), 0.5) ) / numBetaSteps;
			treal bestResidual = std::numeric_limits<treal>::infinity();
			treal bestBeta     = 0;

			cvm::rvector yhatBeta(yActive.size());
			cvm::rvector yhatBeta_orig(yActive.size());
			cvm::rvector residualt(yActive.size());

			// generate the basis matrix
			params.m_Sp = params.m_Gp * params.m_A;

			// construct yhat
			yhatBeta_orig = params.m_Sp * params.m_ahat;

			// for each step 
			for( std::size_t step = 0; step < numBetaSteps; ++step ) {

				treal beta = options.GetLowerLimitBeta(0) + pow((step * betaStepSize) , 2);

				for( integer n = 0; n < activeN; n++ ) {

					treal z = exp(-beta * t(n+nStart) * t(n+nStart));

					// modify sample and put in right place in active matrix
					yhatBeta(n+1)         = yhatBeta_orig(n+1) * z;
					yhatBeta(n+1+activeN) = yhatBeta_orig(n+1+activeN) * z;
				}

				residualt = yActive - yhatBeta;

				//std::cout << "\nTrying beta = " << beta << " residual = " << residualt.norm2() << std::flush;

				if( residualt.norm2() < bestResidual ) {
					bestResidual = residualt.norm2();
					bestBeta = beta;
				}
			}

			//log.DebugMessage(DEBUG_LEVEL_1, "Initial beta estimate: %2.f", bestBeta);
			vParams[nIdxBeta] = pow(options.GetInitBetaUsed(voxel_num - 1),1) / options.GetBetaScale(); // MAGIC BETA SCALE
			log.DebugMessage(DEBUG_LEVEL_1, "Using init_beta: %f", pow(vParams[nIdxBeta],1)*options.GetBetaScale()); // MAGIC BETA SCALE

			// set the maximum number of iterations	
			integer iterations = options.GetMaxIters();
			log.DebugMessage(DEBUG_LEVEL_1, "Doing a maximum of %d iterations.", iterations);

			// optimiser options
			double opts[LM_OPTS_SZ];

			opts[0] = options.GetNUMERICAL_TOL_INIT_MU();
			opts[1] = NUMERICAL_TOL_INF_JE;
			opts[2] = options.GetNUMERICAL_TOL_L2_DP();
			opts[3] = NUMERICAL_TOL_L2_RES;
			opts[4] = NUMERICAL_TOL_DIFF_DELTA;

			// convert opts to stl vector and add to workspace for serialisation
			std::vector<double> STL_opts;
			for (int n = 0; n < LM_OPTS_SZ; n++)
				STL_opts.push_back(opts[n]);
			work.SetLMopts(STL_opts);

			// info about result
			double info[LM_INFO_SZ];

			// works well for 3t case 
			//treal norm = yActive.norm1()*10000;

			treal norm = yActive.norm1()/100;
			yActive = yActive / norm;
			ahat = ahat / norm;
			// call the optimiser analytical Jacobian 
			dlevmar_bc_der(residual_objective_all, jacobian_func, 
					vParams, yActive, P, 2*activeN, vLowerBounds, vUpperBounds, 
					iterations, opts, info, NULL, NULL, (void*)&params);

			// call the optimiser numerical Jacobian 
			//dlevmar_bc_dif(residual_objective_all, vParams, yActive, P, 2*activeN, vLowerBounds, vUpperBounds, iterations, opts, info, NULL, NULL, (void*)&params);

			yActive = yActive * norm;

			// useful for debugging
			// ahat.set(1);

			ahat = ahat * norm;

			// TODO: write a function for outputting optimiser convergence results

			//
			// Populate the workspace with all these nice results.
			// 

			// convert LM convergance info to stl vector and add to workspace for serialisation
			std::vector<double> STL_info;
			for (int n = 0; n < LM_INFO_SZ; n++)
				STL_info.push_back(info[n]);

			std::vector<std::vector<double> >& STL_info_vec = work.GetLMinfo();
			STL_info_vec.push_back(STL_info);

            vParams(nIdxBeta) = pow(vParams(nIdxBeta),1) * options.GetBetaScale(); // MAGIC BETA SCALE

			work.AppendParas(vParams);

			cmat_stdvec& GpOut_vec = work.GetGroupMatrix();
			cvm::cmatrix GpOut;
			cmat_stdvec& SpOut_vec = work.GetBasisMatrix();
			cvm::cmatrix SpOut;
			cvec_stdvec& yhat_vec  = work.GetSignalEstimate();
			cvm::cvector yhat;

			GpOut.resize(N, M);
			SpOut.resize(N, Q);
			yhat.resize(N);


			//std::cout << vParams(1) << std::endl;
			//
			// for each column of the group matrix
			for( integer c = 1; c <= M; c++ ) 
			{
				// indices to parameter vector	
				std::size_t nIdxShift = c;
				std::size_t nIdxAlpha = c+M;

				// the parameters to apply
				tcomplex jomega(0, 2.0 * M_PI * vParams(nIdxShift));
				treal alpha = vParams(nIdxAlpha);
				treal beta  = vParams(nIdxBeta);

				// for each (sample of the active part of the signal == row of G)
				for( integer n = 1; n <= N; n++ ) 
				{
					// the modified sample of the group matrix 
					tcomplex z = G(n, c) * exp(-t(n)*alpha -t(n)*t(n)*beta + t(n)*jomega);

					// modify sample and put in right place in active matrix
					GpOut(n, c) = z;
				}
			}

			std::size_t nIdxPhi0  = 2*M+2; 
			std::size_t nIdxPhi1  = 2*M+3;
			// useful debugging output, what happens for multi-voxel?
			/*
			   std::ofstream fitted_paras;
			   fitted_paras.open ("fitted_paras.csv");
			   fitted_paras << "name," << "group," << "shift (Hz)," << "damping" << std::endl;
			   int basis_cnt = 1;
			   int group_cnt = 0;
			   for( integer c = 1; c <= M; c++ ) 
			   {
			   group_cnt++;
			   std::size_t nIdxShift = c;
			   std::size_t nIdxAlpha = c+M;
			// some output
			integer i = params.m_basis.GetBasisFromGroup(c-1);
			if ( i != basis_cnt )
			{
			basis_cnt = i;
			group_cnt = 1;
			}
			fitted_paras << params.m_metab_names[i-1] << "," << group_cnt << "," << vParams(nIdxShift) << "," << vParams(nIdxAlpha) << std::endl;
			}
			//std::size_t nIdxBeta  = 2*M+1;
			fitted_paras << "beta," << vParams(nIdxBeta) << std::endl;
			fitted_paras << "phi0," << vParams(nIdxPhi0) << std::endl;
			fitted_paras << "phi1," << vParams(nIdxPhi1) << std::endl;
			fitted_paras.close();
			*/

			//
			// Synthesise, fit over whole signal, note: amplitudes come from NNLS 
			// solution over [nStart, nEnd] only.
			//

			// get matrix of metabolites by summing groups
			SpOut = GpOut * basis.GetSummationMatrix();

			// synthesize estimate
			yhat = SpOut * cvm::cvector(ahat);

            cvm::rvector ahat_singlet(ahat.size());
            
            bool singlet_found = false;
            
            if ( NAA_found )
            {
                ahat_singlet(NAA_pos) = ahat(NAA_pos);
                singlet_found = true;
            }

            if ( NAAG_found )
            {
                ahat_singlet(NAAG_pos) = ahat(NAAG_pos);
                singlet_found = true;
            }

            if ( Cr_found )
            {
                ahat_singlet(Cr_pos) = ahat(Cr_pos);
                singlet_found = true;
            }
            
            if ( PCr_found )
            {
                ahat_singlet(PCr_pos) = ahat(PCr_pos);
                singlet_found = true;
            }

            if ( GPC_found )
            {
                ahat_singlet(GPC_pos) = ahat(GPC_pos);
                singlet_found = true;
            }

            if ( PCh_found )
            {
                ahat_singlet(PCh_pos) = ahat(PCh_pos);
                singlet_found = true;
            }

            // find yhat for singlets only TODO add Cho to list
            cvm::cvector yhat_singlet = SpOut * cvm::cvector(ahat_singlet);
            
            cvm::rvector ahat_metab(ahat.size());
            cvm::rvector ahat_broad(ahat.size());
            const std::vector<bool> broad_vec = basis.GetBroadSig();

		    std::vector<std::string> metab_names = basis.GetSignalNames();
		    for( int n = 0; n < broad_vec.size(); n++ )	
            {
                if ( !broad_vec[n] )
                {
			        if ( ( metab_names[n] != "Lac" ) && ( metab_names[n] != "Ala" ) )
                    {
                        ahat_metab(n+1) = ahat(n+1);
                    }
                }
                else
                    ahat_broad(n+1) = ahat(n+1);
            }


            // find yhat for metabs only
            cvm::cvector yhat_metab = SpOut * cvm::cvector(ahat_metab);
            cvm::cvector yhat_broad = SpOut * cvm::cvector(ahat_broad);
            
			//
			// Apply phasing parameters to signal we fitted to for output.
			//
			//std::size_t nIdxPhi0  = 2*M+2; 
			//std::size_t nIdxPhi1  = 2*M+3;

			treal phi0 = vParams(nIdxPhi0); 
			treal phi1 = vParams(nIdxPhi1);

            // update phi0 and phi1 in CFID
            treal phi0_old = fid.GetPhi0(*fit_it);
            treal phi1_old = fid.GetPhi1(*fit_it);
            fid.SetPhi0(*fit_it,-phi0+phi0_old);
            fid.SetPhi1(*fit_it,-phi1+phi1_old);
            

			cvm::cvector Y(y.size());
			fft(y, Y);

			for( integer n = 0; n < Y.size(); n++ ) 
			{
				treal freq = shift_freq_range[n+1];
				Y[n+1] = Y[n+1] * exp( tcomplex(0, -phi0 -phi1*2.0*M_PI*freq) );
			}

			ifft(Y, y);

            // replace first part of FID with first part of y
            if ( options.GetReplaceFp() )
            {
                for (int n = 0; n < options.GetRangeStart(); n++ )
                    y[n+1] = yhat[n+1];
            }

			log.EndTask("done.");

			log.DebugMessage(DEBUG_LEVEL_1, "Final Beta: %.2f", vParams(nIdxBeta));
			log.DebugMessage(DEBUG_LEVEL_1, "Final Phi0: %.2f", vParams(nIdxPhi0));
			log.DebugMessage(DEBUG_LEVEL_1, "Final Phi1: %.2f", vParams(nIdxPhi1));


			// measure of fit quality
			cvm::cvector yz = y; 
			cvm::cvector yhatz = yhat; 
			cvm::cvector yhatz_singlet = yhat_singlet;
			cvm::cvector yhatz_metab = yhat_metab;
			cvm::cvector yhatz_broad = yhat_broad;

			//int zf = 1;
			// zero fill if less than 4096 points
			//if ( y.size() < 4096 )
			//	zf = round(4096.0f/ float(y.size()) );

            int zf = options.GetZF();

			//zf = 1;

			yz.resize(yz.size()*zf);
			yhatz.resize(yhatz.size()*zf);
			yhatz_singlet.resize(yhatz_singlet.size()*zf*2);
			yhatz_metab.resize(yhatz_metab.size()*zf);
			yhatz_broad.resize(yhatz_broad.size()*zf);

			// copy last pts points of real fid to end of zfilled fid
			int pts = 5;
			if ( zf > 1 )
			{
				for ( int n = 1 ; n < pts + 1; n++)
				{
					yz(yz.size()-pts+n) = yz(yz.size()/zf-pts+n);
					yz(yz.size()/zf-pts+n) = 0;
				}
			}

			Y.resize(yz.size());
			cvm::cvector YHAT(yhatz.size());
			cvm::cvector YHAT_singlet(yhatz_singlet.size());
			cvm::cvector YHAT_metab(yhatz_metab.size());
			cvm::cvector YHAT_broad(yhatz_broad.size());

			fft(yz, Y);
			fft(yhatz, YHAT);
			fft(yhatz_singlet, YHAT_singlet);
			fft(yhatz_metab, YHAT_metab);
			fft(yhatz_broad, YHAT_broad);

			Y = fftshift(Y);
			YHAT = fftshift(YHAT);
			YHAT_singlet = fftshift(YHAT_singlet);
			YHAT_metab = fftshift(YHAT_metab);
			YHAT_broad = fftshift(YHAT_broad);


			cvm::cvector residual = yz-yhatz;
			cvm::cvector baseline = residual;

            //plot(residual);

			cvm::cvector RESIDUAL = fft(residual);
			RESIDUAL = fftshift(RESIDUAL);

			//cvm::cvector RES_TEST = fftshift(RESIDUAL);
            //cvm::cvector res_test = ifft(RES_TEST);
            //plot(res_test);
            
            //plot(RESIDUAL);

			cvm::cvector BASELINE;
			td_conv_ws( RESIDUAL, BASELINE, options.GetBL()*options.GetZF(), 10);	

			cvm::rvector freq_scale = fid.GetPPMScale(*fit_it, zf);

			cvm::rvector freq_scale_singlet = fid.GetPPMScale(*fit_it, zf*2);

            // guess metabolite FWHM
            int max_data_pt = YHAT_singlet.real().indofmax();
            double max_val = YHAT_singlet(max_data_pt).real();

            // guess init beta based on the fwhm of prod_shift_cut_metabs
            // find the right base of the peak
            int right_pt = -1;
            if ( singlet_found )
            {
                for ( int n = max_data_pt + 2; n < YHAT_singlet.size(); n++ )
                {
                    if ( YHAT_singlet.real()(n) < (max_val/2.0) )
                    {
                        right_pt = n-1;
                        log.LogMessage(LOG_INFO, "FWHM right ppm = %f", freq_scale_singlet(right_pt));
                        break;
                    }
                }
            }

            int left_pt = -1;
            if ( singlet_found )
            {
                // find the left base of the peak
                for ( int n = max_data_pt - 2; n > 0; n-- )
                {
                    if ( YHAT_singlet.real()(n) < (max_val/2.0) )
                    {
                        left_pt = n+1;
                        log.LogMessage(LOG_INFO, "FWHM left ppm = %f", freq_scale_singlet(left_pt));
                        break;
                    }
                }
            }
            
            if ( left_pt == -1 || right_pt == -1 )
            {
                log.DebugMessage(DEBUG_LEVEL_1, "Warning, metabolite FWHM calc failed, carry on regardless.");
                std::vector<double>& metab_fwhm_vec = work.GetMetabFWHM();
                metab_fwhm_vec.push_back(std::numeric_limits<double>::infinity());
            }
            else
            {
                double metab_fwhm = freq_scale_singlet(left_pt)-freq_scale_singlet(right_pt);
                std::vector<double>& metab_fwhm_vec = work.GetMetabFWHM();
                metab_fwhm_vec.push_back(metab_fwhm);
                log.LogMessage(LOG_INFO, "Metabolite FWHM (ppm) = %f", metab_fwhm);
            }

            //plot(freq_scale_singlet,YHAT_singlet);
            //plot(freq_scale_singlet,YHAT);
            //plot(freq_scale_singlet,YHAT_metab);

            //std::cout << std::endl << "HI: " << options.GetPPMend() << std::endl;

			// find points corresponding to ppm start and ppm end
			int left = 1, right = freq_scale.size()-1;
			for ( int n = 1; n < (freq_scale.size()); n++ ) 
			{
				if ( ( freq_scale(n) > options.GetPPMend() ) && ( freq_scale(n+1) <= options.GetPPMend() ) )
					left = n;
				if ( ( freq_scale(n) > options.GetPPMstart() ) && ( freq_scale(n+1) <= options.GetPPMstart() ) )
					right = n;
			}

            if ( left == 1 || right == (freq_scale.size() - 1) )
                log.LogMessage(LOG_INFO, "Warning, one of the PPM limits is outside the bandwidth.");

            // find points corresponding to metab start end
			int metab_left = 1, metab_right = freq_scale.size()-1;
			for ( int n = 1; n < (freq_scale.size()); n++ ) 
			{
				if ( ( freq_scale(n) > 4.0 ) && ( freq_scale(n+1) <= 4.0 ) )
					metab_left = n;
				if ( ( freq_scale(n) > 1.9 ) && ( freq_scale(n+1) <= 1.9 ) )
					metab_right = n;
			}
            
            if ( metab_left == 1 || metab_right == (freq_scale.size() - 1) )
                log.LogMessage(LOG_INFO, "Warning, one of the metab PPM limits is outside the bandwidth.");

            // find max point  in metab region
            double Ymax_metab = 0;
			for ( int n = metab_left; n < metab_right+1; n++ ) 
				if ( Ymax_metab < std::abs ( Y(n).real() - BASELINE(n).real() ) )
					Ymax_metab = std::abs ( Y(n).real() - BASELINE(n).real() );

			double Ymax = 0;
			cvm::rvector BASELINE_REAL(right-left+1);
			cvm::rvector RESIDUAL_REAL(right-left+1);
			cvm::rvector freq_scale_cut(right-left+1);
			// find max of y for SNR calculation
			// should it be "- BASELINE(n).real()" ?
			for ( int n = left; n < right+1; n++ ) 
			{
				//std::cout << std::abs( Y(n).real() - BASELINE(n).real() ) << std::endl;
				if ( Ymax < std::abs ( Y(n).real() - BASELINE(n).real() ) )
					Ymax = std::abs ( Y(n).real() - BASELINE(n).real() ) ;

                BASELINE_REAL(n-left+1) = BASELINE(n).real();
                RESIDUAL_REAL(n-left+1) = RESIDUAL(n).real();
                freq_scale_cut(n-left+1) = freq_scale(n);
			}
            
			cvm::rvector freq_scale_cut_cut(right-left);
			for ( int n = 1; n < right-left; n++ ) 
                freq_scale_cut_cut(n) = freq_scale_cut(n);

			double spec_noise_min = std::numeric_limits<double>::infinity();
			double spec_noise_temp = 0;
			int block_size = 100;
			for ( int n = 1; n < (RESIDUAL.size()+1) - block_size; n = n + block_size ) 
			{
				spec_noise_temp = stdev(RESIDUAL.real(),n,n+block_size-1);
				if ( spec_noise_temp < spec_noise_min )
					spec_noise_min = spec_noise_temp;
			}
            
            spec_noise_min = spec_noise_min / pow(y.size(),0.5);

			std::vector<double>& SpecNoise = work.GetSpecNoise();
            SpecNoise.push_back(spec_noise_min);
			log.LogMessage(LOG_INFO, "Spec noise            = %f", spec_noise_min);

            //plot(RESIDUAL);

			//    std::cout << std::endl << "stdev of res = " << noise_min << std::endl;
			//   std::cout << "stdev of fit res = " << stdev(RESIDUAL.real() - BASELINE.real(),left,right) << std::endl;

			double Fit_Q = stdev(RESIDUAL.real() - BASELINE.real(),left,right) /spec_noise_min;
			double SNR_res = Ymax/(2*stdev(RESIDUAL.real() - BASELINE.real(),left,right) );
            double max_res = (RESIDUAL_REAL - BASELINE_REAL).norminf();

			//double SNR_true = Ymax/( 2*noise_min );

			log.LogMessage(LOG_INFO, "sdev noise = %f", stdev(RESIDUAL.real() - BASELINE.real(),left,right));

			// Save SNR
			pair_vec& SNR = fid.GetSNR();
			SNR.push_back(std::make_pair(SNR_res, true));

			log.LogMessage(LOG_INFO, "SNR residual = %f", SNR_res);
			//log.LogMessage(LOG_INFO, "SNR true     = %f", SNR_true);
			log.LogMessage(LOG_INFO, "SNR max      = %f", SNR_res*Fit_Q);
            
            double SNR_max_metab = Ymax_metab / ( 2 * spec_noise_min );
            std::vector<double>& MetabSNR_vec = work.GetMetabSNR();
			MetabSNR_vec.push_back(SNR_max_metab);

			log.LogMessage(LOG_INFO, "Ymax       = %f", Ymax);
			log.LogMessage(LOG_INFO, "Ymax metab = %f", Ymax_metab);


			// Save Q
			std::vector<double>& Q_vec = work.GetQ();
			Q_vec.push_back(Fit_Q);
			log.LogMessage(LOG_INFO, "Fit quality = %f", Fit_Q);

			std::vector<double>& Q_rel_vec = work.GetQ_rel();
			Q_rel_vec.push_back(max_res/Ymax*100.0);

            
            double metab_ratio = 100.0*YHAT_metab.norm1()/YHAT.norm1();
            std::vector<double>& metab_rat = work.GetMetabRat();
            metab_rat.push_back(metab_ratio);

            double peak_metab_ratio = 100.0*YHAT_broad.norminf() / YHAT_metab.norminf();
            std::vector<double>& peak_metab_rat = work.GetPeakMetabRat();
            peak_metab_rat.push_back(peak_metab_ratio);

			//log.LogMessage(LOG_INFO, "Metab ratio = %f", metab_ratio);

			cvm::rvector BASELINE_REAL_DIFF;
            diff(BASELINE_REAL, BASELINE_REAL_DIFF);

            
            //double a[] = {1., 2., 3., 4.};
            //double b[] = {1., 2., 3., 4.};
            //double a[] = {4., 3., 2., 1.};
            //double b[] = {4., 3., 2., 1.};
            //const cvm::rvector testx(a, 4); 
            //const cvm::rvector testy(b, 4); 
            
            cvm::rvector mc;
            //lsqf(testx, testy, mc);
            lsqf(freq_scale_cut, BASELINE_REAL, mc);
            //std::cout << std::endl << mc << std::endl;
            
			cvm::rvector line_fit;
            get_fit(freq_scale_cut, line_fit, mc);

            cvm::rvector resids = BASELINE_REAL - line_fit;
            //plot(freq_scale_cut,BASELINE_REAL);
            //plot(freq_scale_cut,resids);

            //double baseline_dev = stdev(BASELINE_REAL,1,BASELINE_REAL.size()) / Ymax;
            double baseline_dev = stdev(resids,1,resids.size()) / Ymax;
			//log.LogMessage(LOG_INFO, "Baseline dev = %f", baseline_dev);
			log.LogMessage(LOG_INFO, "Baseline dev = %f", baseline_dev);
			log.LogMessage(LOG_INFO, "Ymax         = %f", Ymax);
            
            //double baseline_var = BASELINE_REAL_DIFF.norminf();
            //double baseline_var = resids.norm1()/Ymax/resids.size();
			//log.LogMessage(LOG_INFO, "Baseline var = %f", baseline_var);

        	std::vector<double>& BLV_vec = work.GetBLV();
			BLV_vec.push_back(baseline_dev);

            cvm::rvector mc_lims;
            double dlims_freq[] = {freq_scale_cut(1),freq_scale_cut(freq_scale_cut.size())};
            cvm::rvector freq_lims(dlims_freq,2,1);
            double dlims_bl[] = {BASELINE_REAL(1),BASELINE_REAL(BASELINE_REAL.size())};
            cvm::rvector bl_lims(dlims_bl,2,1);
            lsqf(freq_lims, bl_lims, mc_lims);
			cvm::rvector line_fit_lims;
            get_fit(freq_scale_cut, line_fit_lims, mc_lims);
            cvm::rvector BASELINE_LIMS = BASELINE_REAL-line_fit_lims;
            double baseline_shape = 0;
            for ( int n = 1; n < BASELINE_LIMS.size(); n++ )
                baseline_shape += BASELINE_LIMS(n);

            baseline_shape = baseline_shape / BASELINE_REAL.size() / Ymax;
			log.LogMessage(LOG_INFO, "BSL         = %f", baseline_shape);

        	std::vector<double>& BLS_vec = work.GetBLS();
			BLS_vec.push_back(baseline_shape);

            //plot(freq_scale_cut, BASELINE_REAL);
            //plot(freq_scale_cut, BASELINE_LIMS);
            /*plot(freq_scale_cut_cut, resids);
            plot(freq_scale_cut, BASELINE_REAL);
            plot(freq_scale_cut_cut, BASELINE_REAL_DIFF);
            plot(freq_scale_cut_cut, line_fit);*/





			//
			// Attempt to compute the Cramer Rao Lower Bounds of the estimates
			//
			log.BeginTask("Computing CRLBs.");
			cvm::rmatrix JR = params.m_final_jac;

            /*std::cout << std::endl << JR.msize() << std::endl;
            std::cout << JR.nsize() << std::endl;*/


			// Each column of D is the derivative of YHAT w.r.t. each parameter. Since we don't
			// actually compute the derivative of YHAT w.r.t. each amplitude parameter, we just
			// stick the basis vectors as the first few columns of the matrix.
			//cvm::cmatrix D(JR.msize()/2, Q + JR.nsize());
			
            /*cvm::cmatrix D(JR.msize()/2, Q);
			for( int i = 1; i <= activeN; ++i )
				for( int j = 1; j<= Q; ++j )
					D[i][j] = SpOut[i+nStart-1][j];
                    */
            

            // gives identical results to above
			cvm::rmatrix D(JR.msize(), Q+JR.nsize());
			for( int i = 1; i <= activeN; ++i )
            {
                // amps
				for( int j = 1; j<= Q; ++j )
                {
					D[i][j] = SpOut[i+nStart-1][j].real();
					D[i+JR.msize()/2][j] = SpOut[i+nStart-1][j].imag();
                }
            }
                
			for( int i = 1; i <= activeN*2; ++i )
            {
                // other paras
                for( int k = 1; k<= JR.nsize(); ++k )
                {
					D[i][k+Q] = JR[i][k];
                }
            }

            /*std::ofstream dmat;
            dmat.open ("./Dmat.txt");
			for( int i = 1; i <= D.msize(); ++i )
	    		for( int j = 1; j<= D.nsize(); ++j )
                  dmat << D(i,j) << std::endl;

            std::cout << D.msize() << "," << D.nsize() << std::endl;
                  
            //dmat << D(i,j).real() << std::endl << " " << D(i,j).imag() << std::endl;
            
            dmat.close();*/


            // import numpy as np
            // import pylab as pl
            // D = np.matrix(np.reshape(np.loadtxt('Dmat.txt'),(986,73)))
            // IFISHER = np.linalg.pinv(fisher)
            // pl.plot(D[:,60])
            // pl.show()

            /*std::ofstream dmat;
            dmat.open ("./Dmat.txt");
			for( int i = 1; i <= D.msize(); ++i )
	    		for( int j = 1; j<= D.nsize(); ++j )
                  dmat << D(i,j).real() << std::endl << " " << D(i,j).imag() << std::endl;
            
            dmat.close();*/

			// make the fisher information matrix
			cvm::srmatrix FC = (~D)*D;
            
            // estimate noise from TD or FD???
			double sigma = 0;
		    if ( options.GetCRLB_TD() )
            {
                sigma = td_noise_min*td_noise_min;
            }
            else
            {
                sigma = spec_noise_min*spec_noise_min;
            }
            
			//cvm::srmatrix fisher = FC.real();
			cvm::srmatrix fisher = FC;
			fisher *= 1.0 / sigma;
			
			log.DebugMessage(DEBUG_LEVEL_1, "Sigma: %e", sigma);

			cvm::rmatrix IFISHER = fisher.pinv();
			log.EndTask("done.");

			// save in workspace
			rvec_stdvec& crlbs_vec = work.GetCRLBs();
			cvm::rvector crlbs;
			crlbs.resize(Q);

			for( int l = 1; l < Q+1; ++l )
				crlbs[l] = std::sqrt(IFISHER[l][l]);

			crlbs_vec.push_back(crlbs); 

            // Lets find out if we should do an extra CRLB calculation with some
            // signals combined...
            if ( overlapping_signals > 0 )
            {
                cvm::rvector ahat_comb;
                ahat_comb.resize(extra_cols);

                // looks like some overlapping signals were found
                cvm::rmatrix D_comb(JR.msize(), Q-overlapping_signals+extra_cols+JR.nsize());
                // do the non-overlapping signals first
                int k = 1;
                for( int j = 1; j<= Q; ++j )
                {
                    if (std::find(pos_list.begin(), pos_list.end(), j) == pos_list.end()) // if not in pos_list
                    {
                        //std::cout << k << "," << j << std::endl;
                        for( int i = 1; i <= activeN; ++i )
                        {
                            D_comb[i][k] = SpOut[i+nStart-1][j].real();
                            D_comb[i+JR.msize()/2][k] = SpOut[i+nStart-1][j].imag();
                        }
                        k = k + 1;
                    }
                }
                //std::cout << std::endl << Q-overlapping_signals << std::endl;
                //std::cout << TNAA_ind << std::endl;
                //std::cout << Q-overlapping_signals+extra_cols << std::endl;

                // do the combined signals now
                if ( TNAA )
                {
                    double a1 = ahat(NAA_pos)/(ahat(NAA_pos)+ahat(NAAG_pos));
                    double a2 = ahat(NAAG_pos)/(ahat(NAA_pos)+ahat(NAAG_pos));

                    if ((ahat(NAA_pos) + ahat(NAAG_pos)) == 0)
                    {
                        a1 = 0.5;
                        a2 = 0.5;
                    }

                    for( int i = 1; i <= activeN; ++i )
                    {
                        if ((ahat(NAA_pos) + ahat(NAAG_pos)) == 0)
                        {
                            a1 = 0.5;
                            a2 = 0.5;
                        }

                        D_comb[i][TNAA_ind+Q-overlapping_signals] = a1 * SpOut[i+nStart-1][NAA_pos].real() + a2 * SpOut[i+nStart-1][NAAG_pos].real();
                        D_comb[i+JR.msize()/2][TNAA_ind+Q-overlapping_signals] = a1 * SpOut[i+nStart-1][NAA_pos].imag() + a2 * SpOut[i+nStart-1][NAAG_pos].imag();
                        ahat_comb(TNAA_ind) = ahat(NAA_pos) + ahat(NAAG_pos);
                    }
                }
                if ( TCho )
                {
                    double a1 = ahat(PCh_pos)/(ahat(PCh_pos)+ahat(GPC_pos));
                    double a2 = ahat(GPC_pos)/(ahat(PCh_pos)+ahat(GPC_pos));

                    if ((ahat(PCh_pos) + ahat(GPC_pos)) == 0)
                    {
                        a1 = 0.5;
                        a2 = 0.5;
                    }

                    for( int i = 1; i <= activeN; ++i )
                    {
                        D_comb[i][TCho_ind+Q-overlapping_signals] = a1 * SpOut[i+nStart-1][PCh_pos].real() + a2 * SpOut[i+nStart-1][GPC_pos].real();
                        D_comb[i+JR.msize()/2][TCho_ind+Q-overlapping_signals] = a1 * SpOut[i+nStart-1][PCh_pos].imag() + a2 * SpOut[i+nStart-1][GPC_pos].imag();
                        ahat_comb(TCho_ind) = ahat(PCh_pos) + ahat(GPC_pos);
                    }
                }
                if ( TCr )
                {
                    double a1 = ahat(Cr_pos)/(ahat(Cr_pos)+ahat(PCr_pos));
                    double a2 = ahat(PCr_pos)/(ahat(Cr_pos)+ahat(PCr_pos));

                    if ((ahat(Cr_pos) + ahat(PCr_pos)) == 0)
                    {
                        a1 = 0.5;
                        a2 = 0.5;
                    }

                    for( int i = 1; i <= activeN; ++i )
                    {
                        D_comb[i][TCr_ind+Q-overlapping_signals] = a1 * SpOut[i+nStart-1][Cr_pos].real() + a2 * SpOut[i+nStart-1][PCr_pos].real();
                        D_comb[i+JR.msize()/2][TCr_ind+Q-overlapping_signals] = a1 * SpOut[i+nStart-1][Cr_pos].imag() + a2 * SpOut[i+nStart-1][PCr_pos].imag();
                        ahat_comb(TCr_ind) = ahat(Cr_pos) + ahat(PCr_pos);
                    }
                }
                if ( Glx )
                {
                    double a1 = ahat(Glu_pos)/(ahat(Glu_pos)+ahat(Gln_pos));
                    double a2 = ahat(Gln_pos)/(ahat(Glu_pos)+ahat(Gln_pos));

                    if ((ahat(Glu_pos) + ahat(Gln_pos)) == 0)
                    {
                        a1 = 0.5;
                        a2 = 0.5;
                    }
                    for( int i = 1; i <= activeN; ++i )
                    {
                        D_comb[i][Glx_ind+Q-overlapping_signals] = a1 * SpOut[i+nStart-1][Glu_pos].real() + a2 * SpOut[i+nStart-1][Gln_pos].real();
                        D_comb[i+JR.msize()/2][Glx_ind+Q-overlapping_signals] = a1 * SpOut[i+nStart-1][Glu_pos].imag() + a2 * SpOut[i+nStart-1][Gln_pos].imag();
                        ahat_comb(Glx_ind) = ahat(Glu_pos) + ahat(Gln_pos);
                    }
                }
                if ( TLM09 )
                {
                    double a1 = ahat(Lip09_pos)/(ahat(Lip09_pos)+ahat(MM09_pos));
                    double a2 = ahat(MM09_pos)/(ahat(Lip09_pos)+ahat(MM09_pos));

                    if ((ahat(Lip09_pos) + ahat(MM09_pos)) == 0)
                    {
                        a1 = 0.5;
                        a2 = 0.5;
                    }
                    for( int i = 1; i <= activeN; ++i )
                    {
                        D_comb[i][TLM09_ind+Q-overlapping_signals] = a1 * SpOut[i+nStart-1][Lip09_pos].real() + a2 * SpOut[i+nStart-1][MM09_pos].real();
                        D_comb[i+JR.msize()/2][TLM09_ind+Q-overlapping_signals] = a1 * SpOut[i+nStart-1][Lip09_pos].imag() + a2 * SpOut[i+nStart-1][MM09_pos].imag();
                        ahat_comb(TLM09_ind) = ahat(Lip09_pos) + ahat(MM09_pos);
                    }
                }
                /*if ( TL13 )
                {
                    for( int i = 1; i <= activeN; ++i )
                    {
                        D_comb[i][TL13_ind+Q-overlapping_signals] = ahat(Lip13a_pos)/(ahat(Lip13a_pos)+ahat(Lip13b_pos)) * SpOut[i+nStart-1][Lip13a_pos] + ahat(Lip13b_pos)/(ahat(Lip13a_pos)+ahat(Lip13b_pos)) * SpOut[i+nStart-1][Lip13b_pos];
                        ahat_comb(TL13_ind) = ahat(Lip13a_pos) + ahat(Lip13b_pos);
                    }
                }*/
                if ( TLM13 )
                {
                    double a1 = ahat(Lip13a_pos)/(ahat(Lip13a_pos)+ahat(Lip13b_pos)+ahat(MM12_pos)+ahat(MM14_pos));
                    double a2 = ahat(Lip13b_pos)/(ahat(Lip13a_pos)+ahat(Lip13b_pos)+ahat(MM12_pos)+ahat(MM14_pos));
                    double a3 = ahat(MM12_pos)  /(ahat(Lip13a_pos)+ahat(Lip13b_pos)+ahat(MM12_pos)+ahat(MM14_pos));
                    double a4 = ahat(MM14_pos)  /(ahat(Lip13a_pos)+ahat(Lip13b_pos)+ahat(MM12_pos)+ahat(MM14_pos));

                    if ((ahat(Lip13a_pos) + ahat(Lip13b_pos) + ahat(MM12_pos) + ahat(MM14_pos)) == 0)
                    {
                        a1 = 0.25;
                        a2 = 0.25;
                        a3 = 0.25;
                        a4 = 0.25;
                    }

                    for( int i = 1; i <= activeN; ++i )
                    {
                        // HAHAHA Greg will hate this!
                        D_comb[i][TLM13_ind+Q-overlapping_signals] = a1 * SpOut[i+nStart-1][Lip13a_pos].real() + a2 * SpOut[i+nStart-1][Lip13b_pos].real() + a3 * SpOut[i+nStart-1][MM12_pos].real() + a4 * SpOut[i+nStart-1][MM14_pos].real();
                        D_comb[i+JR.msize()/2][TLM13_ind+Q-overlapping_signals] = a1 * SpOut[i+nStart-1][Lip13a_pos].imag() + a2 * SpOut[i+nStart-1][Lip13b_pos].imag() + a3 * SpOut[i+nStart-1][MM12_pos].imag() + a4 * SpOut[i+nStart-1][MM14_pos].imag();

                        ahat_comb(TLM13_ind) = ahat(Lip13a_pos) + ahat(Lip13b_pos) + ahat(MM12_pos) + ahat(MM14_pos);
                    }
                }

                if ( TLM20 )
                {
                    double a1 = ahat(Lip20_pos)/(ahat(Lip20_pos)+ahat(MM20_pos));
                    double a2 = ahat(MM09_pos)/(ahat(Lip20_pos)+ahat(MM20_pos));

                    if ((ahat(Lip20_pos) + ahat(MM20_pos)) == 0)
                    {
                        a1 = 0.5;
                        a2 = 0.5;
                    }
                    for( int i = 1; i <= activeN; ++i )
                    {
                        D_comb[i][TLM20_ind+Q-overlapping_signals] = a1 * SpOut[i+nStart-1][Lip20_pos].real() + a2 * SpOut[i+nStart-1][MM20_pos].real();
                        D_comb[i+JR.msize()/2][TLM20_ind+Q-overlapping_signals] = a1 * SpOut[i+nStart-1][Lip20_pos].imag() + a2 * SpOut[i+nStart-1][MM20_pos].imag();
                        ahat_comb(TLM20_ind) = ahat(Lip20_pos) + ahat(MM20_pos);
                    }
                }

                if ( TGABA )
                {
                    double a1 = ahat(GABAA_pos)/(ahat(GABAA_pos)+ahat(GABAB_pos));
                    double a2 = ahat(GABAB_pos)/(ahat(GABAA_pos)+ahat(GABAB_pos));

                    if ((ahat(GABAA_pos) + ahat(GABAB_pos)) == 0)
                    {
                        a1 = 0.5;
                        a2 = 0.5;
                    }
                    for( int i = 1; i <= activeN; ++i )
                    {
                        D_comb[i][TGABA_ind+Q-overlapping_signals] = a1 * SpOut[i+nStart-1][GABAA_pos].real() + a2 * SpOut[i+nStart-1][GABAB_pos].real();
                        D_comb[i+JR.msize()/2][TGABA_ind+Q-overlapping_signals] = a1 * SpOut[i+nStart-1][GABAA_pos].imag() + a2 * SpOut[i+nStart-1][GABAB_pos].imag();
                        ahat_comb(TGABA_ind) = ahat(GABAA_pos) + ahat(GABAB_pos);
                    }
                }


                for( int i = 1; i <= activeN*2; ++i )
                {
                    // other paras
                    for( int k = 1; k<= JR.nsize(); ++k )
                    {
                        D_comb[i][k+Q-overlapping_signals+extra_cols] = JR[i][k];
                        //std::cout << k+Q-overlapping_signals+extra_cols << std::endl;
                    }
                }
                
                /*for( int k = 1; k<= JR.nsize(); ++k )
                {
                    std::cout << k << std::endl;
                    std::cout << D_comb[1][k+Q-overlapping_signals+extra_cols] << std::endl;
                }*/

                //std::cout << JR << std::endl;


                /*std::ofstream dmat_comb;
                dmat_comb.open ("./Dmat_comb.txt");
                for( int i = 1; i <= D_comb.msize(); ++i )
                    for( int j = 1; j<= D_comb.nsize(); ++j )
                        dmat_comb << D_comb(i,j).real() << std::endl << " " << D_comb(i,j).imag() << std::endl;

                dmat_comb.close();
                */
            
                // import numpy as np
                // import pylab as pl
                // D_comb = np.matrix(np.reshape(np.loadtxt('Dmat_comb.txt'),(986,73)))
                // FC = D_comb.T*D_comb
                // fisher = FC/ 5.214470e-01
                // IFISHER = np.linalg.pinv(fisher)
                // pl.plot(D[:,60])
                // pl.show()

            
                /*std::ofstream dmat_comb;
                dmat_comb.open ("./Dmat_comb.txt");
                for( int i = 1; i <= D_comb.msize(); ++i )
                    for( int j = 1; j<= D_comb.nsize(); ++j )
                        dmat_comb << D_comb(i,j) << std::endl;

                dmat_comb.close();

                std::cout << D_comb.msize() << "," << D_comb.nsize() << std::endl;*/
                  
                // measure CRLBs
                cvm::srmatrix FC_comb = (~D_comb)*D_comb;

                cvm::srmatrix fisher_comb = FC_comb;
                fisher_comb *= 1.0 / sigma;

                cvm::rmatrix IFISHER_comb = fisher_comb.pinv();

                cvm::rvector crlbs_comb;
                crlbs_comb.resize(IFISHER_comb.msize());
                //std::cout << IFISHER_comb.nsize() << " " << IFISHER_comb.msize() << std::flush << std::endl;

                for( int l = 1; l < IFISHER_comb.msize()+1; ++l )
                    crlbs_comb[l] = std::sqrt(IFISHER_comb[l][l]);

                //std::cout << crlbs_comb << std::endl;
                
                // pick out the combined subset
                cvm::rvector crlbs_comb_cut;
                crlbs_comb_cut.resize(extra_cols);

                for( int x = 1; x < extra_cols + 1; x++ )
                    crlbs_comb_cut(x) = crlbs_comb(Q-overlapping_signals+x);
                
                rvec_stdvec& crlbs_vec_comb = work.GetCRLBsComb();
                crlbs_vec_comb.push_back(crlbs_comb_cut);
                
                rvec_stdvec& ahat_vec_comb = work.GetAmplitudesComb();
                ahat_vec_comb.push_back(ahat_comb);
            
                rvec_stdvec& ampl_norm_vec_comb = work.GetAmplitudesNormalisedComb();
                cvm::rvector ampl_norm_comb;
                rvec_stdvec& crlb_norm_vec_comb = work.GetCRLBsNormalisedComb();
                cvm::rvector crlb_norm_comb;

                ampl_norm_comb.resize( ahat_comb.size() );
                crlb_norm_comb.resize( crlbs_comb_cut.size() );
                

                // scale the amplitudes
                for( int i = 1; i <= ahat_comb.size(); ++i )
                    ampl_norm_comb[i] = work.NormaliseValue(voxel_num-1, ahat_comb[i]);

                for( int i = 1; i <= crlbs_comb_cut.size(); ++i )
                    crlb_norm_comb[i] = work.NormaliseValue(voxel_num-1, crlbs_comb_cut[i]);

                // store the results before moving to the next fid
                ampl_norm_vec_comb.push_back(ampl_norm_comb);
                crlb_norm_vec_comb.push_back(crlb_norm_comb);
            }



			//
			// Compute the normalised value
			//

			rvec_stdvec& ampl_norm_vec = work.GetAmplitudesNormalised();
			cvm::rvector ampl_norm;
			rvec_stdvec& crlb_norm_vec = work.GetCRLBsNormalised();
			cvm::rvector crlb_norm;

			ampl_norm.resize( ahat.size() );
			crlb_norm.resize( ahat.size() );

			// scale the amplitudes
			for( int i = 1; i <= ahat.size(); ++i )
				ampl_norm[i] = work.NormaliseValue(voxel_num-1, ahat[i]);

			for( int i = 1; i <= crlbs.size(); ++i )
				crlb_norm[i] = work.NormaliseValue(voxel_num-1, crlbs[i]);

			// store the results before moving to the next fid
			ampl_norm_vec.push_back(ampl_norm);
			crlb_norm_vec.push_back(crlb_norm);
			yhat_vec.push_back(yhat);
			GpOut_vec.push_back(GpOut);
			SpOut_vec.push_back(SpOut);
			ahat_vec.push_back(ahat);
		}
	
		return true;
	}
	catch( const std::exception& e )
	{
		log.LogMessage(LOG_ERROR, "failed during RunTARQUIN: %s", e.what());
		return false;
	}

	return false;
}


	// perform self-deconvolution kept for future development
	/*
    cvm::cvector linesh(yhat.size());
    cvm::cvector linesh_smo(yhat.size());
	cvm::cvector delta(yhat.size());
	cvm::cvector DELTA(yhat.size());

    for (int iters = 1; iters < 1; iters++)
    {
        linesh.set(tcomplex(1,0));
        for (int n = 1; n < yhat.size() + 1; n++)
        {
        //    if ( abs(y(n)/yhat(n)) > 10 )
        //        linesh(n) = 1;
        //    else
        //        linesh(n) = (y(n))/((yhat(n)));
        
            float s = 0.001f;
            tcomplex z(s * std::abs(y(1)), s * std::abs(y(1)));

            linesh(n) = y(n) / (z + yhat(n));
        }
        //plot(linesh);
        //plot(linesh);
        //plot(y);
        td_conv_ws_noext(linesh, linesh_smo, fid.GetSamplingFrequency()/20,10);
        //td_conv_ws_noext(linesh, linesh_smo, 20, 10);
        //for ( int n = 1; n < fid.GetSamplingFrequency()/10 + 1; n++ )
        //    linesh_smo(n) = linesh(n);

        cvm::cmatrix mat(yhat.size(),2);
        mat(1) = linesh;
        mat(2) = linesh_smo;
        //plot(mat);
        //plot(linesh_smo);
        // simulate a delta function
        for (int n = 1; n < yhat.size() + 1; n++)
            delta(n) = tcomplex(1,0);

        treal beta  = vParams(nIdxBeta);
        for (int n = 1; n < yhat.size() + 1; n++)
        {
            yhat(n) = yhat(n) * linesh_smo(n);
            //y(n) = y(n) / linesh_smo(n);
            //delta(n) = delta(n) * exp ( -t(n) * beta * beta );
            delta(n) = delta(n) * linesh_smo(n) * exp (  -t(n) * t(n) * beta );
        }
        fft(delta, DELTA);
        DELTA = fftshift(DELTA);
        //plot(delta);
        //plot(DELTA);
        // apply to basis matrix (plotting)
        for (int n = 1; n < SpOut.nsize() + 1; n++)
            for (int m = 1; m < SpOut.msize() + 1; m++)
                SpOut(m,n) = SpOut(m,n) * linesh_smo(m);

    }
	*/

