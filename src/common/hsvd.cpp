#include "hsvd.hpp"
#include "cvm_util.hpp"
#include <complex> 

// for theory see Barkhuijsen paper or chapter three of Greg's thesis
void tarquin::DecomposeHSVD(
        const cvm::cvector& y,
        treal fs,
		cvm::rvector& freqs, 
		cvm::rvector& dampings, 
		cvm::cmatrix& basis, 
		cvm::cvector& ahat, 
		integer       nPts, 
		integer       K, 
		CBoswell&     log
		)
{
	//const cvm::cvector& y = fid.GetVectorFID(proc_coord);

	// scale y by largest real point
	/*integer maxpt = y.indofmax();
	  for( integer n = 0; n < y.size(); n++ )
	  y(n+1) = y(n+1) / 0.0896; 
	  */

	integer N = nPts;

	// this is the rank we are truncating to
	// (call rank estimation method here)
	// good for 3T Philips
	//integer K = 50;
	//integer K = 100;
	// good for Siemens
	//integer K = 10;

	// output objects are sized appropriately
	freqs.resize(K);
	dampings.resize(K);
	basis.resize(y.size(), K);
	cvm::cmatrix basisSmall(N, K);
	ahat.resize(K);

	// the Hankel matrix must be as square as possible 
	integer L = N / 2;
	integer M = N+1-L;

	// H is the L x M Hankel matrix
	cvm::cmatrix H(L, M);

	// arrange the samples from the signal into the Hankel matrix
	{
		Task progress(log, "building Hankel LP matrix");
		MakeHankelMatrix(y, H);
	}

	// do the SVD of H
	cvm::scmatrix U(L);
	cvm::scmatrix VH(M);
	cvm::rvector s(L);

	{
		Task progress(log, "SVD of LP matrix");
		s = H.svd(U, VH);
	}

	// Ukt is the rank K submatrix of U, without the top rowa (Ukt is L-1 x K)
	cvm::cmatrix Ukt(U, 2, 1, L-1, K);

	// Ukb is the rank K submatrix of U, without the bottom row (Ukt is L-1 x K)
	cvm::cmatrix Ukb(U, 1, 1, L-1, K);

	// eigenvalues of this matrix are the poles of the transfer function
	cvm::scmatrix Zp = Ukb.pinv() * Ukt;

	cvm::cvector q(K);
	{
		Task progress(log, "estimating LP parameters");
		q.eig(Zp);
	}

	// time step
	treal dt = 1.0 / fs;

	// make non-linearly entering parameters
	for( integer i = 0; i < K; i++ ) 
	{
		tcomplex qi = q[i+1];

		// ln(r*e^{j\theta}) = ln(r) + j*theta
		treal r     = std::log(abs(qi));
		treal theta = atan2(qi.imag(), qi.real());

		dampings[i+1] = r / dt;
        // remove any +ve dampings
        if ( dampings[i+1] > 0 )
            dampings[i+1] = 0;

		freqs[i+1]    = theta / (2.0*M_PI*dt);

		// make basis (small basis for amplitude estimation)
		for( integer n = 0; n < N; n++ ) 
		{
			tcomplex phi(dt*n*dampings[i+1], dt*n*2.0*M_PI*freqs[i+1]);
			basisSmall(n+1, i+1) = exp(phi); 
		}
	}

	// amplitude estimates
	cvm::cvector ysub(N);
	{
		Task progress(log, "estimating LP amplitudes");
		for( integer n = 0; n < N;n ++)
			ysub(n+1) = y(n+1);

		ahat = basisSmall.pinv()*ysub;
	}


	// make output basis
	for( integer i = 0; i < K; i++ ) 
	{
		// make basis (small basis for amplitude estimation)
		for( integer n = 0; n < y.size(); n++ ) 
		{
			tcomplex phi(dt*n*dampings[i+1], dt*n*2.0*M_PI*freqs[i+1]);
			basis(n+1, i+1) = exp(phi); 
		}
	}

}

// assumes H is already the right size
void tarquin::MakeHankelMatrix(const cvm::cvector& y, cvm::cmatrix& H)
{
	integer L = H.msize();
	integer M = H.nsize();

	for( int r = 1; r <= L; r++ ) 
	{
		for( int c = 1; c <= M; c++ ) 
		{
			assert( c + r -1 <= y.size() );
			H(r, c) = y(c + r -1);
		}
	}
}

