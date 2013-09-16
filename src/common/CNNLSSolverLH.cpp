#include "CNNLSSolverLH.hpp"
#include <set>
#include <cmath>

namespace tarquin {

    // functions to obtain active set of columns
    // give starting active set

    void CNNLSSolverLH::solve(cvm::rmatrix& mA, cvm::rvector& vb, cvm::rvector& vx, bool bUseStartingValue)
    {
	assert( vx.size() == mA.nsize() );

	//m_mX.resize(mA.nsize(), mA.msize());
	//m_mX = mA.pinv(0.01);
	//vx = m_mX*vb;
	//return;

	// the set of elements in x that are 0 (set of columns of A not used)
	std::set<integer> Z;

	// the set of elements in x that are > 0 (set of columns of A used)
	std::set<integer> P;

	// the submatrix of A with some columns maybe set to zero 
	cvm::rmatrix mAp(mA.msize(), mA.nsize());

	// the candidate solution to the unconstrained problem
	m_vz.resize(vx.size());

	// the gradient of the error
	m_vw.resize(vx.size());

	// the pseudo-inverse of A
	m_mX.resize(mA.nsize(), mA.msize());

	bool bSkipToStep6 = false;

	// ACHTUNG: You may have to make this bigger if convergence of the NNLS algorithm is a problem.
	treal tol = NUMERICAL_TOL_NNLS;

	// initialise Z and x (step 1) if not specified 
	if( false == bUseStartingValue ) {

	    for(integer i = 1; i <= mA.nsize(); i++) {
		Z.insert(i);

		vx[i] = 0.0;
	    }
	}
	// setup index sets based on elements of x
	else {

	    for(integer i = 1; i <= mA.nsize(); i++) {

		if( vx[i] > tol )
		    P.insert(i);
		else {

		    Z.insert(i);
		}
	    }

	    bSkipToStep6 = true;
	}
    
    bool bFirstPass = true; // MW 28th June 13, in rare cases (ie all amps are +ve) this needs to be
                            // added to make sure the new values are re-optimised rather than just getting skipped

    treal wmax_old = 0.0;
	// Continue until "zero" set is empty, i.e. no variables are fixed at zero.
	// There are other ways to stop of course, in the case where some elements have to be zero.
	while( Z.size() > 0 || ( bUseStartingValue && bFirstPass) ) {

        //std::cout << Z.size() << std::endl;
        if ( bFirstPass )
            bFirstPass = !bFirstPass;

	    if( false == bSkipToStep6 ) {

		// compute negative of gradient (step 2) (we could optimise this by precomputing some elements
		m_vw = (~mA)*(vb - mA*vx);

		//std::cout << "\ngradient = " << -m_vw << std::flush;

		// stop if all gradients (with indices in Z) are negative or zero (step 3)
		// also find biggest gradient (step 4)
		integer t = -1;;
		treal wmax = 0.0;

		for( std::set<integer>::iterator itj = Z.begin(); itj != Z.end(); itj++ ) {

		    integer j = *itj;

		    if( m_vw(j) >= wmax ) {
			t = j;
			wmax = m_vw(j);
		    }
		}

        //std::cout << std::cout.precision(15) << "wmax = " << wmax << std::endl;
        //std::cout << "tol = " << tol << std::endl;

		// if biggest gradient is positive in all directions (or really small) then we can stop 
		if( wmax < tol )
		    return;

        // if the gradient is identical to the previous value then we are probally stuck in a loop
        if( wmax == wmax_old )
        {
            std::cout << std::endl << "Warning, NNLS did not fully converge." << std::endl;
		    return;
        }

        wmax_old = wmax;

		// sanity check that t is in Z
		assert( t != -1 );

		// move index t from Z to set P (step 5)
		P.insert(t);
		Z.erase(t);
	    }

	    // construct Ep based on columns of E that have indices in P (step 6)
	    mAp.set(0.0);
	    for( std::set<integer>::iterator itj = P.begin(); itj != P.end(); itj++ )	
		mAp(*itj) = mA(*itj);

	    // find solution z that is solution to Ep*z = f (step 6)
	    // NOTE: if this tolerance is small (e.g. the default) then although the results look
	    // plausible, the gradient never appears to be very accurate.
	    m_mX = mAp.pinv(tol);
	    m_vz = m_mX*vb;

	    assert( m_vz.size() == mA.nsize() );

	    // set all elements of z to zero for indices in Z
	    //std::cout << "\nZ contains:" << std::flush;
	    for( std::set<integer>::iterator itj = Z.begin(); itj != Z.end(); itj++ ) {
		//std::cout << " " << *itj;
		m_vz(*itj) = 0.0;
	    }

	    // (step 7), if z_j > 0 for all j \in P then we jump back to the top and carry on
	    bool bJump = true;

	    //std::cout << "\nP contains:" << std::flush;
	    for( std::set<integer>::iterator itj = P.begin(); itj != P.end(); itj++ )	{

		//std::cout << "P(" << *itj << ") = " << m_vz(*itj) << std::endl;
		if( m_vz(*itj) <= tol ) {
		    bJump = false;
		    break;
		}
	    }

	    // (step 7) initialise x to z then carry on
	    if( true ==  bJump ) {
		vx = m_vz;
		bSkipToStep6 = false;
	    }
	    else {

		// NOTE: the inner loop doesn't happen very often. As such it is a bit difficult to test
		// however, in the cases where it has happened, it appears to have worked...
		//std::cout << "\n INNER LOOP FIRED! " << std::flush;

		// (step 8): find index q \in P s.t. x_q/(x_q-z_q) is min
		treal alpha = std::numeric_limits<treal>::infinity();

		integer q = -1;
		for( std::set<integer>::iterator itj = P.begin(); itj != P.end(); itj++ ) {

		    integer j = *itj;

		    if( m_vz(j) <= tol ) {			
			treal quotient = vx(j)/(vx(j)-m_vz(j));

			if( quotient < alpha ) {
			    alpha = quotient;
			    q = j;
			}
		    }
		}

		// sanity check
		assert( q != -1 );

		// (step 9): already done, stored alpha from when we were looking for it

		// (step 10): update x (I think there is an efficient BLAS type way of doing this)
		vx += alpha*(m_vz-vx);

		// (step 11): move all indices from P back to Z which have zero solutions (x_j == 0.0)
		std::set<integer> Pdel;
		for( std::set<integer>::iterator itj = P.begin(); itj != P.end(); itj++ )	{

		    integer j = *itj;

		    if( std::abs(vx(j)) < tol ) {
			Z.insert(j);
			Pdel.insert(j);
		    }
		}

		// use set difference for this loop
		// remove elements from P
		for( std::set<integer>::iterator itj = Pdel.begin(); itj != Pdel.end(); itj++ )	
		    P.erase(*itj);

		bSkipToStep6 = true;	
	    }
	}
    }

}
