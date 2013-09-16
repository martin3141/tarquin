#ifndef __CNNLSSOLVERLH__
#define __CNNLSSOLVERLH__

#include "CNNLSSolver.hpp"

namespace tarquin {

    /*!
     * Concerete class implementing the Lawson-Hanson NNLS solver.
     */
    class CNNLSSolverLH : public CNNLSSolver {

	public:

	    /*!
	     * Solve the NNLS problem with all new data.
	     * \param mA is the matrix on the right hand side.
	     * \param vb is the vector on the left hand side.
	     * \param x is the result vector.
	     */
	    void solve(cvm::rmatrix& mA, cvm::rvector& vb, cvm::rvector& vx, bool bUseStartingValue = false); 


	private:

	    //! Pseudo-inverse of some submatrix of the right hand side.
	    cvm::rmatrix m_mX;

	    //! Candidate solution to unconstrained sub problem.
	    cvm::rvector m_vz; 

	    //! Gradient of error, e = Ax - b 
	    cvm::rvector m_vw; 
    };

}

#endif
