#ifndef __CNNLSSOLVER__
#define __CNNLSSOLVER__

#include "common.hpp"

namespace tarquin {

    /*!
     * Abstract base class functor for solving the NNLS problem.
     */
    class CNNLSSolver {

	public:

	    //! Constructor.
	    CNNLSSolver()
	    {

	    }

	    //! Destructor must be virtual :)
	    virtual ~CNNLSSolver()
	    {

	    }

	    /*!
	     * Solve the NNLS problem with all new data.
	     * \param mA is the matrix on the left hand side.
	     * \param vb is the vector on the right hand side.
	     * \param x is the result vector.
	     * \param bUseStartingValue uses columns matching non-zero elements of x as initial solution.
	     */
	    virtual void solve(cvm::rmatrix& mA, cvm::rvector& vb, cvm::rvector& vx, bool bUseStartingValue) = 0;
    };

}

#endif
