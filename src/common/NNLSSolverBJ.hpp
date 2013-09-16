#ifndef CNNLSSOLVERBJ
#define CNNLSSOLVERBJ

#include "CNNLSSolver.hpp"

namespace tarquin 
{
    /*!
     * Concerete class implementing a not-quite Bro and de Jong solver.
     */
    class NNLSSolverBJ : public CNNLSSolver 
    {
	public:

	    /*!
	     * Solve the NNLS problem with all new data.
	     * \param mA is the matrix on the right hand side.
	     * \param vb is the vector on the left hand side.
	     * \param x is the result vector.
	     */
	    void solve(cvm::rmatrix& mA, cvm::rvector& vb, cvm::rvector& vx, bool bUseStartingValue = false); 

    };
}

#endif // CNNLSSOLVERBJ
