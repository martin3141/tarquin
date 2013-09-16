#ifndef __HSVD__
#define __HSVD__

#include "CFID.hpp"
#include "CBoswell.hpp"

namespace tarquin {

    /*!
     * Functions for modelling the signal using the HSVD method.
     */

    /*!
     * Decompose the FID passed in into a basis of complex exponentials.
     * \param freqs is a vector of frequencies.
     * \param dampings is a vector of dampings.
     * \param basis is a matrix of basis vectors.
     * \param ahat is an estimate of the amplitudes of each basis vector.
     * \param nPts is the number of points to use.
     * \param log is the logging object where messages go, etc.
     */
    void DecomposeHSVD(const cvm::cvector& y, treal fs,
	    cvm::rvector& freqs, cvm::rvector& dampings, cvm::cmatrix& basis, cvm::cvector& ahat, integer nPts, integer K, CBoswell& log);

    // this could be templated to work for real and complex
    void MakeHankelMatrix(const cvm::cvector& y, cvm::cmatrix& H);

}

#endif
