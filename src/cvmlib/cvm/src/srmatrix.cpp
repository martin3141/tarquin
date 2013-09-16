/* CVM Class Library */
/* Copyright (C), Sergei Nikolaev, 1992-2006, http://cvmlib.com */

#include "cvm.h"
#include "blas.h"

#if defined (_MSC_VER)
#   pragma warning(disable:4786)
#endif

CVM_NAMESPACE_BEG

template<>
CVM_API void
__exp<basic_srmatrix<float>, float>
    (basic_srmatrix<float>& m,
    const basic_srmatrix<float>& mArg,
    float tol) throw (cvmexception)
{
    int nR = 0, nI = 0, nQ = 0, nJ = 0;
    const int nM = m.msize();
    if (nM != mArg.msize()) throw cvmexception (CVM_SIZESMISMATCH);

    basic_srmatrix<float> mTmp;
    const float* pD = mArg._pd();

    if (pD == m.get())
    {
        mTmp << mArg;
        pD = mTmp;
    }

    SMEXPC (&nM, pD, pD == m.get() ? mTmp._pld() : mArg._pldm(), &tol, &nR, &nI, &nQ, &nJ);
    basic_rvector<float> vR (nR);
    basic_array<int> vI (nI);

    const int issymm = 0;
    float work_dummy = 0.F;
    const int lwork_dummy = 0;

    SMEXP (&nM, pD, pD == m.get() ? mTmp._pld() : mArg._pldm(), m, m._pld(), 
           vR, vI, &nR, &nI, &nQ, &nJ, &issymm, &work_dummy, &lwork_dummy);
}

template<>
CVM_API void 
__exp<basic_srmatrix<double>, double>
    (basic_srmatrix<double>& m,
    const basic_srmatrix<double>& mArg,
    double tol) throw (cvmexception)
{
    int nR = 0, nI = 0, nQ = 0, nJ = 0;
    const int nM = m.msize();
    if (nM != mArg.msize()) throw cvmexception (CVM_SIZESMISMATCH);

    basic_srmatrix<double> mTmp;
    const double* pD = mArg._pd();

    if (pD == m.get())
    {
        mTmp << mArg;
        pD = mTmp;
    }

    DMEXPC (&nM, pD, pD == m.get() ? mTmp._pld() : mArg._pldm(), &tol, &nR, &nI, &nQ, &nJ);
    basic_rvector<double> vR (nR);
    basic_array<int> vI (nI);

    const int issymm = 0;
    double work_dummy = 0.;
    const int lwork_dummy = 0;

    DMEXP (&nM, pD, pD == m.get() ? mTmp._pld() : mArg._pldm(), m, m._pld(),
           vR, vI, &nR, &nI, &nQ, &nJ, &issymm, &work_dummy, &lwork_dummy);
}

template<>
CVM_API void
__exp_symm<basic_srsmatrix<float>, float>
    (basic_srsmatrix<float>& m, 
    const basic_srsmatrix<float>& mArg, 
    float tol) throw (cvmexception)
{
    int nR = 0, nI = 0, nQ = 0, nJ = 0;
    const int nM = m.msize();
    if (nM != mArg.msize()) throw cvmexception (CVM_SIZESMISMATCH);

    basic_srmatrix<float> mTmp;
    const float* pD = mArg._pd();

    if (pD == m.get())
    {
        mTmp << mArg;
        pD = mTmp;
    }

    SMEXPC (&nM, pD, pD == m.get() ? mTmp._pld() : mArg._pldm(), &tol, &nR, &nI, &nQ, &nJ);
    basic_rvector<float> vR (nR);
    basic_array<int> vI (nI);

    const int issymm = 1;
    const int lwork  = 64 * nM;
    basic_rvector<float> work (lwork);

    SMEXP (&nM, pD, pD == m.get() ? mTmp._pld() : mArg._pldm(), m, m._pld(),
           vR, vI, &nR, &nI, &nQ, &nJ, &issymm, work, &lwork);
}

template<>
CVM_API void
__exp_symm<basic_srsmatrix<double>, double>
    (basic_srsmatrix<double>& m,
    const basic_srsmatrix<double>& mArg,
    double tol) throw (cvmexception)
{
    int nR = 0, nI = 0, nQ = 0, nJ = 0;
    const int nM = m.msize();
    if (nM != mArg.msize()) throw cvmexception (CVM_SIZESMISMATCH);

    basic_srmatrix<double> mTmp;
    const double* pD = mArg._pd();

    if (pD == m.get())
    {
        mTmp << mArg;
        pD = mTmp;
    }

    DMEXPC (&nM, pD, pD == m.get() ? mTmp._pld() : mArg._pldm(), &tol, &nR, &nI, &nQ, &nJ);
    basic_rvector<double> vR (nR);
    basic_array<int> vI (nI);

    const int issymm = 1;
    const int lwork  = 64 * nM;
    basic_rvector<double> work (lwork);

    DMEXP (&nM, pD, pD == m.get() ? mTmp._pld() : mArg._pldm(), m, m._pld(),
           vR, vI, &nR, &nI, &nQ, &nJ, &issymm, work, &lwork);
}

template<>
CVM_API void
__cond_num<float, basic_srmatrix<float> >
    (const basic_srmatrix<float>& mArg, float& dCond) throw (cvmexception)
{
    dCond = 0.F;
    const int mnM  = mArg.msize();
    int   nOutInfo = 0;
    basic_srmatrix<float> mA (mArg);
    basic_rvector<float> work (mnM * 4);
    basic_array<int> iwork (mnM);

    const float rNorm = mA.norminf();
    SGETRF (&mnM, &mnM, mA, mA._pld(), iwork, &nOutInfo);

    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
    if (nOutInfo == 0)
    {
        SGECON (Chars::pI(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                1,
#endif
                &mnM, mA, mA._pld(), &rNorm, &dCond, work, iwork, &nOutInfo);
    }
}

template<>
CVM_API void
__cond_num<double, basic_srmatrix<double> >
    (const basic_srmatrix<double>& mArg, double& dCond) throw (cvmexception)
{
    dCond = 0.;
    const int mnM   = mArg.msize();
    int    nOutInfo = 0;
    basic_srmatrix<double> mA (mArg);
    basic_rvector<double> work (mnM * 4);
    basic_array<int> iwork (mnM);

    const double rNorm = mA.norminf();
    DGETRF (&mnM, &mnM, mA, mA._pld(), iwork, &nOutInfo);

    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
    if (nOutInfo == 0)
    {
        DGECON (Chars::pI(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                1,
#endif
                &mnM, mA, mA._pld(), &rNorm, &dCond, work, iwork, &nOutInfo);
    }
}

template<>
CVM_API void
__inv<basic_srmatrix<float> >
    (basic_srmatrix<float>& m,
    const basic_srmatrix<float>& mArg) throw (cvmexception)
{
    const int mnM = m.msize();
    if (mnM != mArg.msize()) throw cvmexception (CVM_SIZESMISMATCH);
    if (mnM == 1)
    {
        static const float one(1.F);
        if (_abs (mArg(1,1)) <= basic_cvmMachMin<float>()) throw cvmexception (CVM_SINGULARMATRIX, 1);
        m(1,1) = one / mArg(1,1);
    }
    else
    {
        const int nWorkSize = mnM * 64;
        int nOutInfo = 0;
        basic_array<int> nPivots (mnM);
        basic_rvector<float> vWork (nWorkSize);

        m.low_up (mArg, nPivots);
        SGETRI (&mnM, m, m._pld(), nPivots, vWork, &nWorkSize, &nOutInfo);
        if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
        if (nOutInfo > 0) throw cvmexception (CVM_SINGULARMATRIX, nOutInfo);
    }
}

template<>
CVM_API void
__inv<basic_srmatrix<double> >
    (basic_srmatrix<double>& m,
    const basic_srmatrix<double>& mArg) throw (cvmexception)
{
    const int mnM = m.msize();
    if (mnM != mArg.msize()) throw cvmexception (CVM_SIZESMISMATCH);
    if (mnM == 1)
    {
        static const double one(1.);
        if (_abs (mArg(1,1)) <= basic_cvmMachMin<double>()) throw cvmexception (CVM_SINGULARMATRIX, 1);
        m(1,1) = one / mArg(1,1);
    }
    else
    {
        const int nWorkSize = mnM * 64;
        int nOutInfo  = 0;
        basic_array<int> nPivots (mnM);
        basic_rvector<double> vWork (nWorkSize);

        m.low_up (mArg, nPivots);
        DGETRI (&mnM, m, m._pld(), nPivots, vWork, &nWorkSize, &nOutInfo);
        if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
        if (nOutInfo > 0) throw cvmexception (CVM_SINGULARMATRIX, nOutInfo);
    }
}

template<>
CVM_API void
__polynom<float, basic_rvector<float> >
    (float* mpD, int ldP, int mnM, const float* pD, int ldA, const basic_rvector<float>& v)
{
    basic_rvector<float> vWork (NPOLY (&mnM, v._psize()));
    SPOLY (&mnM, pD, &ldA, v._psize(), v, mpD, &ldP, vWork);
}

template<>
CVM_API void
__polynom<double, basic_rvector<double> >
    (double* mpD, int ldP, int mnM, const double* pD, int ldA, const basic_rvector<double>& v)
{
    basic_rvector<double> vWork (NPOLY (&mnM, v._psize()));
    DPOLY (&mnM, pD, &ldA, v._psize(), v, mpD, &ldP, vWork);
}

// internal solver
// don't forget to make pX equal to pB before call
template<>
CVM_API void
__solve<float, float, basic_srmatrix<float> >
    (const basic_srmatrix<float>& m,
    int nrhs,
    const float* pB, int ldB, 
    float* pX, int ldX,
    float& dErr,
    const float* pLU, const int* pPivots) throw (cvmexception)
{
    const int mnM = m.msize();
    const bool bGivenLU = pLU != NULL && pPivots != NULL;
    int nOutInfo = 0;
    basic_rvector<float> vBerr (nrhs);
    basic_rvector<float> vFerr (nrhs);
    basic_rvector<float> vWork (3 * mnM);
    basic_array<int>   iWork (mnM);
    basic_array<int>   nPivots (mnM);

    if (bGivenLU) nPivots.assign (pPivots);
    basic_srmatrix<float> mLU (mnM);
    if (bGivenLU)
    {
        mLU.assign (pLU);
    }
    else
    {
        mLU = m.low_up (nPivots);
    }
    dErr = 0.F;

    SGETRS (Chars::pN(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
            1,
#endif
            &mnM, &nrhs, mLU, mLU._pld(), nPivots, pX, &ldX, &nOutInfo);
    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

    SGERFS (Chars::pN(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
            1,
#endif
            &mnM, &nrhs,
            m, m._pld(),
            mLU, mLU._pld(),
            nPivots,
            pB, &ldB,
            pX, &ldX,
            vFerr, vBerr, vWork, iWork, &nOutInfo);
    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

    dErr = vFerr.norminf();
}

template<>
CVM_API void
__solve<double, double, basic_srmatrix<double> >
    (const basic_srmatrix<double>& m,
    int nrhs,
    const double* pB, int ldB,
    double* pX, int ldX,
    double& dErr,
    const double* pLU, const int* pPivots) throw (cvmexception)
{
    const int mnM = m.msize();
    const bool bGivenLU = pLU != NULL && pPivots != NULL;
    int nOutInfo = 0;
    basic_rvector<double> vBerr (nrhs);
    basic_rvector<double> vFerr (nrhs);
    basic_rvector<double> vWork (3 * mnM);
    basic_array<int>   iWork (mnM);
    basic_array<int>   nPivots (mnM);

    if (bGivenLU) nPivots.assign (pPivots);
    basic_srmatrix<double> mLU (mnM);
    if (bGivenLU)
    {
        mLU.assign (pLU);
    }
    else
    {
        mLU = m.low_up(nPivots);
    }
    dErr = 0.;

    DGETRS (Chars::pN(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
            1,
#endif
            &mnM, &nrhs, mLU, mLU._pld(), nPivots, pX, &ldX, &nOutInfo);
    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

    DGERFS (Chars::pN(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
            1,
#endif
            &mnM, &nrhs,
            m, m._pld(),
            mLU, mLU._pld(),
            nPivots,
            pB, &ldB,
            pX, &ldX,
            vFerr, vBerr, vWork, iWork, &nOutInfo);
    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

    dErr = vFerr.norminf();
}

template<>
CVM_API void
__solve<float, float, basic_srbmatrix<float> >
    (const basic_srbmatrix<float>& m,
    int nrhs,
    const float* pB, int ldB,
    float* pX, int ldX,
    float& dErr, 
    const float* pLU, const int* pPivots) throw (cvmexception)
{
    const int mnM = m.msize();
    const int mnKL = m.lsize();
    const int mnKU = m.usize();
    const bool bGivenLU = pLU != NULL && pPivots != NULL;
    int nOutInfo = 0;
    basic_rvector<float> vBerr (nrhs);
    basic_rvector<float> vFerr (nrhs);
    basic_rvector<float> vWork (3 * mnM);
    basic_array<int>   iWork (mnM);
    basic_array<int>   nPivots (mnM);

    if (bGivenLU) nPivots.assign (pPivots);
    basic_srbmatrix<float> mLU (mnM, mnKL, mnKL + mnKU);
    if (bGivenLU)
    {
        mLU.assign (pLU);
    }
    else
    {
        mLU = m.low_up(nPivots);
    }
    dErr = 0.F;

    SGBTRS (Chars::pN(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
            1,
#endif
            &mnM, &mnKL, &mnKU, &nrhs, mLU, mLU._pld(), nPivots, pX, &ldX, &nOutInfo);
    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

    SGBRFS (Chars::pN(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
            1,
#endif
            &mnM, &mnKL, &mnKU, &nrhs,
            m, m._pld(),
            mLU, mLU._pld(),
            nPivots,
            pB, &ldB,
            pX, &ldX,
            vFerr, vBerr, vWork, iWork, &nOutInfo);
    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

    dErr = vFerr.norminf();
}

template<>
CVM_API void
__solve<double, double, basic_srbmatrix<double> >
    (const basic_srbmatrix<double>& m,
    int nrhs,
    const double* pB, int ldB,
    double* pX, int ldX,
    double& dErr,
    const double* pLU, const int* pPivots) throw (cvmexception)
{
    const int mnM = m.msize();
    const int mnKL = m.lsize();
    const int mnKU = m.usize();
    const bool bGivenLU = pLU != NULL && pPivots != NULL;
    int nOutInfo = 0;
    basic_rvector<double> vBerr (nrhs);
    basic_rvector<double> vFerr (nrhs);
    basic_rvector<double> vWork (3 * mnM);
    basic_array<int>   iWork (mnM);
    basic_array<int>   nPivots (mnM);

    if (bGivenLU) nPivots.assign (pPivots);
    basic_srbmatrix<double> mLU (mnM, mnKL, mnKL + mnKU);
    if (bGivenLU)
    {
        mLU.assign (pLU);
    }
    else
    {
        mLU = m.low_up(nPivots);
    }
    dErr = 0.;

    DGBTRS (Chars::pN(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
            1,
#endif
            &mnM, &mnKL, &mnKU, &nrhs, mLU, mLU._pld(), nPivots, pX, &ldX, &nOutInfo);
    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

    DGBRFS (Chars::pN(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
            1,
#endif
            &mnM, &mnKL, &mnKU, &nrhs,
            m, m._pld(),
            mLU, mLU._pld(),
            nPivots,
            pB, &ldB,
            pX, &ldX,
            vFerr, vBerr, vWork, iWork, &nOutInfo);
    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

    dErr = vFerr.norminf();
}

// internal solver
// don't forget to make pX equal to pB before call
template<>
CVM_API void
__solve<float, float, basic_srsmatrix<float> >
    (const basic_srsmatrix<float>& m,
    int nrhs,
    const float* pB, int ldB,
    float* pX, int ldX,
    float& dErr,
    const float* pLU, const int* pPivots) throw (cvmexception)
{
    const int nM = m.msize();
    const bool bCholeskyGiven = pLU != NULL && pPivots == NULL;        // no pivots means Cholesky
    const bool bBunchKaufmanGiven = pLU != NULL && pPivots != NULL;
    const bool bCalc = !bCholeskyGiven && !bBunchKaufmanGiven;
    bool bPositiveDefinite = bCholeskyGiven;

    int nOutInfo = 0;
    basic_rvector<float> vBerr (nrhs);
    basic_rvector<float> vFerr (nrhs);
    basic_rvector<float> vWork (3 * nM);
    basic_array<int>   iWork (nM);
    basic_array<int>   nPivots (nM);

    if (bBunchKaufmanGiven) nPivots.assign (pPivots);
    basic_srsmatrix<float> mLU(nM);
    if (bCalc)
    {
        mLU._factorize (m, nPivots, bPositiveDefinite);
    }
    else
    {
        mLU.assign (pLU);
    }
    dErr = 0.F;

    if (bPositiveDefinite)
    {
        SPOTRS (Chars::pU(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                1,
#endif
                &nM, &nrhs, mLU, mLU._pld(), pX, &ldX, &nOutInfo);
        if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

        SPORFS (Chars::pU(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                1,
#endif
                &nM, &nrhs, 
                m, m._pld(),
                mLU, mLU._pld(),
                pB, &ldB,
                pX, &ldX,
                vFerr, vBerr, 
                vWork, iWork, &nOutInfo);

        if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
    }
    else
    {
        SSYTRS (Chars::pU(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                1,
#endif
                &nM, &nrhs, mLU, mLU._pld(), nPivots, pX, &ldX, &nOutInfo);
        if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

        SSYRFS (Chars::pU(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                1,
#endif
                &nM, &nrhs,
                m, m._pld(),
                mLU, mLU._pld(),
                nPivots,
                pB, &ldB,
                pX, &ldX,
                vFerr, vBerr, 
                vWork, iWork, &nOutInfo);
        if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
    }

    dErr = vFerr.norminf();
}

template<>
CVM_API void
__solve<double, double, basic_srsmatrix<double> >
    (const basic_srsmatrix<double>& m,
    int nrhs,
    const double* pB, int ldB,
    double* pX, int ldX,
    double& dErr,
    const double* pLU, const int* pPivots) throw (cvmexception)
{
    const int nM = m.msize();
    const bool bCholeskyGiven = pLU != NULL && pPivots == NULL;        // no pivots means Cholesky
    const bool bBunchKaufmanGiven = pLU != NULL && pPivots != NULL;
    const bool bCalc = !bCholeskyGiven && !bBunchKaufmanGiven;
    bool bPositiveDefinite = bCholeskyGiven;

    int nOutInfo = 0;
    basic_rvector<double> vBerr (nrhs);
    basic_rvector<double> vFerr (nrhs);
    basic_rvector<double> vWork (3 * nM);
    basic_array<int>   iWork (nM);
    basic_array<int>   nPivots (nM);

    if (bBunchKaufmanGiven) nPivots.assign (pPivots);
    basic_srsmatrix<double> mLU(nM);
    if (bCalc)
    {
        mLU._factorize (m, nPivots, bPositiveDefinite);
    }
    else
    {
        mLU.assign (pLU);
    }
    dErr = 0.;

    if (bPositiveDefinite)
    {
        DPOTRS (Chars::pU(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                1,
#endif
                &nM, &nrhs, mLU, mLU._pld(), pX, &ldX, &nOutInfo);
        if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

        DPORFS (Chars::pU(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                1,
#endif
                &nM, &nrhs,
                m, m._pld(),
                mLU, mLU._pld(),
                pB, &ldB,
                pX, &ldX,
                vFerr, vBerr, 
                vWork, iWork, &nOutInfo);
        if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
    }
    else
    {
        DSYTRS (Chars::pU(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                1,
#endif
                &nM, &nrhs, mLU, mLU._pld(), nPivots, pX, &ldX, &nOutInfo);
        if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

        DSYRFS (Chars::pU(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                1,
#endif
                &nM, &nrhs,
                m, m._pld(),
                mLU, mLU._pld(),
                nPivots,
                pB, &ldB,
                pX, &ldX,
                vFerr, vBerr, 
                vWork, iWork, &nOutInfo);
        if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
    }

    dErr = vFerr.norminf();
}

template<>
CVM_API void
__inv<basic_srsmatrix<float> >
    (basic_srsmatrix<float>& m,
    const basic_srsmatrix<float>& mArg) throw (cvmexception)
{
    const int nM = m.msize();
    if (nM != mArg.msize()) throw cvmexception (CVM_SIZESMISMATCH);
    if (nM == 1)
    {
        static const float one(1.F);
        if (_abs (mArg(1,1)) <= basic_cvmMachMin<float>()) throw cvmexception (CVM_SINGULARMATRIX, 1);
        m.set (1, 1, one / mArg(1,1));
    }
    else
    {
        bool bPositiveDefinite = false;
        int nOutInfo = 0;
        basic_array<int> nPivots (nM);

        m._factorize (mArg, nPivots, bPositiveDefinite);

        if (bPositiveDefinite)
        {
            SPOTRI (Chars::pU(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                    1,
#endif
                    &nM, m, m._pld(), &nOutInfo);
            if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
            if (nOutInfo > 0) throw cvmexception (CVM_WRONGCHOLESKYFACTOR, nOutInfo);
        }
        else
        {
            basic_rvector<float> vWork (nM);
            SSYTRI (Chars::pU(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                    1,
#endif
                    &nM, m, m._pld(), nPivots, vWork, &nOutInfo);
            if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
            if (nOutInfo > 0) throw cvmexception (CVM_WRONGBUNCHKAUFMANFACTOR, nOutInfo);
        }
        m._flip();
    }
}

template<>
CVM_API void
__inv<basic_srsmatrix<double> >
    (basic_srsmatrix<double>& m,
    const basic_srsmatrix<double>& mArg) throw (cvmexception)
{
    const int nM = m.msize();
    if (nM != mArg.msize()) throw cvmexception (CVM_SIZESMISMATCH);
    if (nM == 1)
    {
        static const double one(1.);
        if (_abs (mArg(1,1)) <= basic_cvmMachMin<double>()) throw cvmexception (CVM_SINGULARMATRIX, 1);
        m.set (1, 1, one / mArg(1,1));
    }
    else
    {
        bool bPositiveDefinite = false;
        int nOutInfo = 0;
        basic_array<int> nPivots (nM);

        m._factorize (mArg, nPivots, bPositiveDefinite);

        if (bPositiveDefinite)
        {
            DPOTRI (Chars::pU(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                    1,
#endif
                    &nM, m, m._pld(), &nOutInfo);
            if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
            if (nOutInfo > 0) throw cvmexception (CVM_WRONGCHOLESKYFACTOR, nOutInfo);
        }
        else
        {
            basic_rvector<double> vWork (nM);
            DSYTRI (Chars::pU(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                    1,
#endif
                    &nM, m, m._pld(), nPivots, vWork, &nOutInfo);
            if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
            if (nOutInfo > 0) throw cvmexception (CVM_WRONGBUNCHKAUFMANFACTOR, nOutInfo);
        }
        m._flip();
    }
}

CVM_NAMESPACE_END
