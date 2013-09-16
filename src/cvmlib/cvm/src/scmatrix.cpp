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
__exp<basic_scmatrix<float, std::complex<float> >, float>
    (basic_scmatrix<float, std::complex<float> >& m,
    const basic_scmatrix<float, std::complex<float> >& mArg,
    float tol) throw (cvmexception)
{
    int nR = 0, nI = 0, nQ = 0, nJ = 0;
    const int mnM = m.msize();
    if (mnM != mArg.msize()) throw cvmexception (CVM_SIZESMISMATCH);

    basic_scmatrix<float, std::complex<float> > mTmp;
    const std::complex<float>* pD = mArg._pd();

    if (pD == m.get())
    {
        mTmp << mArg;
        pD = mTmp;
    }

    CMEXPC (&mnM, pD, pD == m.get() ? mTmp._pld() : mArg._pldm(), &tol, &nR, &nI, &nQ, &nJ);
    basic_cvector<float, std::complex<float> > vR (nR);
    basic_array<int> vI (nI);

    const int issymm = 0;
    std::complex<float> work_dummy(0.F);
    const int lwork_dummy = 0;

    CMEXP (&mnM, pD, pD == m.get() ? mTmp._pld() : mArg._pldm(), m, m._pld(),
           vR, vI, &nR, &nI, &nQ, &nJ, &issymm, &work_dummy, &lwork_dummy);
}

template<>
CVM_API void
__exp<basic_scmatrix<double, std::complex<double> >, double>
    (basic_scmatrix<double, std::complex<double> >& m,
    const basic_scmatrix<double, std::complex<double> >& mArg,
    double tol) throw (cvmexception)
{
    int nR = 0, nI = 0, nQ = 0, nJ = 0;
    const int mnM = m.msize();
    if (mnM != mArg.msize()) throw cvmexception (CVM_SIZESMISMATCH);

    basic_scmatrix<double, std::complex<double> > mTmp;
    const std::complex<double>* pD = mArg._pd();

    if (pD == m.get())
    {
        mTmp << mArg;
        pD = mTmp;
    }

    ZMEXPC (&mnM, pD, pD == m.get() ? mTmp._pld() : mArg._pldm(), &tol, &nR, &nI, &nQ, &nJ);
    basic_cvector<double, std::complex<double> > vR (nR);
    basic_array<int> vI (nI);

    const int issymm = 0;
    std::complex<double> work_dummy(0.);
    const int lwork_dummy = 0;

    ZMEXP (&mnM, pD, pD == m.get() ? mTmp._pld() : mArg._pldm(), m, m._pld(),
           vR, vI, &nR, &nI, &nQ, &nJ, &issymm, &work_dummy, &lwork_dummy);
}

template<>
CVM_API void
__exp_symm<basic_schmatrix<float, std::complex<float> >, float>
    (basic_schmatrix<float, std::complex<float> >& m, 
    const basic_schmatrix<float, std::complex<float> >& mArg, 
    float tol) throw (cvmexception)
{
    int nR = 0, nI = 0, nQ = 0, nJ = 0;
    const int nM = m.msize();
    if (nM != mArg.msize()) throw cvmexception (CVM_SIZESMISMATCH);

    basic_schmatrix<float, std::complex<float> > mTmp;
    const std::complex<float>* pD = mArg._pd();

    if (pD == m.get())
    {
        mTmp << mArg;
        pD = mTmp;
    }

    CMEXPC (&nM, pD, pD == m.get() ? mTmp._pld() : mArg._pldm(), &tol, &nR, &nI, &nQ, &nJ);
    basic_cvector<float, std::complex<float> > vR (nR);
    basic_array<int> vI (nI);

    const int ishem = 1;
    const int lwork  = 64 * nM;
    basic_cvector<float, std::complex<float> > work (lwork);

    CMEXP (&nM, pD, pD == m.get() ? mTmp._pld() : mArg._pldm(), m, m._pld(),
           vR, vI, &nR, &nI, &nQ, &nJ, &ishem, work, &lwork);
}

template<>
CVM_API void
__exp_symm<basic_schmatrix<double, std::complex<double> >, double>
    (basic_schmatrix<double, std::complex<double> >& m, 
    const basic_schmatrix<double, std::complex<double> >& mArg, 
    double tol) throw (cvmexception)
{
    int nR = 0, nI = 0, nQ = 0, nJ = 0;
    const int nM = m.msize();
    if (nM != mArg.msize()) throw cvmexception (CVM_SIZESMISMATCH);

    basic_schmatrix<double, std::complex<double> > mTmp;
    const std::complex<double>* pD = mArg._pd();

    if (pD == m.get())
    {
        mTmp << mArg;
        pD = mTmp;
    }

    ZMEXPC (&nM, pD, pD == m.get() ? mTmp._pld() : mArg._pldm(), &tol, &nR, &nI, &nQ, &nJ);
    basic_cvector<double, std::complex<double> > vR (nR);
    basic_array<int> vI (nI);

    const int ishem = 1;
    const int lwork  = 64 * nM;
    basic_cvector<double, std::complex<double> > work (lwork);

    ZMEXP (&nM, pD, pD == m.get() ? mTmp._pld() : mArg._pldm(), m, m._pld(),
           vR, vI, &nR, &nI, &nQ, &nJ, &ishem, work, &lwork);
}

template<>
CVM_API void
__cond_num<float, basic_scmatrix<float, std::complex<float> > >
    (const basic_scmatrix<float, std::complex<float> >& mArg, float& dCond) throw (cvmexception)
{
    dCond = 0.F;
    const int mnM  = mArg.msize();
    int   nOutInfo = 0;
    basic_scmatrix<float, std::complex<float> > mA (mArg);
    basic_cvector<float, std::complex<float> > work (mnM * 2);
    basic_rvector<float> rwork (mnM * 2);
    basic_array<int> iwork (mnM);

    const float rNorm = mA.norminf();
    CGETRF (&mnM, &mnM, mA, mA._pld(), iwork, &nOutInfo);

    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
    if (nOutInfo == 0)
    {
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
        CGECON (Chars::pI(), 1,
#else
        CGECON (Chars::pI(),
#endif
                &mnM, mA, mA._pld(), &rNorm, &dCond, work, rwork, &nOutInfo);
    }
}

template<>
CVM_API void
__cond_num<double, basic_scmatrix<double, std::complex<double> > >
    (const basic_scmatrix<double, std::complex<double> >& mArg, double& dCond) throw (cvmexception)
{
    dCond = 0.;
    const int mnM = mArg.msize();
    int nOutInfo  = 0;
    basic_scmatrix<double, std::complex<double> > mA (mArg);
    basic_cvector<double, std::complex<double> > work (mnM * 2);
    basic_rvector<double> rwork (mnM * 2);
    basic_array<int> iwork (mnM);

    const double rNorm = mA.norminf();
    ZGETRF (&mnM, &mnM, mA, mA._pld(), iwork, &nOutInfo);

    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
    if (nOutInfo == 0)
    {
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
        ZGECON (Chars::pI(), 1,
#else
        ZGECON (Chars::pI(),
#endif
                &mnM, mA, mA._pld(), &rNorm, &dCond, work, rwork, &nOutInfo);
    }
}

template <>
CVM_API void
__inv<basic_scmatrix<float, std::complex<float> > >
    (basic_scmatrix<float, std::complex<float> >& m,
    const basic_scmatrix<float, std::complex<float> >& mArg) throw (cvmexception)
{
    const int mnM = m.msize();
    if (mnM != mArg.msize()) throw cvmexception (CVM_SIZESMISMATCH);

    if (mnM == 1)
    {
        static const std::complex<float> one(1.F);
        if (_abs (mArg(1,1)) <= basic_cvmMachMin<float>()) throw cvmexception (CVM_SINGULARMATRIX, 1);
        m(1,1) = one / mArg(1,1);
    }
    else
    {
        const int nWorkSize = mnM * 64;
        int nOutInfo = 0;
        basic_array<int> nPivots (mnM);
        basic_cvector<float, std::complex<float> > vWork (nWorkSize);

        m.low_up (mArg, nPivots);
        CGETRI (&mnM, m, m._pld(), nPivots, vWork, &nWorkSize, &nOutInfo);

        if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
        if (nOutInfo > 0) throw cvmexception (CVM_SINGULARMATRIX, nOutInfo);
    }
}

template<>
CVM_API void
__inv<basic_scmatrix<double, std::complex<double> > >
    (basic_scmatrix<double, std::complex<double> >& m,
    const basic_scmatrix<double, std::complex<double> >& mArg) throw (cvmexception)
{
    const int mnM = m.msize();
    if (mnM != mArg.msize()) throw cvmexception (CVM_SIZESMISMATCH);

    if (mnM == 1)
    {
        static const std::complex<double> one(1.);
        if (_abs (mArg(1,1)) <= basic_cvmMachMin<double>()) throw cvmexception (CVM_SINGULARMATRIX, 1);
        m(1,1) = one / mArg(1,1);
    }
    else
    {
        const int nWorkSize = mnM * 64;
        int nOutInfo  = 0;
        basic_array<int> nPivots (mnM);
        basic_cvector<double, std::complex<double> > vWork (nWorkSize);

        m.low_up (mArg, nPivots);
        ZGETRI (&mnM, m, m._pld(), nPivots, vWork, &nWorkSize, &nOutInfo);

        if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
        if (nOutInfo > 0) throw cvmexception (CVM_SINGULARMATRIX, nOutInfo);
    }
}

template<>
CVM_API void
__inv<basic_schmatrix<float, std::complex<float> > >
    (basic_schmatrix<float, std::complex<float> >& m,
     const basic_schmatrix<float, std::complex<float> >& mArg) throw (cvmexception)
{
    const int nM = m.msize();
    if (nM != mArg.msize()) throw cvmexception (CVM_SIZESMISMATCH);
    if (nM == 1)
    {
        static const std::complex<float> one(1.F);
        if (_abs (mArg(1,1)) <= basic_cvmMachMin<float>()) throw cvmexception (CVM_SINGULARMATRIX, 1);
        m(1,1) = one / mArg(1,1);
    }
    else
    {
        bool bPositiveDefinite = false;
        int nOutInfo = 0;
        basic_array<int> nPivots (nM);

        m._factorize (mArg, nPivots, bPositiveDefinite);

        if (bPositiveDefinite)
        {
            CPOTRI (Chars::pU(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                    1,
#endif
                    &nM, m, m._pld(), &nOutInfo);
            if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
            if (nOutInfo > 0) throw cvmexception (CVM_WRONGCHOLESKYFACTOR, nOutInfo);
        }
        else
        {
            basic_cvector<float, std::complex<float> > vWork (nM);
            CHETRI (Chars::pU(),
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
__inv<basic_schmatrix<double, std::complex<double> > >
    (basic_schmatrix<double, std::complex<double> >& m,
     const basic_schmatrix<double, std::complex<double> >& mArg) throw (cvmexception)
{
    const int nM = m.msize();
    if (nM != mArg.msize()) throw cvmexception (CVM_SIZESMISMATCH);
    if (nM == 1)
    {
        static const std::complex<double> one(1.);
        if (_abs (mArg(1,1)) <= basic_cvmMachMin<double>()) throw cvmexception (CVM_SINGULARMATRIX, 1);
        m(1,1) = one / mArg(1,1);
    }
    else
    {
        bool bPositiveDefinite = false;
        int nOutInfo = 0;
        basic_array<int> nPivots (nM);

        m._factorize (mArg, nPivots, bPositiveDefinite);

        if (bPositiveDefinite)
        {
            ZPOTRI (Chars::pU(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                    1,
#endif
                    &nM, m, m._pld(), &nOutInfo);
            if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
            if (nOutInfo > 0) throw cvmexception (CVM_WRONGCHOLESKYFACTOR, nOutInfo);
        }
        else
        {
            basic_cvector<double, std::complex<double> > vWork (nM);
            ZHETRI (Chars::pU(),
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
__polynom<std::complex<float>, basic_cvector<float, std::complex<float> > >
    (std::complex<float>* mpD, int ldP,
    int mnM,
    const std::complex<float>* pD, int ldA,
    const basic_cvector<float, std::complex<float> >& v)
{
    basic_cvector<float, std::complex<float> > vWork (NPOLY (&mnM, v._psize()));
    CPOLY (&mnM, pD, &ldA, v._psize(), v, mpD, &ldP, vWork);
}

template<>
CVM_API void
__polynom<std::complex<double>, basic_cvector<double, std::complex<double> > >
    (std::complex<double>* mpD, int ldP,
    int mnM,
    const std::complex<double>* pD, int ldA,
    const basic_cvector<double, std::complex<double> >& v)
{
    basic_cvector<double, std::complex<double> > vWork (NPOLY (&mnM, v._psize()));
    ZPOLY (&mnM, pD, &ldA, v._psize(), v, mpD, &ldP, vWork);
}

// internal solver
// do not forget to make pX equal to pB before call
// in case of vectors passed increment MUST be 1 (since those methods assume matrices)
template<>
CVM_API void
__solve<float, std::complex<float>, basic_scmatrix<float, std::complex<float> > >
    (const basic_scmatrix<float, std::complex<float> >& m,
    int nrhs, 
    const std::complex<float>* pB, int ldB,
    std::complex<float>* pX, int ldX,
    float& dErr,
    const std::complex<float>* pLU, const int* pPivots) throw (cvmexception)
{
    const int mnM = m.msize();
    const bool bGivenLU = pLU != NULL && pPivots != NULL;
    int nOutInfo = 0;
    basic_rvector<float> vFerr (nrhs);
    basic_rvector<float> vBerr (nrhs);
    basic_cvector<float, std::complex<float> > vWork (2 * mnM);
    basic_rvector<float> rWork (mnM);
    basic_array<int> nPivots (mnM);

    if (bGivenLU) nPivots.assign (pPivots);
    basic_scmatrix<float, std::complex<float> > mLU (mnM);
    if (bGivenLU)
    {
        mLU.assign (pLU);
    }
    else
    {
        mLU = m.low_up (nPivots);
    }
    dErr = 0.F;

#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
    CGETRS (Chars::pN(), 1,
#else
    CGETRS (Chars::pN(),
#endif
            &mnM, &nrhs, mLU, mLU._pld(), nPivots, pX, &ldX, &nOutInfo);

    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
    CGERFS (Chars::pN(), 1,
#else
    CGERFS (Chars::pN(),
#endif
            &mnM, &nrhs,
            m, m._pld(),
            mLU, mLU._pld(),
            nPivots,
            pB, &ldB,
            pX, &ldX,
            vFerr, vBerr, vWork, rWork, &nOutInfo);

    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
    dErr = vFerr.norminf();
}

template<>
CVM_API void
__solve<double, std::complex<double>, basic_scmatrix<double, std::complex<double> > >
    (const basic_scmatrix<double, std::complex<double> >& m,
    int nrhs,
    const std::complex<double>* pB, int ldB,
    std::complex<double>* pX, int ldX,
    double& dErr,
    const std::complex<double>* pLU, const int* pPivots) throw (cvmexception)
{
    const int mnM = m.msize();
    const bool bGivenLU = pLU != NULL && pPivots != NULL;
    int nOutInfo = 0;
    basic_rvector<double> vBerr (nrhs);
    basic_rvector<double> vFerr (nrhs);
    basic_cvector<double, std::complex<double> > vWork (2 * mnM);
    basic_rvector<double> rWork (mnM);
    basic_array<int> nPivots (mnM);

    if (bGivenLU) nPivots.assign (pPivots);
    basic_scmatrix<double, std::complex<double> > mLU (mnM);
    if (bGivenLU)
    {
        mLU.assign (pLU);
    }
    else
    {
        mLU = m.low_up (nPivots);
    }
    dErr = 0.;

#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
    ZGETRS (Chars::pN(), 1,
#else
    ZGETRS (Chars::pN(),
#endif
            &mnM, &nrhs, mLU, mLU._pld(), nPivots, pX, &ldX, &nOutInfo);

    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
    ZGERFS (Chars::pN(), 1,
#else
    ZGERFS (Chars::pN(),
#endif
            &mnM, &nrhs,
            m, m._pld(),
            mLU, mLU._pld(),
            nPivots,
            pB, &ldB,
            pX, &ldX,
            vFerr, vBerr, vWork, rWork, &nOutInfo);

    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
    dErr = vFerr.norminf();
}

template<>
CVM_API void
__solve<float, std::complex<float>, basic_scbmatrix<float, std::complex<float> > >
    (const basic_scbmatrix<float, std::complex<float> >& m,
    int nrhs,
    const std::complex<float>* pB, int ldB,
    std::complex<float>* pX, int ldX,
    float& dErr,
    const std::complex<float>* pLU, const int* pPivots) throw (cvmexception)
{
    const int mnM = m.msize();
    const int mnKL = m.lsize();
    const int mnKU = m.usize();
    const bool bGivenLU = pLU != NULL && pPivots != NULL;
    int nOutInfo = 0;
    basic_rvector<float> vFerr (nrhs);
    basic_rvector<float> vBerr (nrhs);
    basic_cvector<float, std::complex<float> > vWork (2 * mnM);
    basic_rvector<float> rWork (mnM);
    basic_array<int>nPivots (mnM);

    if (bGivenLU) nPivots.assign (pPivots);
    basic_scbmatrix<float, std::complex<float> > mLU (mnM, mnKL, mnKL + mnKU);
    if (bGivenLU)
    {
        mLU.assign (pLU);
    }
    else
    {
        mLU = m.low_up (nPivots);
    }
    dErr = 0.F;

    CGBTRS (Chars::pN(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
            1,
#endif
            &mnM, &mnKL, &mnKU, &nrhs, mLU, mLU._pld(), nPivots, pX, &ldX, &nOutInfo);

    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

    CGBRFS (Chars::pN(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
            1,
#endif
            &mnM, &mnKL, &mnKU, &nrhs,
            m, m._pld(),
            mLU, mLU._pld(),
            nPivots,
            pB, &ldB,
            pX, &ldX,
            vFerr, vBerr, vWork, rWork, &nOutInfo);

    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
    dErr = vFerr.norminf();
}

template<>
CVM_API void
__solve<double, std::complex<double>, basic_scbmatrix<double, std::complex<double> > >
    (const basic_scbmatrix<double, std::complex<double> >& m,
    int nrhs,
    const std::complex<double>* pB, int ldB,
    std::complex<double>* pX, int ldX,
    double& dErr,
    const std::complex<double>* pLU, const int* pPivots) throw (cvmexception)
{
    const int mnM = m.msize();
    const int mnKL = m.lsize();
    const int mnKU = m.usize();
    const bool bGivenLU = pLU != NULL && pPivots != NULL;
    int nOutInfo = 0;
    basic_rvector<double> vFerr (nrhs);
    basic_rvector<double> vBerr (nrhs);
    basic_cvector<double, std::complex<double> > vWork (2 * mnM);
    basic_rvector<double> rWork (mnM);
    basic_array<int> nPivots (mnM);

    if (bGivenLU) nPivots.assign (pPivots);
    basic_scbmatrix<double, std::complex<double> > mLU (mnM, mnKL, mnKL + mnKU);
    if (bGivenLU)
    {
        mLU.assign (pLU);
    }
    else
    {
        mLU = m.low_up (nPivots);
    }
    dErr = 0.;

    ZGBTRS (Chars::pN(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
            1,
#endif
            &mnM, &mnKL, &mnKU, &nrhs, mLU, mLU._pld(), nPivots, pX, &ldX, &nOutInfo);

    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

    ZGBRFS (Chars::pN(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
            1,
#endif
            &mnM, &mnKL, &mnKU, &nrhs,
            m, m._pld(),
            mLU, mLU._pld(),
            nPivots,
            pB, &ldB,
            pX, &ldB,
            vFerr, vBerr, vWork, rWork, &nOutInfo);

    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
    dErr = vFerr.norminf();
}


// internal solver
// don't forget to make pX equal to pB before call
template<>
CVM_API void
__solve<float, std::complex<float>, basic_schmatrix<float, std::complex<float> > >
    (const basic_schmatrix<float, std::complex<float> >& m,
    int nrhs,
    const std::complex<float>* pB, int ldB,
    std::complex<float>* pX, int ldX,
    float& dErr,
    const std::complex<float>* pLU, const int* pPivots) throw (cvmexception)
{
    const int nM = m.msize();
    const bool bCholeskyGiven = pLU != NULL && pPivots == NULL;        // no pivots means Cholesky
    const bool bBunchKaufmanGiven = pLU != NULL && pPivots != NULL;
    const bool bCalc = !bCholeskyGiven && !bBunchKaufmanGiven;
    bool bPositiveDefinite = bCholeskyGiven;

    int nOutInfo = 0;
    basic_rvector<float> vBerr (nrhs);
    basic_rvector<float> vFerr (nrhs);
    basic_cvector<float, std::complex<float> > vWork (2 * nM);
    basic_rvector<float> vrWork (nM);
    basic_array<int> nPivots (nM);

    if (bBunchKaufmanGiven) nPivots.assign (pPivots);
    basic_schmatrix<float, std::complex<float> > mLU(nM);
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
        CPOTRS (Chars::pU(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                1,
#endif
                &nM, &nrhs, mLU, mLU._pld(), pX, &ldX, &nOutInfo);
        if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

        CPORFS (Chars::pU(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                1,
#endif
                &nM, &nrhs, 
                m, m._pld(),
                mLU, mLU._pld(),
                pB, &ldB,
                pX, &ldX,
                vFerr, vBerr, 
                vWork, vrWork, &nOutInfo);
        if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
    }
    else
    {
        CHETRS (Chars::pU(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                1,
#endif
                &nM, &nrhs, mLU, mLU._pld(), nPivots, pX, &ldX, &nOutInfo);
        if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

        CHERFS (Chars::pU(),
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
                vWork, vrWork, &nOutInfo);
        if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
    }

    dErr = vFerr.norminf();
}

template<>
CVM_API void
__solve<double, std::complex<double>, basic_schmatrix<double, std::complex<double> > >
    (const basic_schmatrix<double, std::complex<double> >& m,
    int nrhs,
    const std::complex<double>* pB, int ldB,
    std::complex<double>* pX, int ldX,
    double& dErr,
    const std::complex<double>* pLU, const int* pPivots) throw (cvmexception)
{
    const int nM = m.msize();
    const bool bCholeskyGiven = pLU != NULL && pPivots == NULL;        // no pivots means Cholesky
    const bool bBunchKaufmanGiven = pLU != NULL && pPivots != NULL;
    const bool bCalc = !bCholeskyGiven && !bBunchKaufmanGiven;
    bool bPositiveDefinite = bCholeskyGiven;

    int nOutInfo = 0;
    basic_rvector<double> vBerr (nrhs);
    basic_rvector<double> vFerr (nrhs);
    basic_cvector<double, std::complex<double> > vWork (2 * nM);
    basic_rvector<double> vrWork (nM);
    basic_array<int> nPivots (nM);

    if (bBunchKaufmanGiven) nPivots.assign (pPivots);
    basic_schmatrix<double, std::complex<double> > mLU(nM);
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
        ZPOTRS (Chars::pU(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                1,
#endif
                &nM, &nrhs, mLU, mLU._pld(), pX, &ldX, &nOutInfo);
        if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

        ZPORFS (Chars::pU(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                1,
#endif
                &nM, &nrhs, 
                m, m._pld(),
                mLU, mLU._pld(),
                pB, &ldB,
                pX, &ldX,
                vFerr, vBerr, 
                vWork, vrWork, &nOutInfo);
        if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
    }
    else
    {
        ZHETRS (Chars::pU(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                1,
#endif
                &nM, &nrhs, mLU, mLU._pld(), nPivots, pX, &ldX, &nOutInfo);
        if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

        ZHERFS (Chars::pU(),
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
                vWork, vrWork, &nOutInfo);
        if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
    }

    dErr = vFerr.norminf();
}

CVM_NAMESPACE_END
