/* CVM Class Library */
/* Copyright (C), Sergei Nikolaev, 1992-2006, http://cvmlib.com */

#include "cvm.h"
#include "blas.h"

#if defined (_MSC_VER)
#   pragma warning(disable:4786)
#endif

CVM_NAMESPACE_BEG

template<>
CVM_API float __dot<float> (const float* mpD, int mnSize, int mnIncr, const float* pD, int nIncr)
{
    return SDOT (&mnSize, mpD, &mnIncr, pD, &nIncr);
}

template<>
CVM_API double __dot<double> (const double* mpD, int mnSize, int mnIncr, const double* pD, int nIncr)
{
    return DDOT (&mnSize, mpD, &mnIncr, pD, &nIncr);
}

template<>
CVM_API void
__gemv<float, basic_rmatrix<float>, basic_rvector<float> >
    (bool bLeft, 
    const basic_rmatrix<float>& m,
    float dAlpha,
    const basic_rvector<float>& v,
    float dBeta,
    basic_rvector<float>& vRes)
{
    SGEMV (bLeft ? Chars::pT() : Chars::pN(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
           1,
#endif
           m._pm(), m._pn(), &dAlpha, m._pd(), m._pldm(), v, v._pincr(), &dBeta, vRes, vRes._pincr());
}

template<>
CVM_API void
__gemv<double, basic_rmatrix<double>, basic_rvector<double> >
    (bool bLeft,
    const basic_rmatrix<double>& m,
    double dAlpha,
    const basic_rvector<double>& v,
    double dBeta,
    basic_rvector<double>& vRes)
{
    DGEMV (bLeft ? Chars::pT() : Chars::pN(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
           1,
#endif
           m._pm(), m._pn(), &dAlpha, m._pd(), m._pldm(), v, v._pincr(), &dBeta, vRes, vRes._pincr());
}

template<>
CVM_API void
__gbmv<float, basic_srbmatrix<float>, basic_rvector<float> >
    (bool bLeft,
    const basic_srbmatrix<float>& m,
    float dAlpha,
    const basic_rvector<float>& v,
    float dBeta,
    basic_rvector<float>& vRes)
{
    SGBMV (bLeft ? Chars::pT() : Chars::pN(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
           1,
#endif
           m._pm(), m._pn(), m._pl(), m._pu(), &dAlpha, m, m._pld(), v, v._pincr(), &dBeta, vRes, vRes._pincr());
}

template<>
CVM_API void
__gbmv<double, basic_srbmatrix<double>, basic_rvector<double> >
    (bool bLeft,
    const basic_srbmatrix<double>& m,
    double dAlpha,
    const basic_rvector<double>& v,
    double dBeta,
    basic_rvector<double>& vRes)
{
    DGBMV (bLeft ? Chars::pT() : Chars::pN(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
           1,
#endif
           m._pm(), m._pn(), m._pl(), m._pu(), &dAlpha, m, m._pld(), v, v._pincr(), &dBeta, vRes, vRes._pincr());
}

template<>
CVM_API void
__symv<float, basic_srsmatrix<float>, basic_rvector<float> >
    (const basic_srsmatrix<float>& m,
    float dAlpha,
    const basic_rvector<float>& v,
    float dBeta,
    basic_rvector<float>& vRes)
{
    SSYMV (Chars::pU(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
           1,
#endif
           m._pm(), &dAlpha, m, m._pld(), v, v._pincr(), &dBeta, vRes, vRes._pincr());
}

template<>
CVM_API void
__symv<double, basic_srsmatrix<double>, basic_rvector<double> >
    (const basic_srsmatrix<double>& m,
    double dAlpha,
    const basic_rvector<double>& v,
    double dBeta,
    basic_rvector<double>& vRes)
{
    DSYMV (Chars::pU(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
           1,
#endif
           m._pm(), &dAlpha, m, m._pld(), v, v._pincr(), &dBeta, vRes, vRes._pincr());
}

template<>
CVM_API void
__svd<float, basic_rmatrix<float>, basic_srmatrix<float> >
    (float* pD, int nSize, int nIncr, 
    const basic_rmatrix<float>& mArg,
    basic_srmatrix<float>* mU,
    basic_srmatrix<float>* mVH) throw (cvmexception)
{
    const bool bSimple  = (mU == NULL || mVH == NULL);
    const int  nM       = mArg.msize();
    const int  nN       = mArg.nsize();
    const int  m        = _cvm_min<int>(nM, nN);
    const int  M        = _cvm_max<int>(nM, nN);
    int lWork = -1; // to calculate lWork
    int nOutInfo = 0;

    if (m != nSize) throw cvmexception (CVM_SIZESMISMATCH);

    basic_rvector<float> mD (nSize);
    mD.assign (pD, nIncr);

    basic_rvector<float> vOffDiag (_cvm_max<int>(1, m - 1));
    basic_rvector<float> vTauQ    (m);
    basic_rvector<float> vTauP    (m);
    basic_rmatrix<float> mA       (mArg);
    float dWork;

    SGEBRD (&nM, &nN, mA, mA._pld(), mD, vOffDiag, vTauQ, vTauP, &dWork, &lWork, &nOutInfo);
    lWork = static_cast<int>(dWork);
    basic_rvector<float> vWork (lWork);

    SGEBRD (&nM, &nN, mA, mA._pld(), mD, vOffDiag, vTauQ, vTauP, vWork, &lWork, &nOutInfo);
    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

    static const int zero(0);
    static const int one(1);
    if (bSimple)
    {
        basic_rvector<float> vWork2 (m * 4);
        SBDSQR (nM >= nN ? Chars::pU() : Chars::pL(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                1,
#endif
                &m,
                &zero, &zero, &zero,
                mD, vOffDiag,
                NULL, &one, NULL, &one, NULL, &one,
                vWork2, &nOutInfo);
    }
    else
    {
        const int lWork3 = m * _cvm_max<int>(M, 64);    // bug fix for cases when N>64*M or M>64*N
        basic_rvector<float> vWork2 (m * 4);
        basic_rvector<float> vWork3 (lWork3);

        basic_rmatrix<float> Q (mA);
        basic_rmatrix<float> P (mA);

        if (nM > nN) Q.resize(nM, nM);
        if (nM < nN) P.resize(nN, nN);

        SORGBR (Chars::pQ(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                1,
#endif
                &nM, &nM, &nN,
                Q, &nM,
                vTauQ,
                vWork3, &lWork3, &nOutInfo);
        if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

        SORGBR (Chars::pP(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                1,
#endif
                &nN, &nN, &nM,
                P, &M,
                vTauP,
                vWork3, &lWork3, &nOutInfo);
        if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

        SBDSQR (nM >= nN ? Chars::pU() : Chars::pL(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                1,
#endif
                &m,
                &nN, &nM,
                &zero,
                mD, vOffDiag,
                P, &M, Q, &nM, NULL, &one,
                vWork2, &nOutInfo);

        (*mU)  = basic_srmatrix<float>(Q.resize(nM, nM));
        (*mVH) = basic_srmatrix<float>(P.resize(nN, nN));
    }

    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
    if (nOutInfo > 0) throw cvmexception (CVM_CONVERGENCE_ERROR);

    __copy<float> (nSize, mD, mD.incr(), pD, nIncr);
}

template<>
CVM_API void
__svd<double, basic_rmatrix<double>, basic_srmatrix<double> >
    (double* pD, int nSize, int nIncr, 
    const basic_rmatrix<double>& mArg,
    basic_srmatrix<double>* mU,
    basic_srmatrix<double>* mVH) throw (cvmexception)
{
    const bool bSimple  = (mU == NULL || mVH == NULL);
    const int  nM       = mArg.msize();
    const int  nN       = mArg.nsize();
    const int  m        = _cvm_min<int>(nM, nN);
    const int  M        = _cvm_max<int>(nM, nN);
    int lWork = -1; // to calculate lWork
    int nOutInfo = 0;

    if (m != nSize) throw cvmexception (CVM_SIZESMISMATCH);

    basic_rvector<double> mD (nSize);
    mD.assign (pD, nIncr);

    basic_rvector<double> vOffDiag (_cvm_max<int>(1, m - 1));
    basic_rvector<double> vTauQ    (m);
    basic_rvector<double> vTauP    (m);
    basic_rmatrix<double> mA       (mArg);
    double dWork;

    DGEBRD (&nM, &nN, mA, mA._pld(), mD, vOffDiag, vTauQ, vTauP, &dWork, &lWork, &nOutInfo);
    lWork = static_cast<int>(dWork);
    basic_rvector<double> vWork(lWork);

    DGEBRD (&nM, &nN, mA, mA._pld(), mD, vOffDiag, vTauQ, vTauP, vWork, &lWork, &nOutInfo);
    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

    static const int zero(0);
    static const int one(1);
    if (bSimple)
    {
        basic_rvector<double> vWork2 (m * 4);
        DBDSQR (nM >= nN ? Chars::pU() : Chars::pL(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                1,
#endif
                &m,
                &zero, &zero, &zero,
                mD, vOffDiag,
                NULL, &one, NULL, &one, NULL, &one,
                vWork2, &nOutInfo);
    }
    else
    {
        const int lWork3 = m * _cvm_max<int>(M, 64);    // bug fix for cases when N>64*M or M>64*N
        basic_rvector<double> vWork2 (m * 4);
        basic_rvector<double> vWork3 (lWork3);

        basic_rmatrix<double> Q (mA);
        basic_rmatrix<double> P (mA);

        if (nM > nN) Q.resize(nM, nM);
        if (nM < nN) P.resize(nN, nN);

        DORGBR (Chars::pQ(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                1,
#endif
                &nM, &nM, &nN,
                Q, &nM,
                vTauQ,
                vWork3, &lWork3, &nOutInfo);
        if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

        DORGBR (Chars::pP(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                1,
#endif
                &nN, &nN, &nM,
                P, &M,
                vTauP,
                vWork3, &lWork3, &nOutInfo);
        if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

        DBDSQR (nM >= nN ? Chars::pU() : Chars::pL(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                1,
#endif
                &m,
                &nN, &nM,
                &zero,
                mD, vOffDiag,
                P, &M, Q, &nM, NULL, &one,
                vWork2, &nOutInfo);

        (*mU)  = basic_srmatrix<double>(Q.resize(nM, nM));
        (*mVH) = basic_srmatrix<double>(P.resize(nN, nN));
    }

    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
    if (nOutInfo > 0) throw cvmexception (CVM_CONVERGENCE_ERROR);

    __copy<double> (nSize, mD, mD.incr(), pD, nIncr);
}

template<>
CVM_API void
__svd<float, basic_cmatrix<float, std::complex<float> >, basic_scmatrix<float, std::complex<float> > >
    (float* pD, int nSize, int nIncr,
    const basic_cmatrix<float, std::complex<float> >& mArg,
    basic_scmatrix<float, std::complex<float> >* mU,
    basic_scmatrix<float, std::complex<float> >* mVH) throw (cvmexception)
{
    const bool bSimple  = (mU == NULL || mVH == NULL);
    const int  nM       = mArg.msize();
    const int  nN       = mArg.nsize();
    const int  m        = _cvm_min<int>(nM, nN);
    const int  M        = _cvm_max<int>(nM, nN);
    int lWork = -1; // to calculate lWork
    int nOutInfo = 0;

    if (m != nSize) throw cvmexception (CVM_SIZESMISMATCH);

    basic_rvector<float> mD (nSize);
    mD.assign (pD, nIncr);

    basic_rvector<float> vOffDiag (_cvm_max<int>(1, m - 1));
    basic_cvector<float, std::complex<float> > vTauQ (m);
    basic_cvector<float, std::complex<float> > vTauP (m);
    basic_cmatrix<float, std::complex<float> > mA (mArg);
    std::complex<float> dWork;

    CGEBRD (&nM, &nN, mA, mA._pld(), mD, vOffDiag, vTauQ, vTauP, &dWork, &lWork, &nOutInfo);
    lWork = static_cast<int>(dWork.real());
    basic_cvector<float, std::complex<float> > vWork (lWork);

    CGEBRD (&nM, &nN, mA, mA._pld(), mD, vOffDiag, vTauQ, vTauP, vWork, &lWork, &nOutInfo);
    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

    static const int zero(0);
    static const int one(1);
    if (bSimple)
    {
        basic_rvector<float> vWork2 (m * 4);
        CBDSQR (nM >= nN ? Chars::pU() : Chars::pL(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                1,
#endif
                &m,
                &zero, &zero, &zero,
                mD, vOffDiag,
                NULL, &one, NULL, &one, NULL, &one,
                vWork2, &nOutInfo);
    }
    else
    {
        const int lWork3 = m * _cvm_max<int>(M, 64);    // bug fix for cases when N>64*M or M>64*N
        basic_rvector<float> vWork2 (m * 4);
        basic_cvector<float, std::complex<float> > vWork3 (lWork3);

        basic_cmatrix<float, std::complex<float> > Q (mA);
        basic_cmatrix<float, std::complex<float> > P (mA);

        if (nM > nN) Q.resize(nM, nM);
        if (nM < nN) P.resize(nN, nN);

        CUNGBR (Chars::pQ(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                1,
#endif
                &nM, &nM, &nN,
                Q, &nM,
                vTauQ,
                vWork3, &lWork3, &nOutInfo);
        if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

        CUNGBR (Chars::pP(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                1,
#endif
                &nN, &nN, &nM,
                P, &M,
                vTauP,
                vWork3, &lWork3, &nOutInfo);
        if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

        CBDSQR (nM >= nN ? Chars::pU() : Chars::pL(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                1,
#endif
                &m,
                &nN, &nM, 
                &zero,
                mD, vOffDiag,
                P, &M, Q, &nM, NULL, &one,
                vWork2, &nOutInfo);

        (*mU)  = basic_scmatrix<float, std::complex<float> >(Q.resize(nM, nM));
        (*mVH) = basic_scmatrix<float, std::complex<float> >(P.resize(nN, nN));
    }

    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
    if (nOutInfo > 0) throw cvmexception (CVM_CONVERGENCE_ERROR);

    __copy<float> (nSize, mD, mD.incr(), pD, nIncr);
}

template<>
CVM_API void
__svd<double, basic_cmatrix<double, std::complex<double> >, basic_scmatrix<double, std::complex<double> > >
    (double* pD, int nSize, int nIncr,
    const basic_cmatrix<double, std::complex<double> >& mArg,
    basic_scmatrix<double, std::complex<double> >* mU, 
    basic_scmatrix<double, std::complex<double> >* mVH) throw (cvmexception)
{
    const bool bSimple  = (mU == NULL || mVH == NULL);
    const int  nM       = mArg.msize();
    const int  nN       = mArg.nsize();
    const int  m        = _cvm_min<int>(nM, nN);
    const int  M        = _cvm_max<int>(nM, nN);
    int lWork = -1; // to calculate lWork
    int nOutInfo = 0;

    if (m != nSize) throw cvmexception (CVM_SIZESMISMATCH);

    basic_rvector<double> mD (nSize);
    mD.assign (pD, nIncr);

    basic_rvector<double> vOffDiag (_cvm_max<int>(1, m - 1));
    basic_cvector<double, std::complex<double> > vTauQ (m);
    basic_cvector<double, std::complex<double> > vTauP (m);
    basic_cmatrix<double, std::complex<double> > mA (mArg);
    std::complex<double> dWork;

    ZGEBRD (&nM, &nN, mA, mA._pld(), mD, vOffDiag, vTauQ, vTauP, &dWork, &lWork, &nOutInfo);
    lWork = static_cast<int>(dWork.real());
    basic_cvector<double, std::complex<double> > vWork (lWork);

    ZGEBRD (&nM, &nN, mA, mA._pld(), mD, vOffDiag, vTauQ, vTauP, vWork, &lWork, &nOutInfo);
    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

    static const int zero(0);
    static const int one(1);
    if (bSimple)
    {
        basic_rvector<double> vWork2 (m * 4);
        ZBDSQR (nM >= nN ? Chars::pU() : Chars::pL(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                1,
#endif
                &m,
                &zero, &zero, &zero,
                mD, vOffDiag,
                NULL, &one, NULL, &one, NULL, &one,
                vWork2, &nOutInfo);
    }
    else
    {
        const int lWork3 = m * _cvm_max<int>(M, 64);    // bug fix for cases when N>64*M or M>64*N
        basic_rvector<double> vWork2 (m * 4);
        basic_cvector<double, std::complex<double> > vWork3 (lWork3);

        basic_cmatrix<double, std::complex<double> > Q (mA);
        basic_cmatrix<double, std::complex<double> > P (mA);

        if (nM > nN) Q.resize(nM, nM);
        if (nM < nN) P.resize(nN, nN);

        ZUNGBR (Chars::pQ(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                1,
#endif
                &nM, &nM, &nN,
                Q, &nM,
                vTauQ,
                vWork3, &lWork3, &nOutInfo);
        if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

        ZUNGBR (Chars::pP(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                1,
#endif
                &nN, &nN, &nM,
                P, &M,
                vTauP,
                vWork3, &lWork3, &nOutInfo);
        if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

        ZBDSQR (nM >= nN ? Chars::pU() : Chars::pL(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                1,
#endif
                &m,
                &nN, &nM, 
                &zero,
                mD, vOffDiag,
                P, &M, Q, &nM, NULL, &one,
                vWork2, &nOutInfo);

        (*mU)  = basic_scmatrix<double, std::complex<double> >(Q.resize(nM, nM));
        (*mVH) = basic_scmatrix<double, std::complex<double> >(P.resize(nN, nN));
    }

    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
    if (nOutInfo > 0) throw cvmexception (CVM_CONVERGENCE_ERROR);

    __copy<double> (nSize, mD, mD.incr(), pD, nIncr);
}

template<>
CVM_API void
__svd<float, basic_srbmatrix<float>, basic_srmatrix<float> >
    (float* pD, int nSize, int nIncr,
    const basic_srbmatrix<float>& mArg,
    basic_srmatrix<float>* mU,
    basic_srmatrix<float>* mVH) throw (cvmexception)
{
    static const int zero(0);
    const int m = mArg.msize();
    if (m != nSize) throw cvmexception (CVM_SIZESMISMATCH);

    const bool bSimple  = (mU == NULL || mVH == NULL);
    int nOutInfo = 0;

    basic_rvector<float> mD (nSize);
    mD.assign (pD, nIncr);

    basic_rvector<float>   vOffDiag (_cvm_max<int>(1, m - 1));
    basic_srmatrix<float>  mQ       (bSimple ? 1 : m);
    basic_srmatrix<float>  mPT      (bSimple ? 1 : m);
    basic_srmatrix<float>  mC       (1);
    basic_rvector<float>   vWork    (2 * m);
    basic_srbmatrix<float> mA       (mArg);

    SGBBRD (bSimple ? Chars::pN() : Chars::pB(), 
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
            1,
#endif
            mA._pm(), mA._pn(), &zero, mA._pl(), mA._pu(), mA, mA._pld(),
            mD, vOffDiag,
            mQ, mQ._pm(),
            mPT, mPT._pm(),
            mC, mC._pm(), 
            vWork, &nOutInfo);
    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

    basic_rvector<float> vWork2 (m * 4);
    SBDSQR (Chars::pU(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
            1,
#endif
            &m,
            bSimple ? &zero : &m,
            bSimple ? &zero : &m,
            &zero,
            mD, vOffDiag,
            mPT, mPT._pm(),
            mQ, mQ._pm(),
            mC, mC._pm(), 
            vWork2, &nOutInfo);
    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
    if (nOutInfo > 0) throw cvmexception (CVM_CONVERGENCE_ERROR);

    if (!bSimple)
    {
        (*mU)  = mQ;
        (*mVH) = mPT;
    }

    __copy<float> (nSize, mD, mD.incr(), pD, nIncr);
}

template<>
CVM_API void
__svd<double, basic_srbmatrix<double>, basic_srmatrix<double> >
    (double* pD, int nSize, int nIncr,
    const basic_srbmatrix<double>& mArg,
    basic_srmatrix<double>* mU,
    basic_srmatrix<double>* mVH) throw (cvmexception)
{
    static const int zero(0);
    const int m = mArg.msize();
    if (m != nSize) throw cvmexception (CVM_SIZESMISMATCH);

    const bool bSimple  = (mU == NULL || mVH == NULL);
    int nOutInfo = 0;

    basic_rvector<double> mD (nSize);
    mD.assign (pD, nIncr);

    basic_rvector<double>   vOffDiag (_cvm_max<int>(1, m - 1));
    basic_srmatrix<double>  mQ       (bSimple ? 1 : m);
    basic_srmatrix<double>  mPT      (bSimple ? 1 : m);
    basic_srmatrix<double>  mC       (1);
    basic_rvector<double>   vWork    (2 * m);
    basic_srbmatrix<double> mA       (mArg);

    DGBBRD (bSimple ? Chars::pN() : Chars::pB(), 
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
            1,
#endif
            mA._pm(), mA._pn(), &zero, mA._pl(), mA._pu(), mA, mA._pld(),
            mD, vOffDiag,
            mQ, mQ._pm(),
            mPT, mPT._pm(),
            mC, mC._pm(), 
            vWork, &nOutInfo);
    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

    basic_rvector<double> vWork2 (m * 4);
    DBDSQR (Chars::pU(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
            1,
#endif
            &m,
            bSimple ? &zero : &m,
            bSimple ? &zero : &m,
            &zero,
            mD, vOffDiag,
            mPT, mPT._pm(),
            mQ, mQ._pm(),
            mC, mC._pm(), 
            vWork2, &nOutInfo);
    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
    if (nOutInfo > 0) throw cvmexception (CVM_CONVERGENCE_ERROR);

    if (!bSimple)
    {
        (*mU)  = mQ;
        (*mVH) = mPT;
    }

    __copy<double> (nSize, mD, mD.incr(), pD, nIncr);
}

template<>
CVM_API void
__svd<float, basic_scbmatrix<float, std::complex<float> >, basic_scmatrix<float, std::complex<float> > >
    (float* pD, int nSize, int nIncr,
    const basic_scbmatrix<float, std::complex<float> >& mArg,
    basic_scmatrix<float, std::complex<float> >* mU,
    basic_scmatrix<float, std::complex<float> >* mVH) throw (cvmexception)
{
    static const int zero(0);
    const int m = mArg.msize();
    if (m != nSize) throw cvmexception (CVM_SIZESMISMATCH);

    const bool bSimple  = (mU == NULL || mVH == NULL);
    int nOutInfo = 0;

    basic_rvector<float> mD (nSize);
    mD.assign (pD, nIncr);

    basic_rvector<float>   vOffDiag (_cvm_max<int>(1, m - 1));
    basic_scmatrix<float, std::complex<float> >  mQ       (bSimple ? 1 : m);
    basic_scmatrix<float, std::complex<float> >  mPT      (bSimple ? 1 : m);
    basic_scmatrix<float, std::complex<float> >  mC       (1);
    basic_cvector<float, std::complex<float> >   vWork    (m);
    basic_rvector<float>   vRWork   (m);
    basic_scbmatrix<float, std::complex<float> > mA       (mArg);

    CGBBRD (bSimple ? Chars::pN() : Chars::pB(), 
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
            1,
#endif
            mA._pm(), mA._pn(), &zero, mA._pl(), mA._pu(), mA, mA._pld(),
            mD, vOffDiag,
            mQ, mQ._pm(),
            mPT, mPT._pm(),
            mC, mC._pm(), 
            vWork, vRWork, &nOutInfo);
    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

    basic_rvector<float> vWork2 (m * 4);
    CBDSQR (Chars::pU(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
            1,
#endif
            &m,
            bSimple ? &zero : &m,
            bSimple ? &zero : &m,
            &zero,
            mD, vOffDiag,
            mPT, mPT._pm(),
            mQ, mQ._pm(),
            mC, mC._pm(), 
            vWork2, &nOutInfo);
    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
    if (nOutInfo > 0) throw cvmexception (CVM_CONVERGENCE_ERROR);

    if (!bSimple)
    {
        (*mU)  = mQ;
        (*mVH) = mPT;
    }

    __copy<float> (nSize, mD, mD.incr(), pD, nIncr);
}

template<>
CVM_API void
__svd<double, basic_scbmatrix<double, std::complex<double> >, basic_scmatrix<double, std::complex<double> > >
    (double* pD, int nSize, int nIncr,
    const basic_scbmatrix<double, std::complex<double> >& mArg,
    basic_scmatrix<double, std::complex<double> >* mU,
    basic_scmatrix<double, std::complex<double> >* mVH) throw (cvmexception)
{
    static const int zero(0);
    const int m = mArg.msize();
    if (m != nSize) throw cvmexception (CVM_SIZESMISMATCH);

    const bool bSimple  = (mU == NULL || mVH == NULL);
    int nOutInfo = 0;

    basic_rvector<double> mD (nSize);
    mD.assign (pD, nIncr);

    basic_rvector<double>   vOffDiag (_cvm_max<int>(1, m - 1));
    basic_scmatrix<double, std::complex<double> >  mQ       (bSimple ? 1 : m);
    basic_scmatrix<double, std::complex<double> >  mPT      (bSimple ? 1 : m);
    basic_scmatrix<double, std::complex<double> >  mC       (1);
    basic_cvector<double, std::complex<double> >   vWork    (m);
    basic_rvector<double>   vRWork   (m);
    basic_scbmatrix<double, std::complex<double> > mA       (mArg);

    ZGBBRD (bSimple ? Chars::pN() : Chars::pB(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
            1,
#endif
            mA._pm(), mA._pn(), &zero, mA._pl(), mA._pu(), mA, mA._pld(),
            mD, vOffDiag,
            mQ, mQ._pm(),
            mPT, mPT._pm(),
            mC, mC._pm(),
            vWork, vRWork, &nOutInfo);
    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

    basic_rvector<double> vWork2 (m * 4);
    ZBDSQR (Chars::pU(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
            1,
#endif
            &m,
            bSimple ? &zero : &m,
            bSimple ? &zero : &m,
            &zero,
            mD, vOffDiag,
            mPT, mPT._pm(),
            mQ, mQ._pm(),
            mC, mC._pm(),
            vWork2, &nOutInfo);
    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
    if (nOutInfo > 0) throw cvmexception (CVM_CONVERGENCE_ERROR);

    if (!bSimple)
    {
        (*mU)  = mQ;
        (*mVH) = mPT;
    }

    __copy<double> (nSize, mD, mD.incr(), pD, nIncr);
}


static int _ssyevd_lwork (int n)
{
    int k = 1, nn = n - 1;
    while (nn >>= 1) ++k;
    return 3 * n * n + (5 + 2 * k) * n + 1;
}

template<>
CVM_API void
__eig<basic_rvector<float>, basic_srsmatrix<float>, basic_srmatrix<float> >
    (basic_rvector<float>& vRes,
    const basic_srsmatrix<float>& mArg,
    basic_srmatrix<float>* mEigVect,
    bool /*bRightVect*/) throw (cvmexception)
{
    const int nM = mArg.msize();
    if (nM != vRes.size()) throw cvmexception (CVM_SIZESMISMATCH);
    const bool bEigVect = (mEigVect != NULL);

    if (nM == 1)
    {
        vRes[1] = mArg(1,1);
        if (bEigVect)
        {
            static const float one(1.F);
            mEigVect -> resize (1);
            (*mEigVect)[1].set (one);
        }
    }
    else
    {
        const char* pcJob = bEigVect ? Chars::pV() : Chars::pN();
        basic_srsmatrix<float> mA (mArg);
        int lwork = bEigVect ? _ssyevd_lwork (nM) : (2 * nM + 1);               // LAPACK wants (1 + 6*N + 2*N**2) though
        int liwork = bEigVect ? (5 * nM + 3) : 1;                               // MKL wants 5*n + 2 though
        basic_rvector<float> work (lwork);
        basic_array<int> iwork (liwork);
        int nOutInfo = 0;

        SSYEVD (pcJob,
    #if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                1,
    #endif
                Chars::pU(),
    #if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                1,
    #endif
                &nM, mA, mA._pld(), vRes, work, &lwork, iwork, &liwork, &nOutInfo);
        if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
        if (nOutInfo > 0) throw cvmexception (CVM_CONVERGENCE_ERROR);

        if (bEigVect)
        {
            (*mEigVect) << mA;
        }
    }
}

template<>
CVM_API void
__eig<basic_rvector<double>, basic_srsmatrix<double>, basic_srmatrix<double> >
    (basic_rvector<double>& vRes,
    const basic_srsmatrix<double>& mArg,
    basic_srmatrix<double>* mEigVect,
    bool /*bRightVect*/) throw (cvmexception)
{
    const int nM = mArg.msize();
    if (nM != vRes.size()) throw cvmexception (CVM_SIZESMISMATCH);
    const bool bEigVect = (mEigVect != NULL);

    if (nM == 1)
    {
        vRes[1] = mArg(1,1);
        if (bEigVect)
        {
            static const double one(1.);
            mEigVect -> resize (1);
            (*mEigVect)[1].set (one);
        }
    }
    else
    {
        const char* pcJob = bEigVect ? Chars::pV() : Chars::pN();
        basic_srsmatrix<double> mA (mArg);
        int lwork = bEigVect ? _ssyevd_lwork (nM) : (2 * nM + 1);               // LAPACK wants (1 + 6*N + 2*N**2) though
        int liwork = bEigVect ? (5 * nM + 3) : 1;                               // MKL wants 5*n + 2 though
        basic_rvector<double> work (lwork);
        basic_array<int> iwork (liwork);
        int nOutInfo = 0;

        DSYEVD (pcJob,
    #if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                1,
    #endif
                Chars::pU(),
    #if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                1,
    #endif
                &nM, mA, mA._pld(), vRes, work, &lwork, iwork, &liwork, &nOutInfo);
        if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
        if (nOutInfo > 0) throw cvmexception (CVM_CONVERGENCE_ERROR);

        if (bEigVect)
        {
            (*mEigVect) << mA;
        }
    }
}

static int _cheevd_lwork (int n)
{
    int k = 1, nn = n - 1;
    while (nn >>= 1) ++k;
    return 3 * n * n + (4 + 2 * k) * n + 1;
}

template<>
CVM_API void
__eig<basic_rvector<float>, basic_schmatrix<float, std::complex<float> >, basic_scmatrix<float, std::complex<float> > >
    (basic_rvector<float>& vRes,
    const basic_schmatrix<float, std::complex<float> >& mArg,
    basic_scmatrix<float, std::complex<float> >* mEigVect,
    bool /*bRightVect*/) throw (cvmexception)
{
    const int nM = mArg.msize();
    if (nM != vRes.size()) throw cvmexception (CVM_SIZESMISMATCH);
    const bool bEigVect = (mEigVect != NULL);

    if (nM == 1)
    {
        vRes[1] = 1.F;
        if (bEigVect)
        {
            mEigVect -> resize (1);
            (*mEigVect)[1].set (mArg(1,1));
        }
    }
    else
    {
        const char* pcJob = bEigVect ? Chars::pV() : Chars::pN();
        basic_schmatrix<float, std::complex<float> > mA (mArg);
        int lwork  = bEigVect ? nM * (nM + 2) : nM + 1;
        int lrwork = bEigVect ? _cheevd_lwork (nM) : (2 * nM + 1);  // LAPACK wants (1 + 5*N + 2*N**2)  (or N for 'N' mode) though
        int liwork = bEigVect ? (5 * nM + 3) : 1;                   // MKL wants 5*n + 2 though
        basic_cvector<float, std::complex<float> > work (lwork);
        basic_rvector<float> rwork (lrwork);
        basic_array<int> iwork (liwork);
        int nOutInfo = 0;

        CHEEVD (pcJob,
    #if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                1,
    #endif
                Chars::pU(),
    #if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                1,
    #endif
                &nM, mA, mA._pld(), vRes, work, &lwork, rwork, &lrwork, iwork, &liwork, &nOutInfo);
        if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
        if (nOutInfo > 0) throw cvmexception (CVM_CONVERGENCE_ERROR);

        if (bEigVect)
        {
            (*mEigVect) << mA;
        }
    }
}

template<>
CVM_API void
__eig<basic_rvector<double>, basic_schmatrix<double, std::complex<double> >, basic_scmatrix<double, std::complex<double> > >
    (basic_rvector<double>& vRes,
    const basic_schmatrix<double, std::complex<double> >& mArg,
    basic_scmatrix<double, std::complex<double> >* mEigVect,
    bool /*bRightVect*/) throw (cvmexception)
{
    const int nM = mArg.msize();
    if (nM != vRes.size()) throw cvmexception (CVM_SIZESMISMATCH);
    const bool bEigVect = (mEigVect != NULL);

    if (nM == 1)
    {
        vRes[1] = 1.;
        if (bEigVect)
        {
            mEigVect -> resize (1);
            (*mEigVect)[1].set (mArg(1,1));
        }
    }
    else
    {
        const char* pcJob = bEigVect ? Chars::pV() : Chars::pN();
        basic_schmatrix<double, std::complex<double> > mA (mArg);
        int lwork  = bEigVect ? nM * (nM + 2) : nM + 1;
        int lrwork = bEigVect ? _cheevd_lwork (nM) : (2 * nM + 1);  // LAPACK wants (1 + 5*N + 2*N**2)  (or N for 'N' mode) though
        int liwork = bEigVect ? (5 * nM + 3) : 1;                   // MKL wants 5*n + 2 though
        basic_cvector<double, std::complex<double> > work (lwork);
        basic_rvector<double> rwork (lrwork);
        basic_array<int> iwork (liwork);
        int nOutInfo = 0;

        ZHEEVD (pcJob,
    #if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                1,
    #endif
                Chars::pU(),
    #if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                1,
    #endif
                &nM, mA, mA._pld(), vRes, work, &lwork, rwork, &lrwork, iwork, &liwork, &nOutInfo);
        if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
        if (nOutInfo > 0) throw cvmexception (CVM_CONVERGENCE_ERROR);

        if (bEigVect)
        {
            (*mEigVect) << mA;
        }
    }
}

// pseudo (generalized) inversion routines
template<>
CVM_API void
__pinv<float, basic_rmatrix<float>, basic_rmatrix<float> >
    (basic_rmatrix<float>& mX, 
    const basic_rmatrix<float>& mArg, float threshold) throw (cvmexception)
{
    const int nM = mArg.msize();
    const int nN = mArg.nsize();
    const int m  = _cvm_min<int>(nM, nN);
    const int M  = _cvm_max<int>(nM, nN);
    int lWork    = -1; // to calculate lWork
    int nOutInfo = 0;

    basic_rvector<float> mD       (m);
    basic_rvector<float> vOffDiag (_cvm_max<int>(1, m - 1));
    basic_rvector<float> vTauQ    (m);
    basic_rvector<float> vTauP    (m);
    basic_rmatrix<float> mA       (mArg);
    float dWork;

    SGEBRD (&nM, &nN, mA, mA._pld(), mD, vOffDiag, vTauQ, vTauP, &dWork, &lWork, &nOutInfo);
    lWork = static_cast<int>(dWork);
    basic_rvector<float> vWork(lWork);

    SGEBRD (&nM, &nN, mA, mA._pld(), mD, vOffDiag, vTauQ, vTauP, vWork, &lWork, &nOutInfo);
    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

    static const int zero(0);
    static const int one(1);
    const int lWork3 = m * _cvm_max<int>(M, 64);    // bug fix for cases when N>64*M or M>64*N
    basic_rvector<float> vWork2 (m * 4);
    basic_rvector<float> vWork3 (lWork3);
    
    // few words about economy:
    // for m > n case we care about m-by-n matrix U and n-by-n matrix VH
    // for m < n case we care about m-by-m matrix U and m-by-n matrix VH
    // however, the whole matrix A is needed to start computations
    basic_rmatrix<float> Q (mA);
    basic_rmatrix<float> P (mA);

    SORGBR (Chars::pQ(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
            1,
#endif
            &nM, &m, &nN,
            Q, Q._pld(),
            vTauQ,
            vWork3, &lWork3, &nOutInfo);
    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

    SORGBR (Chars::pP(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
            1,
#endif
            &m, &nN, &nM,
            P, P._pld(),
            vTauP,
            vWork3, &lWork3, &nOutInfo);
    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

    SBDSQR (nM >= nN ? Chars::pU() : Chars::pL(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
            1,
#endif
            &m,
            &nN, &nM,
            &zero,
            mD, vOffDiag,
            P, P._pld(), Q, Q._pld(), NULL, &one,
            vWork2, &nOutInfo);
    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
    if (nOutInfo > 0) throw cvmexception (CVM_CONVERGENCE_ERROR);

    if (nM > nN) P.resize(nN, nN);   // VH
    if (nM < nN) Q.resize(nM, nM);   // U
    for (int i = 1; i <= P.msize(); ++i) {
        if (mD[i] > threshold) {
            P[i] /= mD[i];  // multiplying V by S^-1
        }
        else {
            P[i].set(0.);
        }
    }
    mX.mult (~P, ~Q);
}

template<>
CVM_API void
__pinv<double, basic_rmatrix<double>, basic_rmatrix<double> >
    (basic_rmatrix<double>& mX, 
    const basic_rmatrix<double>& mArg, double threshold) throw (cvmexception)
{
    const int nM = mArg.msize();
    const int nN = mArg.nsize();
    const int m  = _cvm_min<int>(nM, nN);
    const int M  = _cvm_max<int>(nM, nN);
    int lWork    = -1; // to calculate lWork
    int nOutInfo = 0;

    basic_rvector<double> mD       (m);
    basic_rvector<double> vOffDiag (_cvm_max<int>(1, m - 1));
    basic_rvector<double> vTauQ    (m);
    basic_rvector<double> vTauP    (m);
    basic_rmatrix<double> mA       (mArg);
    double dWork;

    DGEBRD (&nM, &nN, mA, mA._pld(), mD, vOffDiag, vTauQ, vTauP, &dWork, &lWork, &nOutInfo);
    lWork = static_cast<int>(dWork);
    basic_rvector<double> vWork(lWork);

    DGEBRD (&nM, &nN, mA, mA._pld(), mD, vOffDiag, vTauQ, vTauP, vWork, &lWork, &nOutInfo);
    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

    static const int zero(0);
    static const int one(1);
    const int lWork3 = m * _cvm_max<int>(M, 64);    // bug fix for cases when N>64*M or M>64*N
    basic_rvector<double> vWork2 (m * 4);
    basic_rvector<double> vWork3 (lWork3);
    
    // few words about economy:
    // for m > n case we care about m-by-n matrix U and n-by-n matrix VH
    // for m < n case we care about m-by-m matrix U and m-by-n matrix VH
    // however, the whole matrix A is needed to start computations
    basic_rmatrix<double> Q (mA);
    basic_rmatrix<double> P (mA);

    DORGBR (Chars::pQ(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
            1,
#endif
            &nM, &m, &nN,
            Q, Q._pld(),
            vTauQ,
            vWork3, &lWork3, &nOutInfo);
    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

    DORGBR (Chars::pP(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
            1,
#endif
            &m, &nN, &nM,
            P, P._pld(),
            vTauP,
            vWork3, &lWork3, &nOutInfo);
    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

    DBDSQR (nM >= nN ? Chars::pU() : Chars::pL(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
            1,
#endif
            &m,
            &nN, &nM,
            &zero,
            mD, vOffDiag,
            P, P._pld(), Q, Q._pld(), NULL, &one,
            vWork2, &nOutInfo);
    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
    if (nOutInfo > 0) throw cvmexception (CVM_CONVERGENCE_ERROR);

    if (nM > nN) P.resize(nN, nN);   // VH
    if (nM < nN) Q.resize(nM, nM);   // U
    for (int i = 1; i <= P.msize(); ++i) {
        if (mD[i] > threshold) {
            P[i] /= mD[i];  // multiplying V by S^-1
        }
        else {
            P[i].set(0.);
        }
    }
    mX.mult (~P, ~Q);
}

template<>
CVM_API void
__pinv<float, basic_cmatrix<float, std::complex<float> >, basic_cmatrix<float, std::complex<float> > >
    (basic_cmatrix<float, std::complex<float> >& mX, 
    const basic_cmatrix<float, std::complex<float> >& mArg, float threshold) throw (cvmexception)
{
    const int nM = mArg.msize();
    const int nN = mArg.nsize();
    const int m  = _cvm_min<int>(nM, nN);
    const int M  = _cvm_max<int>(nM, nN);
    int lWork    = -1; // to calculate lWork
    int nOutInfo = 0;

    basic_rvector<float> mD       (m);
    basic_rvector<float> vOffDiag (_cvm_max<int>(1, m - 1));
    basic_cvector<float, std::complex<float> > vTauQ    (m);
    basic_cvector<float, std::complex<float> > vTauP    (m);
    basic_cmatrix<float, std::complex<float> > mA       (mArg);
    std::complex<float>  dWork;

    CGEBRD (&nM, &nN, mA, mA._pld(), mD, vOffDiag, vTauQ, vTauP, &dWork, &lWork, &nOutInfo);
    lWork = static_cast<int>(dWork.real());
    basic_cvector<float, std::complex<float> > vWork(lWork);

    CGEBRD (&nM, &nN, mA, mA._pld(), mD, vOffDiag, vTauQ, vTauP, vWork, &lWork, &nOutInfo);
    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

    static const int zero(0);
    static const int one(1);
    const int lWork3 = m * _cvm_max<int>(M, 64);    // bug fix for cases when N>64*M or M>64*N
    basic_rvector<float> vWork2 (m * 4);
    basic_cvector<float, std::complex<float> > vWork3 (lWork3);
    
    // few words about economy:
    // for m > n case we care about m-by-n matrix U and n-by-n matrix VH
    // for m < n case we care about m-by-m matrix U and m-by-n matrix VH
    // however, the whole matrix A is needed to start computations
    basic_cmatrix<float, std::complex<float> > Q (mA);
    basic_cmatrix<float, std::complex<float> > P (mA);

    CUNGBR (Chars::pQ(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
            1,
#endif
            &nM, &m, &nN,
            Q, Q._pld(),
            vTauQ,
            vWork3, &lWork3, &nOutInfo);
    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

    CUNGBR (Chars::pP(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
            1,
#endif
            &m, &nN, &nM,
            P, P._pld(),
            vTauP,
            vWork3, &lWork3, &nOutInfo);
    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

    CBDSQR (nM >= nN ? Chars::pU() : Chars::pL(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
            1,
#endif
            &m,
            &nN, &nM,
            &zero,
            mD, vOffDiag,
            P, P._pld(), Q, Q._pld(), NULL, &one,
            vWork2, &nOutInfo);
    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
    if (nOutInfo > 0) throw cvmexception (CVM_CONVERGENCE_ERROR);

    if (nM > nN) P.resize(nN, nN);   // VH
    if (nM < nN) Q.resize(nM, nM);   // U
    for (int i = 1; i <= P.msize(); ++i) {
        if (mD[i] > threshold) {
            P[i] /= mD[i];  // multiplying V by S^-1
        }
        else {
            P[i].set(0.);
        }
    }
    mX.mult (~P, ~Q);
}

template<>
CVM_API void
__pinv<double, basic_cmatrix<double, std::complex<double> >, basic_cmatrix<double, std::complex<double> > >
    (basic_cmatrix<double, std::complex<double> >& mX, 
    const basic_cmatrix<double, std::complex<double> >& mArg, double threshold) throw (cvmexception)
{
    const int nM = mArg.msize();
    const int nN = mArg.nsize();
    const int m  = _cvm_min<int>(nM, nN);
    const int M  = _cvm_max<int>(nM, nN);
    int lWork    = -1; // to calculate lWork
    int nOutInfo = 0;

    basic_rvector<double> mD       (m);
    basic_rvector<double> vOffDiag (_cvm_max<int>(1, m - 1));
    basic_cvector<double, std::complex<double> > vTauQ    (m);
    basic_cvector<double, std::complex<double> > vTauP    (m);
    basic_cmatrix<double, std::complex<double> > mA       (mArg);
    std::complex<double>  dWork;

    ZGEBRD (&nM, &nN, mA, mA._pld(), mD, vOffDiag, vTauQ, vTauP, &dWork, &lWork, &nOutInfo);
    lWork = static_cast<int>(dWork.real());
    basic_cvector<double, std::complex<double> > vWork(lWork);

    ZGEBRD (&nM, &nN, mA, mA._pld(), mD, vOffDiag, vTauQ, vTauP, vWork, &lWork, &nOutInfo);
    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

    static const int zero(0);
    static const int one(1);
    const int lWork3 = m * _cvm_max<int>(M, 64);    // bug fix for cases when N>64*M or M>64*N
    basic_rvector<double> vWork2 (m * 4);
    basic_cvector<double, std::complex<double> > vWork3 (lWork3);
    
    // few words about economy:
    // for m > n case we care about m-by-n matrix U and n-by-n matrix VH
    // for m < n case we care about m-by-m matrix U and m-by-n matrix VH
    // however, the whole matrix A is needed to start computations
    basic_cmatrix<double, std::complex<double> > Q (mA);
    basic_cmatrix<double, std::complex<double> > P (mA);

    ZUNGBR (Chars::pQ(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
            1,
#endif
            &nM, &m, &nN,
            Q, Q._pld(),
            vTauQ,
            vWork3, &lWork3, &nOutInfo);
    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

    ZUNGBR (Chars::pP(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
            1,
#endif
            &m, &nN, &nM,
            P, P._pld(),
            vTauP,
            vWork3, &lWork3, &nOutInfo);
    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

    ZBDSQR (nM >= nN ? Chars::pU() : Chars::pL(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
            1,
#endif
            &m,
            &nN, &nM,
            &zero,
            mD, vOffDiag,
            P, P._pld(), Q, Q._pld(), NULL, &one,
            vWork2, &nOutInfo);
    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
    if (nOutInfo > 0) throw cvmexception (CVM_CONVERGENCE_ERROR);

    if (nM > nN) P.resize(nN, nN);   // VH
    if (nM < nN) Q.resize(nM, nM);   // U
    for (int i = 1; i <= P.msize(); ++i) {
        if (mD[i] > threshold) {
            P[i] /= mD[i];  // multiplying V by S^-1
        }
        else {
            P[i].set(0.);
        }
    }
    mX.mult (~P, ~Q);
}

template<>
CVM_API void
__pinv<float, basic_srbmatrix<float>, basic_rmatrix<float> >
    (basic_rmatrix<float>& mX, 
    const basic_srbmatrix<float>& mArg, float threshold) throw (cvmexception)
{
    static const int zero(0);
    const int m = mArg.msize();
    int nOutInfo = 0;

    basic_rvector<float>   mD       (m);
    basic_rvector<float>   vOffDiag (_cvm_max<int>(1, m - 1));
    basic_srmatrix<float>  mQ       (m);
    basic_srmatrix<float>  mPT      (m);
    basic_srmatrix<float>  mC       (1);
    basic_rvector<float>   vWork    (2 * m);
    basic_srbmatrix<float> mA       (mArg);

    SGBBRD (Chars::pB(), 
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
            1,
#endif
            mA._pm(), mA._pn(), &zero, mA._pl(), mA._pu(), mA, mA._pld(),
            mD, vOffDiag,
            mQ, mQ._pm(),
            mPT, mPT._pm(),
            mC, mC._pm(), 
            vWork, &nOutInfo);
    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

    basic_rvector<float> vWork2 (m * 4);
    SBDSQR (Chars::pU(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
            1,
#endif
            &m,
            &m,
            &m,
            &zero,
            mD, vOffDiag,
            mPT, mPT._pm(),
            mQ, mQ._pm(),
            mC, mC._pm(), 
            vWork2, &nOutInfo);
    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
    if (nOutInfo > 0) throw cvmexception (CVM_CONVERGENCE_ERROR);

    for (int i = 1; i <= m; ++i) {
        if (mD[i] > threshold) {
            mPT[i] /= mD[i];  // multiplying V by S^-1
        }
        else {
            mPT[i].set(0.F);
        }
    }
    mX.mult (~mPT, ~mQ);
}

template<>
CVM_API void
__pinv<double, basic_srbmatrix<double>, basic_rmatrix<double> >
    (basic_rmatrix<double>& mX, 
    const basic_srbmatrix<double>& mArg, double threshold) throw (cvmexception)
{
    static const int zero(0);
    const int m = mArg.msize();
    int nOutInfo = 0;

    basic_rvector<double>   mD       (m);
    basic_rvector<double>   vOffDiag (_cvm_max<int>(1, m - 1));
    basic_srmatrix<double>  mQ       (m);
    basic_srmatrix<double>  mPT      (m);
    basic_srmatrix<double>  mC       (1);
    basic_rvector<double>   vWork    (2 * m);
    basic_srbmatrix<double> mA       (mArg);

    DGBBRD (Chars::pB(), 
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
            1,
#endif
            mA._pm(), mA._pn(), &zero, mA._pl(), mA._pu(), mA, mA._pld(),
            mD, vOffDiag,
            mQ, mQ._pm(),
            mPT, mPT._pm(),
            mC, mC._pm(), 
            vWork, &nOutInfo);
    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

    basic_rvector<double> vWork2 (m * 4);
    DBDSQR (Chars::pU(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
            1,
#endif
            &m,
            &m,
            &m,
            &zero,
            mD, vOffDiag,
            mPT, mPT._pm(),
            mQ, mQ._pm(),
            mC, mC._pm(), 
            vWork2, &nOutInfo);
    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
    if (nOutInfo > 0) throw cvmexception (CVM_CONVERGENCE_ERROR);

    for (int i = 1; i <= m; ++i) {
        if (mD[i] > threshold) {
            mPT[i] /= mD[i];  // multiplying V by S^-1
        }
        else {
            mPT[i].set(0.);
        }
    }
    mX.mult (~mPT, ~mQ);
}


template<>
CVM_API void
__pinv<float, basic_scbmatrix<float, std::complex<float> >, basic_cmatrix<float, std::complex<float> > >
    (basic_cmatrix<float, std::complex<float> >& mX, 
    const basic_scbmatrix<float, std::complex<float> >& mArg, float threshold) throw (cvmexception)
{
    static const int zero(0);
    const int m = mArg.msize();
    int nOutInfo = 0;

    basic_rvector<float> mD (m);
    basic_rvector<float> vOffDiag (_cvm_max<int>(1, m - 1));
    basic_scmatrix<float, std::complex<float> >  mQ    (m);
    basic_scmatrix<float, std::complex<float> >  mPT   (m);
    basic_scmatrix<float, std::complex<float> >  mC    (1);
    basic_cvector<float, std::complex<float> >   vWork (2 * m);
    basic_rvector<float> vRWork (m);
    basic_scbmatrix<float, std::complex<float> > mA    (mArg);

    CGBBRD (Chars::pB(), 
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
            1,
#endif
            mA._pm(), mA._pn(), &zero, mA._pl(), mA._pu(), mA, mA._pld(),
            mD, vOffDiag,
            mQ, mQ._pm(),
            mPT, mPT._pm(),
            mC, mC._pm(), 
            vWork, vRWork, &nOutInfo);
    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

    basic_rvector<float> vWork2 (m * 4);
    CBDSQR (Chars::pU(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
            1,
#endif
            &m,
            &m,
            &m,
            &zero,
            mD, vOffDiag,
            mPT, mPT._pm(),
            mQ, mQ._pm(),
            mC, mC._pm(), 
            vWork2, &nOutInfo);
    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
    if (nOutInfo > 0) throw cvmexception (CVM_CONVERGENCE_ERROR);

    for (int i = 1; i <= m; ++i) {
        if (mD[i] > threshold) {
            mPT[i] /= mD[i];  // multiplying V by S^-1
        }
        else {
            mPT[i].set(0.F);
        }
    }
    mX.mult (~mPT, ~mQ);
}

template<>
CVM_API void
__pinv<double, basic_scbmatrix<double, std::complex<double> >, basic_cmatrix<double, std::complex<double> > >
    (basic_cmatrix<double, std::complex<double> >& mX, 
    const basic_scbmatrix<double, std::complex<double> >& mArg, double threshold) throw (cvmexception)
{
    static const int zero(0);
    const int m = mArg.msize();
    int nOutInfo = 0;

    basic_rvector<double> mD (m);
    basic_rvector<double> vOffDiag (_cvm_max<int>(1, m - 1));
    basic_scmatrix<double, std::complex<double> >  mQ    (m);
    basic_scmatrix<double, std::complex<double> >  mPT   (m);
    basic_scmatrix<double, std::complex<double> >  mC    (1);
    basic_cvector<double, std::complex<double> >   vWork (2 * m);
    basic_rvector<double> vRWork (m);
    basic_scbmatrix<double, std::complex<double> > mA (mArg);

    ZGBBRD (Chars::pB(), 
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
            1,
#endif
            mA._pm(), mA._pn(), &zero, mA._pl(), mA._pu(), mA, mA._pld(),
            mD, vOffDiag,
            mQ, mQ._pm(),
            mPT, mPT._pm(),
            mC, mC._pm(), 
            vWork, vRWork, &nOutInfo);
    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

    basic_rvector<double> vWork2 (m * 4);
    ZBDSQR (Chars::pU(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
            1,
#endif
            &m,
            &m,
            &m,
            &zero,
            mD, vOffDiag,
            mPT, mPT._pm(),
            mQ, mQ._pm(),
            mC, mC._pm(), 
            vWork2, &nOutInfo);
    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
    if (nOutInfo > 0) throw cvmexception (CVM_CONVERGENCE_ERROR);

    for (int i = 1; i <= m; ++i) {
        if (mD[i] > threshold) {
            mPT[i] /= mD[i];  // multiplying V by S^-1
        }
        else {
            mPT[i].set(0.);
        }
    }
    mX.mult (~mPT, ~mQ);
}

CVM_NAMESPACE_END
