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
__gemm<std::complex<float>, basic_cmatrix<float, std::complex<float> > >
    (const basic_cmatrix<float, std::complex<float> >& ml, bool bTrans1,
     const basic_cmatrix<float, std::complex<float> >& mr, bool bTrans2,
     std::complex<float> dAlpha, 
     basic_cmatrix<float, std::complex<float> >& mRes, 
     std::complex<float> dBeta)
{
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
    CGEMM (bTrans1 ? Chars::pC() : Chars::pN(), 1, bTrans2 ? Chars::pC() : Chars::pN(), 1,
#else
    CGEMM (bTrans1 ? Chars::pC() : Chars::pN(),    bTrans2 ? Chars::pC() : Chars::pN(),
#endif
           bTrans1 ? ml._pn() : ml._pm(),
           bTrans2 ? mr._pm() : mr._pn(),
           bTrans1 ? ml._pm() : ml._pn(),
           &dAlpha,
           ml._pd(), ml._pldm(),
           mr._pd(), mr._pldm(),
           &dBeta,
           mRes, mRes._pld());
}

template<>
CVM_API void
__gemm<std::complex<double>, basic_cmatrix<double, std::complex<double> > >
    (const basic_cmatrix<double, std::complex<double> >& ml, bool bTrans1,
     const basic_cmatrix<double, std::complex<double> >& mr, bool bTrans2,
     std::complex<double> dAlpha,
     basic_cmatrix<double, std::complex<double> >& mRes,
     std::complex<double> dBeta)
{
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
    ZGEMM (bTrans1 ? Chars::pC() : Chars::pN(), 1, bTrans2 ? Chars::pC() : Chars::pN(), 1,
#else
    ZGEMM (bTrans1 ? Chars::pC() : Chars::pN(),    bTrans2 ? Chars::pC() : Chars::pN(),
#endif
           bTrans1 ? ml._pn() : ml._pm(),
           bTrans2 ? mr._pm() : mr._pn(),
           bTrans1 ? ml._pm() : ml._pn(),
           &dAlpha,
           ml._pd(), ml._pldm(),
           mr._pd(), mr._pldm(),
           &dBeta,
           mRes, mRes._pld());
}

template<>
CVM_API void
__hemm<std::complex<float>, basic_schmatrix<float, std::complex<float> >, basic_cmatrix<float, std::complex<float> > >
    (bool bLeft,
     const basic_schmatrix<float, std::complex<float> >& ml,
     const basic_cmatrix<float, std::complex<float> >& mr, 
     std::complex<float> dAlpha,
     basic_cmatrix<float, std::complex<float> >& mRes,
     std::complex<float> dBeta)
{
    CHEMM (bLeft ? Chars::pL() : Chars::pR(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
           1,
#endif
           Chars::pU(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
           1,
#endif
           mRes._pm(), mRes._pn(),
           &dAlpha,
           ml._pd(), ml._pld(),
           mr._pd(), mr._pld(),
           &dBeta,
           mRes, mRes._pld());
}

template<>
CVM_API void
__hemm<std::complex<double>, basic_schmatrix<double, std::complex<double> >, basic_cmatrix<double, std::complex<double> > >
    (bool bLeft,
     const basic_schmatrix<double, std::complex<double> >& ml,
     const basic_cmatrix<double, std::complex<double> >& mr, 
     std::complex<double> dAlpha,
     basic_cmatrix<double, std::complex<double> >& mRes,
     std::complex<double> dBeta)
{
    ZHEMM (bLeft ? Chars::pL() : Chars::pR(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
           1,
#endif
           Chars::pU(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
           1,
#endif
           mRes._pm(), mRes._pn(),
           &dAlpha,
           ml._pd(), ml._pld(),
           mr._pd(), mr._pld(),
           &dBeta,
           mRes, mRes._pld());
}

template <>
CVM_API void
__herk<std::complex<float>, basic_schmatrix<float, std::complex<float> > >
    (bool bTransp, 
    std::complex<float> alpha, int k,
    const std::complex<float>* pA, int ldA,
    std::complex<float> beta, basic_schmatrix<float, std::complex<float> >& m)
{
    CHERK (Chars::pU(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
           1,
#endif
           bTransp ? Chars::pC() : Chars::pN(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
           1,
#endif
           m._pm(), &k, &alpha, pA, &ldA, &beta, m, m._pld());
}

template <>
CVM_API void
__herk<std::complex<double>, basic_schmatrix<double, std::complex<double> > >
    (bool bTransp, 
    std::complex<double> alpha, int k,
    const std::complex<double>* pA, int ldA,
    std::complex<double> beta, basic_schmatrix<double, std::complex<double> >& m)
{
    ZHERK (Chars::pU(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
           1,
#endif
           bTransp ? Chars::pC() : Chars::pN(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
           1,
#endif
           m._pm(), &k, &alpha, pA, &ldA, &beta, m, m._pld());
}

template <>
CVM_API void
__her2k<std::complex<float>, basic_schmatrix<float, std::complex<float> > >
    (bool bTransp, 
    std::complex<float> alpha, int k,
    const std::complex<float>* pA, int ldA,
    const std::complex<float>* pB, int ldB,
    std::complex<float> beta, basic_schmatrix<float, std::complex<float> >& m)
{
    CHER2K (Chars::pU(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
           1,
#endif
           bTransp ? Chars::pC() : Chars::pN(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
           1,
#endif
           m._pm(), &k, &alpha, pA, &ldA, pB, &ldB, &beta, m, m._pld());
}

template <>
CVM_API void
__her2k<std::complex<double>, basic_schmatrix<double, std::complex<double> > >
    (bool bTransp, 
    std::complex<double> alpha, int k,
    const std::complex<double>* pA, int ldA,
    const std::complex<double>* pB, int ldB,
    std::complex<double> beta, basic_schmatrix<double, std::complex<double> >& m)
{
    ZHER2K (Chars::pU(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
           1,
#endif
           bTransp ? Chars::pC() : Chars::pN(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
           1,
#endif
           m._pm(), &k, &alpha, pA, &ldA, pB, &ldB, &beta, m, m._pld());
}



// Case 1:
// call ZGEQRF to get R and TAU
// call ZORGQR to get Q (using TAU)

// Case 1: economy mode, A is (m x n) and Q is (m x n) and R is (n x n)
template <>
CVM_API void 
__qre<basic_cmatrix<float, std::complex<float> >, basic_scmatrix<float, std::complex<float> > >
    (const basic_cmatrix<float, std::complex<float> >& mArg, 
    basic_cmatrix<float, std::complex<float> >& mQ,
    basic_scmatrix<float, std::complex<float> >& mR) throw (cvmexception)
{
    const int nM = mArg.msize();
    const int nN = mArg.nsize();
    const int nK = _cvm_min<int>(nM, nN);

    // we will eventually overwrite mQ to be the output matrix
    mQ = mArg;
    basic_cvector<float, std::complex<float> > vTau (nK);

    int lWork = -1;
    int nOutInfo = 0;
    std::complex<float> dWork;

    // calculate size of workspace
    CGEQRF (&nM, &nN, mQ._pd(), mQ._pld(), vTau, &dWork, &lWork, &nOutInfo);
    lWork = static_cast<int>(dWork.real());
    basic_cvector<float, std::complex<float> > vWork(lWork);

    // now do most of hardwork, find R and TAU and Householderish bits of Q
    CGEQRF (&nM, &nN, mQ._pd(), mQ._pld(), vTau, vWork, &lWork, &nOutInfo);
    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

    // get upper-triangular R from overwritten A
    mR.vanish();
    for (int row = 1; row <= nK; ++row)
        for (int col = row; col <= nN; ++col)
            mR(row,col) = mQ(row,col);

    // calculate size of workspace for finding Q
    lWork = -1;
    CUNGQR (&nM, &nK, &nK, mQ._pd(), mQ._pld(), vTau, &dWork, &lWork, &nOutInfo);
    lWork = static_cast<int>(dWork.real());
    if (lWork > vWork.size()) vWork.resize(lWork);

    // find Q
    CUNGQR (&nM, &nK, &nK, mQ._pd(), mQ._pld(), vTau, vWork, &lWork, &nOutInfo);
    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
}

template <>
CVM_API void 
__qre<basic_cmatrix<double, std::complex<double> >, basic_scmatrix<double, std::complex<double> > >
    (const basic_cmatrix<double, std::complex<double> >& mArg, 
    basic_cmatrix<double, std::complex<double> >& mQ,
    basic_scmatrix<double, std::complex<double> >& mR) throw (cvmexception)
{
    const int nM = mArg.msize();
    const int nN = mArg.nsize();
    const int nK = _cvm_min<int>(nM, nN);

    // we will eventually overwrite mQ to be the output matrix
    mQ = mArg;
    basic_cvector<double, std::complex<double> > vTau (nK);

    int lWork = -1;
    int nOutInfo = 0;
    std::complex<double> dWork;

    // calculate size of workspace
    ZGEQRF (&nM, &nN, mQ._pd(), mQ._pld(), vTau, &dWork, &lWork, &nOutInfo);
    lWork = static_cast<int>(dWork.real());
    basic_cvector<double, std::complex<double> > vWork(lWork);

    // now do most of hardwork, find R and TAU and Householderish bits of Q
    ZGEQRF (&nM, &nN, mQ._pd(), mQ._pld(), vTau, vWork, &lWork, &nOutInfo);
    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

    // get upper-triangular R from overwritten A
    mR.vanish();
    for (int row = 1; row <= nK; ++row)
        for (int col = row; col <= nN; ++col)
            mR(row,col) = mQ(row,col);

    // calculate size of workspace for finding Q
    lWork = -1;
    ZUNGQR (&nM, &nK, &nK, mQ._pd(), mQ._pld(), vTau, &dWork, &lWork, &nOutInfo);
    lWork = static_cast<int>(dWork.real());
    if (lWork > vWork.size()) vWork.resize(lWork);

    // find Q
    ZUNGQR (&nM, &nK, &nK, mQ._pd(), mQ._pld(), vTau, vWork, &lWork, &nOutInfo);
    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
}

// Case 2: full mode, A is (m x n) and Q is (m x m) and R is (m x n)
template <>
CVM_API void 
__qrf<basic_cmatrix<float, std::complex<float> >, basic_scmatrix<float, std::complex<float> > >
    (const basic_cmatrix<float, std::complex<float> >& mArg, 
    basic_scmatrix<float, std::complex<float> >& mQ,
    basic_cmatrix<float, std::complex<float> >& mR) throw (cvmexception)
{
    const int nM = mArg.msize();
    const int nN = mArg.nsize();
    const int nK = _cvm_min<int>(nM, nN);

    // unlike economy mode, we need a copy here since Q will be m x m and this may bigger than original A
    basic_cmatrix<float, std::complex<float> > mA (nM, nN <= nM ? nM : nN);

    // copy over argument matrix
    mA.assign (1, 1, mArg);

    basic_cvector<float, std::complex<float> > vTau (nK);

    int row, col;
    int lWork = -1;
    int nOutInfo = 0;
    std::complex<float> dWork;

    // calculate size of workspace
    CGEQRF (&nM, &nN, mA, mA._pld(), vTau, &dWork, &lWork, &nOutInfo);
    lWork = static_cast<int>(dWork.real());
    basic_cvector<float, std::complex<float> > vWork(lWork);

    // now do most of hardwork, find R and TAU and Householderish bits of Q
    CGEQRF (&nM, &nN, mA, mA._pld(), vTau, vWork, &lWork, &nOutInfo);
    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

    // get upper-triangular R which is now m x n from overwritten A
    mR.vanish();
    for (row = 1; row <= nK; ++row)
        for (col = row; col <= nN; ++col)
            mR(row,col) = mA(row,col);

    // calculate size of workspace for finding Q that is m x m
    lWork = -1;
    CUNGQR (&nM, &nM, &nK, mA, mA._pld(), vTau, &dWork, &lWork, &nOutInfo);
    lWork = static_cast<int>(dWork.real());
    if (lWork > vWork.size()) vWork.resize(lWork);

    // find Q that is m x m
    CUNGQR (&nM, &nM, &nK, mA, mA._pld(), vTau, vWork, &lWork, &nOutInfo);
    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

    // is Q big enough to have all conents of mA ?
    if (nN <= nM)
        mQ.assign (1, 1, mA);
    else
        for (row = 1; row <= nM; ++row)
            for (col = 1; col <= nM; ++col)
                mQ(row,col) = mA(row,col);
}

template <>
CVM_API void 
__qrf<basic_cmatrix<double, std::complex<double> >, basic_scmatrix<double, std::complex<double> > >
    (const basic_cmatrix<double, std::complex<double> >& mArg, 
    basic_scmatrix<double, std::complex<double> >& mQ,
    basic_cmatrix<double, std::complex<double> >& mR) throw (cvmexception)
{
    const int nM = mArg.msize();
    const int nN = mArg.nsize();
    const int nK = _cvm_min<int>(nM, nN);

    // unlike economy mode, we need a copy here since Q will be m x m and this may bigger than original A
    basic_cmatrix<double, std::complex<double> > mA (nM, nN <= nM ? nM : nN);

    // copy over argument matrix
    mA.assign (1, 1, mArg);

    basic_cvector<double, std::complex<double> > vTau (nK);

    int row, col;
    int lWork = -1;
    int nOutInfo = 0;
    std::complex<double> dWork;

    // calculate size of workspace
    ZGEQRF (&nM, &nN, mA, mA._pld(), vTau, &dWork, &lWork, &nOutInfo);
    lWork = static_cast<int>(dWork.real());
    basic_cvector<double, std::complex<double> > vWork(lWork);

    // now do most of hardwork, find R and TAU and Householderish bits of Q
    ZGEQRF (&nM, &nN, mA, mA._pld(), vTau, vWork, &lWork, &nOutInfo);
    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

    // get upper-triangular R which is now m x n from overwritten A
    mR.vanish();
    for (row = 1; row <= nK; ++row)
        for (col = row; col <= nN; ++col)
            mR(row,col) = mA(row,col);

    // calculate size of workspace for finding Q that is m x m
    lWork = -1;
    ZUNGQR (&nM, &nM, &nK, mA, mA._pld(), vTau, &dWork, &lWork, &nOutInfo);
    lWork = static_cast<int>(dWork.real());
    if (lWork > vWork.size()) vWork.resize(lWork);

    // find Q that is m x m
    ZUNGQR (&nM, &nM, &nK, mA, mA._pld(), vTau, vWork, &lWork, &nOutInfo);
    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

    // is Q big enough to have all conents of mA ?
    if (nN <= nM)
        mQ.assign (1, 1, mA);
    else
        for (row = 1; row <= nM; ++row)
            for (col = 1; col <= nM; ++col)
                mQ(row,col) = mA(row,col);
}

CVM_NAMESPACE_END
