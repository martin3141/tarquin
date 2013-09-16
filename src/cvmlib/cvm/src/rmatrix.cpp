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
__gemm<float, basic_rmatrix<float> >
    (const basic_rmatrix<float>& ml, bool bTrans1,
     const basic_rmatrix<float>& mr, bool bTrans2,
     float dAlpha,
     basic_rmatrix<float>& mRes,
     float dBeta)
{
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
    SGEMM (bTrans1 ? Chars::pT() : Chars::pN(), 1, bTrans2 ? Chars::pT() : Chars::pN(), 1,
#else
    SGEMM (bTrans1 ? Chars::pT() : Chars::pN(),    bTrans2 ? Chars::pT() : Chars::pN(),
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
__gemm<double, basic_rmatrix<double> >
    (const basic_rmatrix<double>& ml, bool bTrans1,
     const basic_rmatrix<double>& mr, bool bTrans2,
     double dAlpha,
     basic_rmatrix<double>& mRes,
     double dBeta)
{
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
    DGEMM (bTrans1 ? Chars::pT() : Chars::pN(), 1, bTrans2 ? Chars::pT() : Chars::pN(), 1,
#else
    DGEMM (bTrans1 ? Chars::pT() : Chars::pN(),    bTrans2 ? Chars::pT() : Chars::pN(),
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
__symm<float, basic_srsmatrix<float>, basic_rmatrix<float> >
    (bool bLeft,
     const basic_srsmatrix<float>& ml,
     const basic_rmatrix<float>& mr, 
     float dAlpha,
     basic_rmatrix<float>& mRes,
     float dBeta)
{
    SSYMM (bLeft ? Chars::pL() : Chars::pR(),
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
__symm<double, basic_srsmatrix<double>, basic_rmatrix<double>  >
    (bool bLeft,
     const basic_srsmatrix<double>& ml,
     const basic_rmatrix<double>& mr, 
     double dAlpha,
     basic_rmatrix<double>& mRes,
     double dBeta)
{
    DSYMM (bLeft ? Chars::pL() : Chars::pR(),
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
__syrk<float, basic_srsmatrix<float> >
    (bool bTransp, 
    float alpha, int k,
    const float* pA, int ldA,
    float beta, basic_srsmatrix<float>& m)
{
    SSYRK (Chars::pU(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
           1,
#endif
           bTransp ? Chars::pT() : Chars::pN(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
           1,
#endif
           m._pm(), &k, &alpha, pA, &ldA, &beta, m, m._pld());
}

template <>
CVM_API void
__syrk<double, basic_srsmatrix<double> >
    (bool bTransp, 
    double alpha, int k,
    const double* pA, int ldA,
    double beta, basic_srsmatrix<double>& m)
{
    DSYRK (Chars::pU(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
           1,
#endif
           bTransp ? Chars::pT() : Chars::pN(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
           1,
#endif
           m._pm(), &k, &alpha, pA, &ldA, &beta, m, m._pld());
}

template <>
CVM_API void
__syr2k<float, basic_srsmatrix<float> >
    (bool bTransp, 
    float alpha, int k,
    const float* pA, int ldA,
    const float* pB, int ldB,
    float beta, basic_srsmatrix<float>& m)
{
    SSYR2K (Chars::pU(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
           1,
#endif
           bTransp ? Chars::pT() : Chars::pN(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
           1,
#endif
           m._pm(), &k, &alpha, pA, &ldA, pB, &ldB, &beta, m, m._pld());
}

template <>
CVM_API void
__syr2k<double, basic_srsmatrix<double> >
    (bool bTransp, 
    double alpha, int k,
    const double* pA, int ldA,
    const double* pB, int ldB,
    double beta, basic_srsmatrix<double>& m)
{
    DSYR2K (Chars::pU(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
           1,
#endif
           bTransp ? Chars::pT() : Chars::pN(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
           1,
#endif
           m._pm(), &k, &alpha, pA, &ldA, pB, &ldB, &beta, m, m._pld());
}


// Case 1:
// call *GEQRF to get R and TAU
// call *ORGQR to get Q (using TAU)

// Case 1: economy mode, A is (m x n) and Q is (m x n) and R is (n x n)
template <>
CVM_API void 
__qre<basic_rmatrix<float>, basic_srmatrix<float> >
    (const basic_rmatrix<float>& mArg, 
    basic_rmatrix<float>& mQ, 
    basic_srmatrix<float>& mR) throw (cvmexception)
{
    const int nM = mArg.msize();
    const int nN = mArg.nsize();
    const int nK = _cvm_min<int>(nM, nN);

    // we will eventually overwrite mQ to be the output matrix
    mQ = mArg;
    basic_rvector<float> vTau (nK);

    int lWork = -1;
    int nOutInfo = 0;
    float dWork;

    // calculate size of workspace
    SGEQRF (&nM, &nN, mQ._pd(), mQ._pld(), vTau, &dWork, &lWork, &nOutInfo);
    lWork = static_cast<int>(dWork);
    basic_rvector<float> vWork(lWork);

    // now do most of hardwork, find R and TAU and Householderish bits of Q
    SGEQRF (&nM, &nN, mQ._pd(), mQ._pld(), vTau, vWork, &lWork, &nOutInfo);
    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

    // get upper-triangular R from overwritten A
    mR.vanish();
    for (int row = 1; row <= nK; ++row)
        for (int col = row; col <= nN; ++col)
            mR(row,col) = mQ(row,col);

    // calculate size of workspace for finding Q
    lWork = -1;
    SORGQR (&nM, &nK, &nK, mQ._pd(), mQ._pld(), vTau, &dWork, &lWork, &nOutInfo);
    lWork = static_cast<int>(dWork);
    if (lWork > vWork.size()) vWork.resize(lWork);

    // find Q
    SORGQR (&nM, &nK, &nK, mQ._pd(), mQ._pld(), vTau, vWork, &lWork, &nOutInfo);
    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
}

template <>
CVM_API void 
__qre<basic_rmatrix<double>, basic_srmatrix<double> >
    (const basic_rmatrix<double>& mArg, 
    basic_rmatrix<double>& mQ, 
    basic_srmatrix<double>& mR) throw (cvmexception)
{
    const int nM = mArg.msize();
    const int nN = mArg.nsize();
    const int nK = _cvm_min<int>(nM, nN);

    // we will eventually overwrite mQ to be the output matrix
    mQ = mArg;
    basic_rvector<double> vTau (nK);

    int lWork = -1;
    int nOutInfo = 0;
    double dWork;

    // calculate size of workspace
    DGEQRF (&nM, &nN, mQ._pd(), mQ._pld(), vTau, &dWork, &lWork, &nOutInfo);
    lWork = static_cast<int>(dWork);
    basic_rvector<double> vWork(lWork);

    // now do most of hardwork, find R and TAU and Householderish bits of Q
    DGEQRF (&nM, &nN, mQ._pd(), mQ._pld(), vTau, vWork, &lWork, &nOutInfo);
    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

    // get upper-triangular R from overwritten A
    mR.vanish();
    for (int row = 1; row <= nK; ++row)
        for (int col = row; col <= nN; ++col)
            mR(row,col) = mQ(row,col);

    // calculate size of workspace for finding Q
    lWork = -1;
    DORGQR (&nM, &nK, &nK, mQ._pd(), mQ._pld(), vTau, &dWork, &lWork, &nOutInfo);
    lWork = static_cast<int>(dWork);
    if (lWork > vWork.size()) vWork.resize(lWork);

    // find Q
    DORGQR (&nM, &nK, &nK, mQ._pd(), mQ._pld(), vTau, vWork, &lWork, &nOutInfo);
    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
}

// Case 2: full mode, A is (m x n) and Q is (m x m) and R is (m x n)
template <>
CVM_API void 
__qrf<basic_rmatrix<float>, basic_srmatrix<float> >
    (const basic_rmatrix<float>& mArg, 
    basic_srmatrix<float>& mQ, 
    basic_rmatrix<float>& mR) throw (cvmexception)
{
    const int nM = mArg.msize();
    const int nN = mArg.nsize();
    const int nK = _cvm_min<int>(nM, nN);

    // unlike economy mode, we need a copy here since Q will be m x m and this may bigger than original A
    basic_rmatrix<float> mA (nM, nN <= nM ? nM : nN);

    // copy over argument matrix
    mA.assign (1, 1, mArg);

    basic_rvector<float> vTau (nK);

    int row, col;
    int lWork = -1;
    int nOutInfo = 0;
    float dWork;

    // calculate size of workspace
    SGEQRF (&nM, &nN, mA, mA._pld(), vTau, &dWork, &lWork, &nOutInfo);
    lWork = static_cast<int>(dWork);
    basic_rvector<float> vWork(lWork);

    // now do most of hardwork, find R and TAU and Householderish bits of Q
    SGEQRF (&nM, &nN, mA, mA._pld(), vTau, vWork, &lWork, &nOutInfo);
    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

    // get upper-triangular R which is now m x n from overwritten A
    mR.vanish();
    for (row = 1; row <= nK; ++row)
        for (col = row; col <= nN; ++col)
            mR(row,col) = mA(row,col);

    // calculate size of workspace for finding Q that is m x m
    lWork = -1;
    SORGQR (&nM, &nM, &nK, mA, mA._pld(), vTau, &dWork, &lWork, &nOutInfo);
    lWork = static_cast<int>(dWork);
    if (lWork > vWork.size()) vWork.resize(lWork);

    // find Q that is m x m
    SORGQR (&nM, &nM, &nK, mA, mA._pld(), vTau, vWork, &lWork, &nOutInfo);
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
__qrf<basic_rmatrix<double>, basic_srmatrix<double> >
    (const basic_rmatrix<double>& mArg, 
    basic_srmatrix<double>& mQ, 
    basic_rmatrix<double>& mR) throw (cvmexception)
{
    const int nM = mArg.msize();
    const int nN = mArg.nsize();
    const int nK = _cvm_min<int>(nM, nN);

    // unlike economy mode, we need a copy here since Q will be m x m and this may bigger than original A
    basic_rmatrix<double> mA (nM, nN <= nM ? nM : nN);

    // copy over argument matrix
    mA.assign (1, 1, mArg);

    basic_rvector<double> vTau (nK);

    int row, col;
    int lWork = -1;
    int nOutInfo = 0;
    double dWork;

    // calculate size of workspace
    DGEQRF (&nM, &nN, mA, mA._pld(), vTau, &dWork, &lWork, &nOutInfo);
    lWork = static_cast<int>(dWork);
    basic_rvector<double> vWork(lWork);

    // now do most of hardwork, find R and TAU and Householderish bits of Q
    DGEQRF (&nM, &nN, mA, mA._pld(), vTau, vWork, &lWork, &nOutInfo);
    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

    // get upper-triangular R which is now m x n from overwritten A
    mR.vanish();
    for (row = 1; row <= nK; ++row)
        for (col = row; col <= nN; ++col)
            mR(row,col) = mA(row,col);

    // calculate size of workspace for finding Q that is m x m
    lWork = -1;
    DORGQR (&nM, &nM, &nK, mA, mA._pld(), vTau, &dWork, &lWork, &nOutInfo);
    lWork = static_cast<int>(dWork);
    if (lWork > vWork.size()) vWork.resize(lWork);

    // find Q that is m x m
    DORGQR (&nM, &nM, &nK, mA, mA._pld(), vTau, vWork, &lWork, &nOutInfo);
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
