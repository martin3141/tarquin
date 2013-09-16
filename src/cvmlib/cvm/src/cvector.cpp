/* CVM Class Library */
/* Copyright (C), Sergei Nikolaev, 1992-2006, http://cvmlib.com */

#include "cvm.h"
#include "blas.h"

#if defined (_MSC_VER)
#   pragma warning(disable:4786)
#endif

CVM_NAMESPACE_BEG

template<>
CVM_API std::complex<float> __dotu<std::complex<float> > (const std::complex<float>* mpD, int mnSize, int mnIncr, const std::complex<float>* pD, int nIncr)
{
    std::complex<float> cRes;
    VCDOTU (&cRes, &mnSize, mpD, &mnIncr, pD, &nIncr);
    return cRes;
}

template<>
CVM_API std::complex<double> __dotu<std::complex<double> > (const std::complex<double>* mpD, int mnSize, int mnIncr, const std::complex<double>* pD, int nIncr)
{
    std::complex<double> cRes;
    VZDOTU (&cRes, &mnSize, mpD, &mnIncr, pD, &nIncr);
    return cRes;
}

template<>
CVM_API std::complex<float> __dotc<std::complex<float> > (const std::complex<float>* mpD, int mnSize, int mnIncr, const std::complex<float>* pD, int nIncr)
{
    std::complex<float> cRes;
    VCDOTC (&cRes, &mnSize, mpD, &mnIncr, pD, &nIncr);
    return cRes;
}

template<>
CVM_API std::complex<double> __dotc<std::complex<double> > (const std::complex<double>* mpD, int mnSize, int mnIncr, const std::complex<double>* pD, int nIncr)
{
    std::complex<double> cRes;
    VZDOTC (&cRes, &mnSize, mpD, &mnIncr, pD, &nIncr);
    return cRes;
}

template<>
CVM_API void
__gemv<std::complex<float>, basic_cmatrix<float, std::complex<float> >, basic_cvector<float, std::complex<float> > >
    (bool bLeft,
    const basic_cmatrix<float, std::complex<float> >& m,
    std::complex<float> dAlpha,
    const basic_cvector<float, std::complex<float> >& v,
    std::complex<float> dBeta,
    basic_cvector<float, std::complex<float> >& vRes)
{
    CGEMV (bLeft ? Chars::pT() : Chars::pN(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
           1,
#endif
           m._pm(), m._pn(), &dAlpha, m._pd(), m._pldm(), v, v._pincr(), &dBeta, vRes, vRes._pincr());
}

template<>
CVM_API void
__gemv<std::complex<double>, basic_cmatrix<double, std::complex<double> >, basic_cvector<double, std::complex<double> > >
    (bool bLeft,
    const basic_cmatrix<double, std::complex<double> >& m,
    std::complex<double> dAlpha,
    const basic_cvector<double, std::complex<double> >& v,
    std::complex<double> dBeta,
    basic_cvector<double, std::complex<double> >& vRes)
{
    ZGEMV (bLeft ? Chars::pT() : Chars::pN(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
           1,
#endif
           m._pm(), m._pn(), &dAlpha, m._pd(), m._pldm(), v, v._pincr(), &dBeta, vRes, vRes._pincr());
}

template<>
CVM_API void
__gbmv<std::complex<float>, basic_scbmatrix<float, std::complex<float> >, basic_cvector<float, std::complex<float> > >
    (bool bLeft,
    const basic_scbmatrix<float, std::complex<float> >& m,
    std::complex<float> dAlpha,
    const basic_cvector<float, std::complex<float> >& v,
    std::complex<float> dBeta,
    basic_cvector<float, std::complex<float> >& vRes)
{
    CGBMV (bLeft ? Chars::pT() : Chars::pN(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
           1,
#endif
           m._pm(), m._pn(), m._pl(), m._pu(), &dAlpha, m, m._pld(), v, v._pincr(), &dBeta, vRes, vRes._pincr());
}

template<>
CVM_API void
__gbmv<std::complex<double>, basic_scbmatrix<double, std::complex<double> >, basic_cvector<double, std::complex<double> > >
    (bool bLeft,
    const basic_scbmatrix<double, std::complex<double> >& m,
    std::complex<double> dAlpha,
    const basic_cvector<double, std::complex<double> >& v,
    std::complex<double> dBeta,
    basic_cvector<double, std::complex<double> >& vRes)
{
    ZGBMV (bLeft ? Chars::pT() : Chars::pN(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
           1,
#endif
           m._pm(), m._pn(), m._pl(), m._pu(), &dAlpha, m, m._pld(), v, v._pincr(), &dBeta, vRes, vRes._pincr());
}

template<>
CVM_API void
__shmv<std::complex<float>, basic_schmatrix<float, std::complex<float> >, basic_cvector<float, std::complex<float> > >
    (const basic_schmatrix<float, std::complex<float> >& m,
     std::complex<float> cAlpha,
     const basic_cvector<float, std::complex<float> >& v,
     std::complex<float> cBeta,
     basic_cvector<float, std::complex<float> >& vRes)
{
    CHEMV (Chars::pU(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
           1,
#endif
           m._pm(), &cAlpha, m, m._pld(), v, v._pincr(), &cBeta, vRes, vRes._pincr());
}

template<>
CVM_API void
__shmv<std::complex<double>, basic_schmatrix<double, std::complex<double> >, basic_cvector<double, std::complex<double> > >
    (const basic_schmatrix<double, std::complex<double> >& m,
     std::complex<double> cAlpha,
     const basic_cvector<double, std::complex<double> >& v,
     std::complex<double> cBeta,
     basic_cvector<double, std::complex<double> >& vRes)
{
    ZHEMV (Chars::pU(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
           1,
#endif
           m._pm(), &cAlpha, m, m._pld(), v, v._pincr(), &cBeta, vRes, vRes._pincr());
}


template<>
CVM_API void
__eig<basic_cvector<float, std::complex<float> >, basic_srmatrix<float>, basic_scmatrix<float, std::complex<float> > >
    (basic_cvector<float, std::complex<float> >& vRes,
    const basic_srmatrix<float>& mArg,
    basic_scmatrix<float, std::complex<float> >* mEigVect,
    bool bRightVect) throw (cvmexception)
{
    const bool bEigVect = (mEigVect != NULL);
    const int  nM       = mArg.msize();
    const int  lWork    = nM * 64;
          int  ilo      = 0;
          int  ihi      = 0;
          int  nOutInfo = 0;

    if (nM != vRes.size()) throw cvmexception (CVM_SIZESMISMATCH);
    if (nM == 1)
    {
        static const std::complex<float> one(1.);
        vRes[1] = std::complex<float>(mArg(1,1), 0.F);
        if (bEigVect)
        {
            mEigVect -> resize (1);
            (*mEigVect)[1].set(one);
        }
    }
    else
    {
        basic_srmatrix<float> mA     (mArg);
        basic_rvector<float>  vScale (nM);
        basic_rvector<float>  vTau   (_cvm_max<int>(1, nM - 1));
        basic_rvector<float>  vWork  (lWork);
        basic_rvector<float>  vR     (nM);
        basic_rvector<float>  vI     (nM);

        SGEBAL (Chars::pB(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                1,
#endif
                &nM, mA, &nM, &ilo, &ihi, vScale, &nOutInfo);
        if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

        SGEHRD (&nM, &ilo, &ihi, mA, &nM, vTau, vWork, &lWork, &nOutInfo);
        if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

        if (bEigVect)
        {
            static const int one(1);
            int       m       = 0;
            int       lSelect = 0;
            const int ldvl    = bRightVect ? 1 : nM;
            const int ldvr    = bRightVect ? nM : 1;
            basic_srmatrix<float> vl   (ldvl);
            basic_srmatrix<float> vr   (ldvr);
            basic_rvector <float> work (3 * nM);
            basic_srmatrix<float> mH   (mA);
            const char* pRL = bRightVect ? Chars::pR() : Chars::pL();

            if (bRightVect)
            {
                vr = mA;
            }
            else
            {
                vl = mA;
            }

            SORGHR (&nM, &one, &nM, bRightVect ? vr : vl, &nM, vTau, vWork, &lWork, &nOutInfo);
            if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
            SHSEQR (Chars::pS(), 1, Chars::pV(), 1,
#else
            SHSEQR (Chars::pS(),    Chars::pV(),
#endif
                    &nM, &ilo, &ihi, mH, &nM, vR, vI, bRightVect ? vr : vl, &nM, vWork, &lWork, &nOutInfo);
            if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
            if (nOutInfo > 0) throw cvmexception (CVM_CONVERGENCE_ERROR);

#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
            STREVC (pRL, 1, Chars::pB(), 1,
#else
            STREVC (pRL,    Chars::pB(),
#endif
                    &lSelect, &nM, mH, &nM, vl, &ldvl, vr, &ldvr, &nM, &m, work, &nOutInfo);
            if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

            basic_srmatrix<float>& v = bRightVect ? vr   : vl;
            const int ldv = bRightVect ? ldvr : ldvl;

#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
            SGEBAK (Chars::pB(), 1, pRL, 1,
#else
            SGEBAK (Chars::pB(),    pRL,
#endif
                    &nM, &ilo, &ihi, vScale, &nM, v, &ldv, &nOutInfo);
            if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

            bool bPair = false;
            m = 1;
            mEigVect -> resize (nM);
            for (int i = 1; i <= nM; i++)
            {
                if (_abs(vI(i)) > basic_cvmMachMin<float>())
                {
                    (*mEigVect)(i).assign_real (v (m));

                    if (bPair)
                    {
                        (*mEigVect)(i).assign_imag (- v (m + 1));
                        m += 2;
                        bPair = false;
                    }
                    else
                    {
                        (*mEigVect)(i).assign_imag (v (m + 1));
                        bPair = true;
                    }
                }
                else
                {
                    static const float zero(0.);
                    (*mEigVect)(i).assign_real (v (m));
                    (*mEigVect)(i).set_imag (zero);
                    m++;
                }
            }
        }
        else
        {
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
            SHSEQR (Chars::pE(), 1, Chars::pN(), 1,
#else
            SHSEQR (Chars::pE(),    Chars::pN(),
#endif
                    &nM, &ilo, &ihi, mA, &nM, vR, vI, NULL, &nM, vWork, &lWork, &nOutInfo);
            if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
            if (nOutInfo > 0) throw cvmexception (CVM_CONVERGENCE_ERROR);
        }

        vRes.assign_real (vR);
        vRes.assign_imag (vI);
    }
}

template<>
CVM_API void
__eig<basic_cvector<double, std::complex<double> >, basic_srmatrix<double>, basic_scmatrix<double, std::complex<double> > >
    (basic_cvector<double, std::complex<double> >& vRes,
    const basic_srmatrix<double>& mArg,
    basic_scmatrix<double, std::complex<double> >* mEigVect,
    bool bRightVect) throw (cvmexception)
{
    const bool bEigVect = (mEigVect != NULL);
    const int  nM       = mArg.msize();
    const int  lWork    = nM * 64;
          int  ilo      = 0;
          int  ihi      = 0;
          int  nOutInfo = 0;

    if (nM != vRes.size()) throw cvmexception (CVM_SIZESMISMATCH);
    if (nM == 1)
    {
        vRes[1] = std::complex<double>(mArg(1,1), 0.);
        if (bEigVect)
        {
            static const std::complex<double> one(1.);
            mEigVect -> resize (1);
            (*mEigVect)[1].set(one);
        }
    }
    else
    {
        basic_srmatrix<double> mA     (mArg);
        basic_rvector<double>  vScale (nM);
        basic_rvector<double>  vTau   (_cvm_max<int>(1, nM - 1));
        basic_rvector<double>  vWork  (lWork);
        basic_rvector<double>  vR     (nM);
        basic_rvector<double>  vI     (nM);

        DGEBAL (Chars::pB(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                1,
#endif
                &nM, mA, &nM, &ilo, &ihi, vScale, &nOutInfo);
        if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

        DGEHRD (&nM, &ilo, &ihi, mA, &nM, vTau, vWork, &lWork, &nOutInfo);
        if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

        if (bEigVect)
        {
            static const int one(1);
            int       m       = 0;
            int       lSelect = 0;
            const int ldvl    = bRightVect ? 1 : nM;
            const int ldvr    = bRightVect ? nM : 1;
            basic_srmatrix<double> vl   (ldvl);
            basic_srmatrix<double> vr   (ldvr);
            basic_rvector <double> work (3 * nM);
            basic_srmatrix<double> mH   (mA);
            const char* pRL = bRightVect ? Chars::pR() : Chars::pL();

            if (bRightVect)
            {
                vr = mA;
            }
            else
            {
                vl = mA;
            }

            DORGHR (&nM, &one, &nM, bRightVect ? vr : vl, &nM, vTau, vWork, &lWork, &nOutInfo);
            if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
            DHSEQR (Chars::pS(), 1, Chars::pV(), 1,
#else
            DHSEQR (Chars::pS(),    Chars::pV(),
#endif
                    &nM, &ilo, &ihi, mH, &nM, vR, vI, bRightVect ? vr : vl, &nM, vWork, &lWork, &nOutInfo);
            if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
            if (nOutInfo > 0) throw cvmexception (CVM_CONVERGENCE_ERROR);

#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
            DTREVC (pRL, 1, Chars::pB(), 1,
#else
            DTREVC (pRL,    Chars::pB(),
#endif
                    &lSelect, &nM, mH, &nM, vl, &ldvl, vr, &ldvr, &nM, &m, work, &nOutInfo);
            if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

            basic_srmatrix<double>& v = bRightVect ? vr   : vl;
            const int ldv = bRightVect ? ldvr : ldvl;

#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
            DGEBAK (Chars::pB(), 1, pRL, 1,
#else
            DGEBAK (Chars::pB(),    pRL,
#endif
                    &nM, &ilo, &ihi, vScale, &nM, v, &ldv, &nOutInfo);
            if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

            bool bPair = false;
            m = 1;
            mEigVect -> resize (nM);
            for (int i = 1; i <= nM; i++)
            {
                if (_abs(vI(i)) > basic_cvmMachMin<double>())
                {
                    (*mEigVect)(i).assign_real (v (m));

                    if (bPair)
                    {
                        (*mEigVect)(i).assign_imag (- v (m + 1));
                        m += 2;
                        bPair = false;
                    }
                    else
                    {
                        (*mEigVect)(i).assign_imag (v (m + 1));
                        bPair = true;
                    }
                }
                else
                {
                    static const double zero(0.);
                    (*mEigVect)(i).assign_real (v (m));
                    (*mEigVect)(i).set_imag (zero);
                    m++;
                }
            }
        }
        else
        {
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
            DHSEQR (Chars::pE(), 1, Chars::pN(), 1,
#else
            DHSEQR (Chars::pE(), Chars::pN(),
#endif
                    &nM, &ilo, &ihi, mA, &nM, vR, vI, NULL, &nM, vWork, &lWork, &nOutInfo);
            if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
            if (nOutInfo > 0) throw cvmexception (CVM_CONVERGENCE_ERROR);
        }

        vRes.assign_real (vR);
        vRes.assign_imag (vI);
    }
}

template<>
CVM_API void
__eig<basic_cvector<float, std::complex<float> >, basic_scmatrix<float, std::complex<float> >, basic_scmatrix<float, std::complex<float> > >
    (basic_cvector<float, std::complex<float> >& vRes,
    const basic_scmatrix<float, std::complex<float> >& mArg,
    basic_scmatrix<float, std::complex<float> >* mEigVect,
    bool bRightVect) throw (cvmexception)
{
    const bool bEigVect = (mEigVect != NULL);
    const int  nM       = mArg.msize();
    int  lWork    = -1;
    int  ilo      = 0;
    int  ihi      = 0;
    int  nOutInfo = 0;

    if (nM != vRes.size()) throw cvmexception (CVM_SIZESMISMATCH);
    if (nM == 1)
    {
        vRes[1] = mArg(1,1);
        if (bEigVect)
        {
            static const std::complex<float> one(1.F);
            mEigVect -> resize (1);
            (*mEigVect)[1].set(one);
        }
    }
    else
    {
        basic_scmatrix<float, std::complex<float> > mA (mArg);
        basic_rvector<float> vScale (nM);
        basic_cvector<float, std::complex<float> > vTau (_cvm_max<int>(1, nM - 1));
        basic_cvector<float, std::complex<float> > vWork (nM * 64);
        basic_cvector<float, std::complex<float> > vW (nM);

        CGEBAL (Chars::pB(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                1,
#endif
                &nM, mA, &nM, &ilo, &ihi, vScale, &nOutInfo);
        if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

        CGEHRD (&nM, &ilo, &ihi, mA, &nM, vTau, vWork, &lWork, &nOutInfo);
        lWork = static_cast<int> (vWork[1].real());
        if (lWork > vWork.size()) vWork.resize(lWork);

        CGEHRD (&nM, &ilo, &ihi, mA, &nM, vTau, vWork, &lWork, &nOutInfo);
        if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

        if (bEigVect)
        {
            static const int one(1);
            const int ldvl    = bRightVect ? 1 : nM;
            const int ldvr    = bRightVect ? nM : 1;
            int       m       = 0;
            int       lSelect = 0;
            basic_scmatrix<float, std::complex<float> > vl    (ldvl);
            basic_scmatrix<float, std::complex<float> > vr    (ldvr);
            basic_cvector <float, std::complex<float> > work  (2 * nM);
            basic_rvector <float> rwork (nM);
            const char* pRL = bRightVect ? Chars::pR() : Chars::pL();

            if (bRightVect)
            {
                vr = mA;
            }
            else
            {
                vl = mA;
            }

            lWork = -1;
            CUNGHR (&nM, &one, &nM, bRightVect ? vr : vl, &nM, vTau, vWork, &lWork, &nOutInfo);
            lWork = static_cast<int> (vWork[1].real());
            if (lWork > vWork.size()) vWork.resize(lWork);

            CUNGHR (&nM, &one, &nM, bRightVect ? vr : vl, &nM, vTau, vWork, &lWork, &nOutInfo);
            if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

            lWork = -1;
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
            CHSEQR (Chars::pS(), 1, Chars::pV(), 1,
#else
            CHSEQR (Chars::pS(),    Chars::pV(),
#endif
                    &nM, &ilo, &ihi, mA, &nM, vW, bRightVect ? vr : vl, &nM, vWork, &lWork, &nOutInfo);
            lWork = static_cast<int> (vWork[1].real());
            if (lWork > vWork.size()) vWork.resize(lWork);

#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
            CHSEQR (Chars::pS(), 1, Chars::pV(), 1,
#else
            CHSEQR (Chars::pS(),    Chars::pV(),
#endif
                    &nM, &ilo, &ihi, mA, &nM, vW, bRightVect ? vr : vl, &nM, vWork, &lWork, &nOutInfo);
            if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
            if (nOutInfo > 0) throw cvmexception (CVM_CONVERGENCE_ERROR);

#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
            CTREVC (pRL, 1, Chars::pB(), 1,
#else
            CTREVC (pRL,    Chars::pB(),
#endif
                    &lSelect, &nM, mA, &nM, vl, &ldvl, vr, &ldvr, &nM, &m, work, rwork, &nOutInfo);
            if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

            basic_scmatrix<float, std::complex<float> >& v = bRightVect ? vr : vl;
            const int ldv = bRightVect ? ldvr : ldvl;

#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
            CGEBAK (Chars::pB(), 1, pRL, 1,
#else
            CGEBAK (Chars::pB(),    pRL,
#endif
                    &nM, &ilo, &ihi, vScale, &nM, v, &ldv, &nOutInfo);
            if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

            (*mEigVect) << v;
        }
        else
        {
            lWork = -1;
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
            CHSEQR (Chars::pE(), 1, Chars::pN(), 1,
#else
            CHSEQR (Chars::pE(),    Chars::pN(),
#endif
                    &nM, &ilo, &ihi, mA, &nM, vW, NULL, &nM, vWork, &lWork, &nOutInfo);
            lWork = static_cast<int> (vWork[1].real());
            if (lWork > vWork.size()) vWork.resize(lWork);


#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
            CHSEQR (Chars::pE(), 1, Chars::pN(), 1,
#else
            CHSEQR (Chars::pE(),    Chars::pN(),
#endif
                    &nM, &ilo, &ihi, mA, &nM, vW, NULL, &nM, vWork, &lWork, &nOutInfo);
            if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
            if (nOutInfo > 0) throw cvmexception (CVM_CONVERGENCE_ERROR);
        }

        vRes.assign (vW, vW.incr());
    }
}

template<>
CVM_API void
__eig<basic_cvector<double, std::complex<double> >, basic_scmatrix<double, std::complex<double> >, basic_scmatrix<double, std::complex<double> > >
    (basic_cvector<double, std::complex<double> >& vRes,
    const basic_scmatrix<double, std::complex<double> >& mArg,
    basic_scmatrix<double, std::complex<double> >* mEigVect,
    bool bRightVect) throw (cvmexception)
{
    const bool bEigVect = (mEigVect != NULL);
    const int  nM       = mArg.msize();
    int  lWork    = -1;
    int  ilo      = 0;
    int  ihi      = 0;
    int  nOutInfo = 0;

    if (nM != vRes.size()) throw cvmexception (CVM_SIZESMISMATCH);
    if (nM == 1)
    {
        vRes[1] = mArg(1,1);
        if (bEigVect)
        {
            static const std::complex<double> one(1.);
            mEigVect -> resize (1);
            (*mEigVect)[1].set(one);
        }
    }
    else
    {
        basic_scmatrix<double, std::complex<double> > mA (mArg);
        basic_rvector<double> vScale (nM);
        basic_cvector<double, std::complex<double> > vTau (_cvm_max<int>(1, nM - 1));
        basic_cvector<double, std::complex<double> > vWork (nM * 64);
        basic_cvector<double, std::complex<double> > vW (nM);

        ZGEBAL (Chars::pB(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                1,
#endif
                &nM, mA, &nM, &ilo, &ihi, vScale, &nOutInfo);
        if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

        ZGEHRD (&nM, &ilo, &ihi, mA, &nM, vTau, vWork, &lWork, &nOutInfo);
        lWork = static_cast<int> (vWork[1].real());
        if (lWork > vWork.size()) vWork.resize(lWork);

        ZGEHRD (&nM, &ilo, &ihi, mA, &nM, vTau, vWork, &lWork, &nOutInfo);
        if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

        if (bEigVect)
        {
            static const int one(1);
            const int ldvl    = bRightVect ? 1 : nM;
            const int ldvr    = bRightVect ? nM : 1;
            int       m       = 0;
            int       lSelect = 0;
            basic_scmatrix<double, std::complex<double> > vl    (ldvl);
            basic_scmatrix<double, std::complex<double> > vr    (ldvr);
            basic_cvector <double, std::complex<double> > work  (2 * nM);
            basic_rvector <double> rwork (nM);
            const char* pRL = bRightVect ? Chars::pR() : Chars::pL();

            if (bRightVect)
            {
                vr = mA;
            }
            else
            {
                vl = mA;
            }

            lWork = -1;
            ZUNGHR (&nM, &one, &nM, bRightVect ? vr : vl, &nM, vTau, vWork, &lWork, &nOutInfo);
            lWork = static_cast<int> (vWork[1].real());
            if (lWork > vWork.size()) vWork.resize(lWork);

            ZUNGHR (&nM, &one, &nM, bRightVect ? vr : vl, &nM, vTau, vWork, &lWork, &nOutInfo);
            if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

            lWork = -1;
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
            ZHSEQR (Chars::pS(), 1, Chars::pV(), 1,
#else
            ZHSEQR (Chars::pS(),    Chars::pV(),
#endif
                    &nM, &ilo, &ihi, mA, &nM, vW, bRightVect ? vr : vl, &nM, vWork, &lWork, &nOutInfo);
            lWork = static_cast<int> (vWork[1].real());
            if (lWork > vWork.size()) vWork.resize(lWork);

#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
            ZHSEQR (Chars::pS(), 1, Chars::pV(), 1,
#else
            ZHSEQR (Chars::pS(),    Chars::pV(),
#endif
                    &nM, &ilo, &ihi, mA, &nM, vW, bRightVect ? vr : vl, &nM, vWork, &lWork, &nOutInfo);
            if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
            if (nOutInfo > 0) throw cvmexception (CVM_CONVERGENCE_ERROR);

#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
            ZTREVC (pRL, 1, Chars::pB(), 1,
#else
            ZTREVC (pRL,    Chars::pB(),
#endif
                    &lSelect, &nM, mA, &nM, vl, &ldvl, vr, &ldvr, &nM, &m, work, rwork, &nOutInfo);
            if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

            basic_scmatrix<double, std::complex<double> >& v = bRightVect ? vr : vl;
            const int ldv = bRightVect ? ldvr : ldvl;

#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
            ZGEBAK (Chars::pB(), 1, pRL, 1,
#else
            ZGEBAK (Chars::pB(),    pRL,
#endif
                    &nM, &ilo, &ihi, vScale, &nM, v, &ldv, &nOutInfo);
            if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

            (*mEigVect) << v;
        }
        else
        {
            lWork = -1;
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
            ZHSEQR (Chars::pE(), 1, Chars::pN(), 1,
#else
            ZHSEQR (Chars::pE(),    Chars::pN(),
#endif
                    &nM, &ilo, &ihi, mA, &nM, vW, NULL, &nM, vWork, &lWork, &nOutInfo);
            lWork = static_cast<int> (vWork[1].real());
            if (lWork > vWork.size()) vWork.resize(lWork);

#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
            ZHSEQR (Chars::pE(), 1, Chars::pN(), 1,
#else
            ZHSEQR (Chars::pE(),    Chars::pN(),
#endif
                    &nM, &ilo, &ihi, mA, &nM, vW, NULL, &nM, vWork, &lWork, &nOutInfo);
            if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
            if (nOutInfo > 0) throw cvmexception (CVM_CONVERGENCE_ERROR);
        }

        vRes.assign (vW, vW.incr());
    }
}

CVM_NAMESPACE_END
