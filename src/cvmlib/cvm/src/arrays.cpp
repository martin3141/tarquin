/* CVM Class Library */
/* Copyright (C), Sergei Nikolaev, 1992-2006, http://cvmlib.com */

#include "cvm.h"
#include "blas.h"

#if defined (_MSC_VER)
#   pragma warning(disable:4786)
#endif


CVM_NAMESPACE_BEG

template<>
CVM_API float __norm<float, float>(const float* pD, int nSize, int nIncr)
{
    CVM_ASSERT(pD, ((nSize - 1) * nIncr + 1) * sizeof(float))
    return SNRM2 (&nSize, pD, &nIncr);
}
template<>
CVM_API double __norm<double, double>(const double* pD, int nSize, int nIncr)
{
    CVM_ASSERT(pD, ((nSize - 1) * nIncr + 1) * sizeof(double))
    return DNRM2 (&nSize, pD, &nIncr);
}
template<>
CVM_API float __norm<float, std::complex<float> >(const std::complex<float>* pD, int nSize, int nIncr)
{
    CVM_ASSERT(pD, ((nSize - 1) * nIncr + 1) * sizeof(std::complex<float>))
    return SCNRM2 (&nSize, pD, &nIncr);
}
template<>
CVM_API double __norm<double, std::complex<double> >(const std::complex<double>* pD, int nSize, int nIncr)
{
    CVM_ASSERT(pD, ((nSize - 1) * nIncr + 1) * sizeof(std::complex<double>))
    return DZNRM2 (&nSize, pD, &nIncr);
}

template<>
CVM_API int __idamax<float>(const float* pD, int nSize, int nIncr)
{
    CVM_ASSERT(pD, ((nSize - 1) * nIncr + 1) * sizeof(float))
    return ISAMAX (&nSize, pD, &nIncr);
}
template<>
CVM_API int __idamax<double>(const double* pD, int nSize, int nIncr)
{
    CVM_ASSERT(pD, ((nSize - 1) * nIncr + 1) * sizeof(double))
    return IDAMAX (&nSize, pD, &nIncr);
}
template<>
CVM_API int __idamin<float>(const float* pD, int nSize, int nIncr)
{
    CVM_ASSERT(pD, ((nSize - 1) * nIncr + 1) * sizeof(float))
    return ISAMIN (&nSize, pD, &nIncr);
}
template<>
CVM_API int __idamin<double>(const double* pD, int nSize, int nIncr)
{
    CVM_ASSERT(pD, ((nSize - 1) * nIncr + 1) * sizeof(double))
    return IDAMIN (&nSize, pD, &nIncr);
}
template<>
CVM_API int __idamax<std::complex<float> >(const std::complex<float>* pD, int nSize, int nIncr)
{
    CVM_ASSERT(pD, ((nSize - 1) * nIncr + 1) * sizeof(std::complex<float>))
    return ICAMAX (&nSize, pD, &nIncr);
}
template<>
CVM_API int __idamax<std::complex<double> >(const std::complex<double>* pD, int nSize, int nIncr)
{
    CVM_ASSERT(pD, ((nSize - 1) * nIncr + 1) * sizeof(std::complex<double>))
    return IZAMAX (&nSize, pD, &nIncr);
}
template<>
CVM_API int __idamin<std::complex<float> >(const std::complex<float>* pD, int nSize, int nIncr)
{
    CVM_ASSERT(pD, ((nSize - 1) * nIncr + 1) * sizeof(std::complex<float>))
    return ICAMIN (&nSize, pD, &nIncr);
}
template<>
CVM_API int __idamin<std::complex<double> >(const std::complex<double>* pD, int nSize, int nIncr)
{
    CVM_ASSERT(pD, ((nSize - 1) * nIncr + 1) * sizeof(std::complex<double>))
    return IZAMIN (&nSize, pD, &nIncr);
}

template<>
CVM_API void __add<float>(float* mpD, int mnSize, int mnIncr, const float* pv, int nIncr)
{
    static const float one(1.F);
    CVM_ASSERT(mpD, ((mnSize - 1) * mnIncr + 1) * sizeof(float))
    CVM_ASSERT(pv, ((mnSize - 1) * nIncr + 1) * sizeof(float))
    SAXPY (&mnSize, &one, pv, &nIncr, mpD, &mnIncr);
}

template<>
CVM_API void __add<double>(double* mpD, int mnSize, int mnIncr, const double* pv, int nIncr)
{
    static const double one(1.);
    CVM_ASSERT(mpD, ((mnSize - 1) * mnIncr + 1) * sizeof(double))
    CVM_ASSERT(pv, ((mnSize - 1) * nIncr + 1) * sizeof(double))
    DAXPY (&mnSize, &one, pv, &nIncr, mpD, &mnIncr);
}

template<>
CVM_API void __subtract<float>(float* mpD, int mnSize, int mnIncr, const float* pv, int nIncr)
{
    static const float mone(-1.F);
    CVM_ASSERT(mpD, ((mnSize - 1) * mnIncr + 1) * sizeof(float))
    CVM_ASSERT(pv, ((mnSize - 1) * nIncr + 1) * sizeof(float))
    SAXPY (&mnSize, &mone, pv, &nIncr, mpD, &mnIncr);
}

template<>
CVM_API void __subtract<double>(double* mpD, int mnSize, int mnIncr, const double* pv, int nIncr)
{
    static const double mone(-1.F);
    CVM_ASSERT(mpD, ((mnSize - 1) * mnIncr + 1) * sizeof(double))
    CVM_ASSERT(pv, ((mnSize - 1) * nIncr + 1) * sizeof(double))
    DAXPY (&mnSize, &mone, pv, &nIncr, mpD, &mnIncr);
}

template<>
CVM_API void __add<std::complex<float> >(std::complex<float>* mpD, int mnSize, int mnIncr, const std::complex<float>* pv, int nIncr)
{
    static const std::complex<float> one(1.F);
    CVM_ASSERT(mpD, ((mnSize - 1) * mnIncr + 1) * sizeof(std::complex<float>))
    CVM_ASSERT(pv, ((mnSize - 1) * nIncr + 1) * sizeof(std::complex<float>))
    CAXPY (&mnSize, &one, pv, &nIncr, mpD, &mnIncr);
}

template<>
CVM_API void __add<std::complex<double> >(std::complex<double>* mpD, int mnSize, int mnIncr, const std::complex<double>* pv, int nIncr)
{
    static const std::complex<double> one(1.);
    CVM_ASSERT(mpD, ((mnSize - 1) * mnIncr + 1) * sizeof(std::complex<double>))
    CVM_ASSERT(pv, ((mnSize - 1) * nIncr + 1) * sizeof(std::complex<double>))
    ZAXPY (&mnSize, &one, pv, &nIncr, mpD, &mnIncr);
}

template<>
CVM_API void __subtract<std::complex<float> >(std::complex<float>* mpD, int mnSize, int mnIncr, const std::complex<float>* pv, int nIncr)
{
    static const std::complex<float> mone(-1.F);
    CVM_ASSERT(mpD, ((mnSize - 1) * mnIncr + 1) * sizeof(std::complex<float>))
    CVM_ASSERT(pv, ((mnSize - 1) * nIncr + 1) * sizeof(std::complex<float>))
    CAXPY (&mnSize, &mone, pv, &nIncr, mpD, &mnIncr);
}

template<>
CVM_API void __subtract<std::complex<double> >(std::complex<double>* mpD, int mnSize, int mnIncr, const std::complex<double>* pv, int nIncr)
{
    static const std::complex<double> mone(-1.);
    CVM_ASSERT(mpD, ((mnSize - 1) * mnIncr + 1) * sizeof(std::complex<double>))
    CVM_ASSERT(pv, ((mnSize - 1) * nIncr + 1) * sizeof(std::complex<double>))
    ZAXPY (&mnSize, &mone, pv, &nIncr, mpD, &mnIncr);
}

template<>
CVM_API void __scal<float, float>(float* mpD, int mnSize, int mnIncr, float dScal)
{
    CVM_ASSERT(mpD, ((mnSize - 1) * mnIncr + 1) * sizeof(float))
    SSCAL (&mnSize, &dScal, mpD, &mnIncr);
}
template<>
CVM_API void __scal<double, double>(double* mpD, int mnSize, int mnIncr, double dScal)
{
    CVM_ASSERT(mpD, ((mnSize - 1) * mnIncr + 1) * sizeof(double))
    DSCAL (&mnSize, &dScal, mpD, &mnIncr);
}

template<>
CVM_API void __scal<float, std::complex<float> >(std::complex<float>* mpD, int mnSize, int mnIncr, float dScal)
{
    CVM_ASSERT(mpD, ((mnSize - 1) * mnIncr + 1) * sizeof(std::complex<float>))
    CSSCAL (&mnSize, &dScal, mpD, &mnIncr);
}

template<>
CVM_API void __scal<double, std::complex<double> >(std::complex<double>* mpD, int mnSize, int mnIncr, double dScal)
{
    CVM_ASSERT(mpD, ((mnSize - 1) * mnIncr + 1) * sizeof(std::complex<double>))
    ZDSCAL (&mnSize, &dScal, mpD, &mnIncr);
}

template<>
CVM_API void __scal<std::complex<float>, std::complex<float> >(std::complex<float>* mpD, int mnSize, int mnIncr, std::complex<float> cScal)
{
    CVM_ASSERT(mpD, ((mnSize - 1) * mnIncr + 1) * sizeof(std::complex<float>))
    CSCAL (&mnSize, &cScal, mpD, &mnIncr);
}

template<>
CVM_API void __scal<std::complex<double>, std::complex<double> >(std::complex<double>* mpD, int mnSize, int mnIncr, std::complex<double> cScal)
{
    CVM_ASSERT(mpD, ((mnSize - 1) * mnIncr + 1) * sizeof(std::complex<double>))
    ZSCAL (&mnSize, &cScal, mpD, &mnIncr);
}

template<>
CVM_API void __copy2<float, std::complex<float> > (std::complex<float>* mpD, int mnSize, int mnIncr,
                                                   const float* pRe, const float* pIm, int nReIncr, int nImIncr)
{
    const int nIncr2 = mnIncr * 2;
    float* pR = __get_real_p<float>(mpD);
    float* pI = __get_imag_p<float>(mpD);

    if (pRe)
    {
        CVM_ASSERT(pRe, ((mnSize - 1) * nReIncr + 1) * sizeof(float))
        CVM_ASSERT(pR, ((mnSize - 1) * nIncr2 + 1) * sizeof(float))
        SCOPY (&mnSize, pRe, &nReIncr, pR, &nIncr2);
    }
    else
    {
        static const float zero(0.F);
        CVM_ASSERT(pR, ((mnSize - 1) * nIncr2 + 1) * sizeof(float))
        SSCAL (&mnSize, &zero, pR, &nIncr2);
    }

    if (pIm)
    {
        CVM_ASSERT(pIm, ((mnSize - 1) * nImIncr + 1) * sizeof(float))
        CVM_ASSERT(pI, ((mnSize - 1) * nIncr2 + 1) * sizeof(float))
        SCOPY (&mnSize, pIm, &nImIncr, pI, &nIncr2);
    }
    else
    {
        static const float zero(0.F);
        CVM_ASSERT(pI, ((mnSize - 1) * nIncr2 + 1) * sizeof(float))
        SSCAL (&mnSize, &zero, pI, &nIncr2);
    }
}

template<>
CVM_API void __copy2<double, std::complex<double> > (std::complex<double>* mpD, int mnSize, int mnIncr,
                                                     const double* pRe, const double* pIm, int nReIncr, int nImIncr)
{
    const int nIncr2 = mnIncr * 2;
    double* pR = __get_real_p<double>(mpD);
    double* pI = __get_imag_p<double>(mpD);

    if (pRe)
    {
        CVM_ASSERT(pRe, ((mnSize - 1) * nReIncr + 1) * sizeof(double))
        CVM_ASSERT(pR, ((mnSize - 1) * nIncr2 + 1) * sizeof(double))
        DCOPY (&mnSize, pRe, &nReIncr, pR, &nIncr2);
    }
    else
    {
        static const double zero(0.);
        CVM_ASSERT(pR, ((mnSize - 1) * nIncr2 + 1) * sizeof(double))
        DSCAL (&mnSize, &zero, pR, &nIncr2);
    }

    if (pIm)
    {
        CVM_ASSERT(pIm, ((mnSize - 1) * nImIncr + 1) * sizeof(double))
        CVM_ASSERT(pI, ((mnSize - 1) * nIncr2 + 1) * sizeof(double))
        DCOPY (&mnSize, pIm, &nImIncr, pI, &nIncr2);
    }
    else
    {
        static const double zero(0.);
        CVM_ASSERT(pI, ((mnSize - 1) * nIncr2 + 1) * sizeof(double))
        DSCAL (&mnSize, &zero, pI, &nIncr2);
    }
}

template<>
CVM_API void __copy_real<float, std::complex<float> > (std::complex<float>* mpD, int mnSize, int mnIncr, const float* pRe, int nReIncr)
{
    float* pD = __get_real_p<float>(mpD);
    if (pD != pRe)
    {
        const int nIncr2 = mnIncr * 2;
        CVM_ASSERT(pRe, ((mnSize - 1) * nReIncr + 1) * sizeof(float))
        CVM_ASSERT(pD,  ((mnSize - 1) * nIncr2 + 1) * sizeof(float))
        SCOPY (&mnSize, pRe, &nReIncr, pD, &nIncr2);
    }
}

template<>
CVM_API void __copy_real<double, std::complex<double> > (std::complex<double>* mpD, int mnSize, int mnIncr, const double* pRe, int nReIncr)
{
    double* pD = __get_real_p<double>(mpD);
    if (pD != pRe)
    {
        const int nIncr2 = mnIncr * 2;
        CVM_ASSERT(pRe, ((mnSize - 1) * nReIncr + 1) * sizeof(double))
        CVM_ASSERT(pD,  ((mnSize - 1) * nIncr2 + 1) * sizeof(double))
        DCOPY (&mnSize, pRe, &nReIncr, pD, &nIncr2);
    }
}

template<>
CVM_API void __copy_imag<float, std::complex<float> > (std::complex<float>* mpD, int mnSize, int mnIncr, const float* pIm, int nImIncr)
{
    float* pD = __get_imag_p<float>(mpD);
    if (pD != pIm)
    {
        const int nIncr2 = mnIncr * 2;
        CVM_ASSERT(pIm, ((mnSize - 1) * nImIncr + 1) * sizeof(float))
        CVM_ASSERT(pD,  ((mnSize - 1) * nIncr2 + 1) * sizeof(float))
        SCOPY (&mnSize, pIm, &nImIncr, pD, &nIncr2);
    }
}

template<>
CVM_API void __copy_imag<double, std::complex<double> > (std::complex<double>* mpD, int mnSize, int mnIncr, const double* pIm, int nImIncr)
{
    double* pD = __get_imag_p<double>(mpD);
    if (pD != pIm)
    {
        const int nIncr2 = mnIncr * 2;
        CVM_ASSERT(pIm, ((mnSize - 1) * nImIncr + 1) * sizeof(double))
        CVM_ASSERT(pD,  ((mnSize - 1) * nIncr2 + 1) * sizeof(double))
        DCOPY (&mnSize, pIm, &nImIncr, pD, &nIncr2);
    }
}

template<>
CVM_API void __conj<std::complex<float> > (std::complex<float>* mpD, int mnSize, int mnIncr)
{
    static const float mone(-1.F);
    const int nIncr2 = mnIncr * 2;
    SSCAL (&mnSize, &mone, __get_imag_p<float>(mpD), &nIncr2);
}

template<>
CVM_API void __conj<std::complex<double> > (std::complex<double>* mpD, int mnSize, int mnIncr)
{
    static const double mone(-1.);
    const int nIncr2 = mnIncr * 2;
    DSCAL (&mnSize, &mone, __get_imag_p<double>(mpD), &nIncr2);
}

template<>
CVM_API void __randomize<float> (float* mpD, int mnSize, int mnIncr, float dFrom, float dTo)
{
    const int nSize = mnSize * mnIncr;
    for (int i = 0; i < nSize; i += mnIncr)
    {
        mpD[i] = Randomizer<float>::get (dFrom, dTo);
    }
}

template<>
CVM_API void __randomize<double> (double* mpD, int mnSize, int mnIncr, double dFrom, double dTo)
{
    const int nSize = mnSize * mnIncr;
    for (int i = 0; i < nSize; i += mnIncr)
    {
        mpD[i] = Randomizer<double>::get (dFrom, dTo);
    }
}

template<>
CVM_API void __randomize_real<std::complex<float>, float> (std::complex<float>* mpD, int mnSize, int mnIncr,
                                                           float dFrom, float dTo)
{
    const int nIncr = 2 * mnIncr;
    const int nSize = mnSize * nIncr;
    float* pD = __get_real_p<float>(mpD);

    for (int i = 0; i < nSize; i += nIncr)
    {
        pD[i] = Randomizer<float>::get (dFrom, dTo);
    }
}

template<>
CVM_API void __randomize_real<std::complex<double>, double> (std::complex<double>* mpD, int mnSize, int mnIncr,
                                                             double dFrom, double dTo)
{
    const int nIncr = 2 * mnIncr;
    const int nSize = mnSize * nIncr;
    double* pD = __get_real_p<double>(mpD);

    for (int i = 0; i < nSize; i += nIncr)
    {
        pD[i] = Randomizer<double>::get (dFrom, dTo);
    }
}

template<>
CVM_API void __randomize_imag<std::complex<float>, float> (std::complex<float>* mpD, int mnSize, int mnIncr,
                                                           float dFrom, float dTo)
{
    const int nIncr = 2 * mnIncr;
    const int nSize = mnSize * nIncr;
    float* pD = __get_imag_p<float>(mpD);

    for (int i = 0; i < nSize; i += nIncr)
    {
        pD[i] = Randomizer<float>::get (dFrom, dTo);
    }
}

template<>
CVM_API void __randomize_imag<std::complex<double>, double> (std::complex<double>* mpD, int mnSize, int mnIncr,
                                                             double dFrom, double dTo)
{
    const int nIncr = 2 * mnIncr;
    const int nSize = mnSize * nIncr;
    double* pD = __get_imag_p<double>(mpD);

    for (int i = 0; i < nSize; i += nIncr)
    {
        pD[i] = Randomizer<double>::get (dFrom, dTo);
    }
}

CVM_NAMESPACE_END
