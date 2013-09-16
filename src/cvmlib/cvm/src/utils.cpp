    /* CVM Class Library */
/* Copyright (C), Sergei Nikolaev, 1992-2006, http://cvmlib.com */

#include "cvm.h"
#include "blas.h"

#if defined (_MSC_VER)
#   pragma warning(disable:4786)
#   pragma warning(disable:4503)
#endif

CVM_NAMESPACE_BEG

char Chars::mchars[15] = {'T','N','U','L','P','Q','B','E','R','A','S','V','O','I','C'};

// global error messages holder
CVM_API ErrMessages& ErrMessages::ErrMessagesInstance()
{
    static ErrMessages _ErrMessages;
    return _ErrMessages;
}

CVM_API ErrMessages::ErrMessages()
    : msUnknown ("Unknown exception")
{
    mmMsg.insert (pair_Msg (CVM_OK,                      "All OK"));
    mmMsg.insert (pair_Msg (CVM_OUTOFMEMORY,             "Out of memory"));
    mmMsg.insert (pair_Msg (CVM_OUTOFRANGE,              "Index %d is out of range"));
    mmMsg.insert (pair_Msg (CVM_OUTOFRANGE1,             "First index %d is out of range"));
    mmMsg.insert (pair_Msg (CVM_OUTOFRANGE2,             "Second index %d is out of range"));
    mmMsg.insert (pair_Msg (CVM_WRONGSIZE,               "Wrong size %d"));
    mmMsg.insert (pair_Msg (CVM_SIZESMISMATCH,           "Sizes mismatch"));
    mmMsg.insert (pair_Msg (CVM_WRONGMKLARG,             "Wrong argument passed to BLAS or LAPACK subroutine"));
    mmMsg.insert (pair_Msg (CVM_WRONGMKLARG2,            "Wrong argument %d passed to BLAS or LAPACK subroutine %s"));
    mmMsg.insert (pair_Msg (CVM_SINGULARMATRIX,          "The diagonal element (or main minor) %d of the matrix is zero (or singular)"));
    mmMsg.insert (pair_Msg (CVM_NOTPOSITIVEDEFINITE,     "The leading minor of order %d (and hence the matrix itself) is not positive-definite"));
    mmMsg.insert (pair_Msg (CVM_WRONGCHOLESKYFACTOR,     "The diagonal element %d of the Cholesky factor (and hence the factor itself) is zero"));
    mmMsg.insert (pair_Msg (CVM_WRONGBUNCHKAUFMANFACTOR, "The diagonal element %d of the Bunch-Kaufman factor (and hence the factor itself) is zero"));
    mmMsg.insert (pair_Msg (CVM_NOTPOSITIVEDIAG,         "The diagonal element %d of the matrix is nonpositive. Equilibration failed"));
    mmMsg.insert (pair_Msg (CVM_CONVERGENCE_ERROR,       "Method failed to converge"));
    mmMsg.insert (pair_Msg (CVM_DIVISIONBYZERO,          "Division by zero"));

#if defined (WIN32) || defined (_WIN32)
    mmMsg.insert (pair_Msg (CVM_SEMAPHOREERROR,          "Critical Section access error"));
#else
    mmMsg.insert (pair_Msg (CVM_SEMAPHOREERROR,          "Semaphore access error"));
#endif

    mmMsg.insert (pair_Msg (CVM_READ_ONLY_ACCESS,        "Attempt to change a read-only element"));
    mmMsg.insert (pair_Msg (CVM_SUBMATRIXACCESSERROR,    "Attempt to access non-continuous submatrix as continuous array, see manual for details"));
    mmMsg.insert (pair_Msg (CVM_SUBMATRIXNOTAVAILABLE,   "Submatrix instantiation is not available for class \'%s\', see manual for details"));
    mmMsg.insert (pair_Msg (CVM_MATRIXNOTSYMMETRIC,      "The matrix passed doesn't appear to be symmetric"));
    mmMsg.insert (pair_Msg (CVM_MATRIXNOTHERMITIAN,      "The matrix passed doesn't appear to be hermitian"));
    mmMsg.insert (pair_Msg (CVM_BREAKS_HERMITIANITY,     "This operation could make the matrix non-hermitian. Use %s instead"));
    mmMsg.insert (pair_Msg (CVM_METHODNOTAVAILABLE,      "Method \'%s\' is not available for class \'%s\'. See programmer\'s reference for further details"));
    mmMsg.insert (pair_Msg (CVM_NOTIMPLEMENTED,          "Function is not implemented"));
}

CVM_API const std::string& ErrMessages::_get (int nException)
{
    citr_Msg i = mmMsg.size() > 0 ? mmMsg.find (nException) : mmMsg.end();
    return i == mmMsg.end() ? msUnknown : (*i).second;
}

CVM_API bool ErrMessages::_add (int nNewCause, const char* szNewMessage)
{
    bool bRes = true;
    itr_Msg i = mmMsg.find (nNewCause);
    if (i != mmMsg.end())
    {
        (*i).second = (*i).second + " | " + szNewMessage;       // Defenition is overlapped. This is not a good idea
        bRes = false;                                           // to do so, use CVM_THE_LAST_ERROR_CODE + 1 as an error code.
    }
    else
    {
        mmMsg.insert (pair_Msg (nNewCause, szNewMessage));      // new error definition
    }
    return bRes;
}


#ifdef CVM_USE_POOL_MANAGER

void MemoryBlocks::AddBlock (tbyte* pBlock, size_t nBytes, bool bOccupied)
{
    if (!bOccupied)                                             // Add freed block
    {
        itr_FreeIt j;
        itr_Blocks i = mBlocks.upper_bound (pBlock);
        itr_Blocks i_next = i;
                                                                // Is there upper neighboring memory block?
        if (i != mBlocks.end())
        {
            tbyte* pUpperBlock = (*i).first;
            j = mFreeIt.find (pUpperBlock);
            if (j != mFreeIt.end() && pBlock + nBytes == pUpperBlock)           // Yes. It's free and will be concatenated
            {
                nBytes += (*i).second.mnSize;
                ++i_next;
                mBlocks.erase (i);
                i = i_next;
                mFreeBs.erase ((*j).second);
                mFreeIt.erase (j);
            }
        }
                                                                // Is there lower neighboring memory block?
        if (i != mBlocks.begin() && mBlocks.size() > 0)
        {
            --i;
            tbyte* pLowerBlock = (*i).first;
            const size_t nLowerBytes = (*i).second.mnSize;
            j = mFreeIt.find (pLowerBlock);
            if (j != mFreeIt.end() && pLowerBlock + nLowerBytes == pBlock)      // Yes. It's free and will be concatenated
            {
                pBlock = pLowerBlock;
                nBytes += nLowerBytes;
                mBlocks.erase (i);
                mFreeBs.erase ((*j).second);
                mFreeIt.erase (j);
            }
        }
        mFreeIt[pBlock] = mFreeBs.insert (std::pair<int, tbyte*>(nBytes, pBlock));
    }

    mBlocks.insert (std::pair<tbyte*, BlockProperty>(pBlock, BlockProperty(nBytes, 1)));
}

int MemoryBlocks::FreeBlock (tbyte* pBlock)
{
    int nRefCounter = 0;
    itr_Blocks i = mBlocks.find (pBlock);

    if (i != mBlocks.end())
    {
        if (mFreeIt.find (pBlock) == mFreeIt.end())
        {
            nRefCounter = -- (*i).second.mnRefCount;
            if (nRefCounter <= 0)
            {
                const int nBytes = (*i).second.mnSize;
                mBlocks.erase (i);
                AddBlock (pBlock, nBytes, false);                               // return free block to the pool
            }
        }
#ifdef CVM_DEBUG
        else
        {
            assert (mFreeIt.find (pBlock) == mFreeIt.end());
        }
#endif
    }
    else
    {
        nRefCounter = -1;                                                       // foreign array.
    }

    return nRefCounter;
}

tbyte* MemoryBlocks::GetFreeBlock (size_t nBytes)
{
    tbyte* pBlock = NULL;
    if (mFreeBs.size() > 0)
    {
        // Is there a suitable memory block?
        itr_FreeBs i = mFreeBs.lower_bound (nBytes);

        // Yes. Let's use it.
        if (i != mFreeBs.end())
        {
            const size_t nRest = (*i).first - nBytes;
            pBlock = (*i).second;

            mFreeBs.erase (i);
            if (mFreeIt.size() > 0 && mFreeIt.find(pBlock) != mFreeIt.end())
            {
                mFreeIt.erase (pBlock);
            }
            if (mBlocks.size() > 0 && mBlocks.find(pBlock) != mBlocks.end())
            {
                mBlocks.erase (pBlock);
            }

            AddPair (pBlock, nBytes, nRest);
        }
    }
    return pBlock;
}

tbyte* MemoryBlocks::AddRef (const tbyte* pcBlock)
{
    tbyte* pBlock = const_cast<tbyte*>(pcBlock);
    itr_Blocks i = mBlocks.find (pBlock);
    if (i != mBlocks.end())
    {
        ++ (*i).second.mnRefCount;
    }
    else
    {
        // This is a foreign array. Leave it alone.
    }
    return pBlock;
}


#ifdef CVM_DEBUG
void MemoryBlocks::Assert (const void* pvBlock, size_t nBytes)
{
    tbyte* pBlock = (tbyte*) const_cast<void*>(pvBlock);
    itr_Blocks i = mBlocks.find (pBlock);
    if (i != mBlocks.end())
    {
        const size_t nSize = (*i).second.mnSize;
        assert (nSize >= nBytes);
    }
    else
    {
        tbyte* pB;
        size_t nB;
        itr_Blocks end = mBlocks.end();
        for (i = mBlocks.begin(); i != end; ++i)
        {
            pB = (*i).first;
            nB = (*i).second.mnSize;
            if (pBlock >= pB && pBlock < pB + nB)
            {
                tbyte* pBase = pB + nB;
                tbyte* pTest = pBlock + nBytes;
                assert (pTest <= pBase);
            }
        }
    }
}
#endif

void MemoryBlocks::AddPair (tbyte* pBlock, size_t nBytes, size_t nRest)
{
    AddBlock (pBlock, nBytes, true);                            // occupied block...
    if (nRest > 0)
    {
        AddBlock (pBlock + nBytes, nRest, false);               // ...and the rest free block
    }
}

#else // CVM_USE_POOL_MANAGER

tbyte* MemoryBlocks::AddRef (const tbyte* pcBlock)
{
    CVM_LONGEST_UINT nBlock = reinterpret_cast<CVM_LONGEST_UINT>(pcBlock);
    itr_Blocks i = mBlocks.find (nBlock);
    if (i != mBlocks.end())
    {
        ++(*i).second.mnRefCount;
    }
    return reinterpret_cast<tbyte*>(nBlock);
}

int MemoryBlocks::FreeBlock (tbyte* pBlock)
{
    int nRefCounter;
    itr_Blocks i = mBlocks.find (reinterpret_cast<CVM_LONGEST_UINT>(pBlock));
    if (i != mBlocks.end())
    {
        nRefCounter = -- (*i).second.mnRefCount;
        if (nRefCounter <= 0)
        {
            AllocatorInstance<tbyte>().deallocate(pBlock, (*i).second.mnSize);
            mBlocks.erase (i);
        }
    }
    else
    {
        nRefCounter = -1;
    }
    return nRefCounter;
}

#endif // !CVM_USE_POOL_MANAGER

template <>
CVM_API float _real<std::complex<float>, float> (const std::complex<float>& mT)
{
    return mT.real();
}
template <>
CVM_API double _real<std::complex<double>, double> (const std::complex<double>& mT)
{
    return mT.real();
}
template <>
CVM_API float _imag<std::complex<float>, float> (const std::complex<float>& mT)
{
    return mT.imag();
}
template <>
CVM_API double _imag<std::complex<double>, double> (const std::complex<double>& mT)
{
    return mT.imag();
}

template <>
CVM_API void __copy<float> (int nSize, const float* pFrom, int nFromIncr, float* pTo, int nToIncr)
{
    CVM_ASSERT(pFrom, ((nFromIncr) * (nSize - 1) + 1) * sizeof(float))
    CVM_ASSERT(pTo,   ((nToIncr) * (nSize - 1) + 1) * sizeof(float))
    SCOPY (&nSize, pFrom, &nFromIncr, pTo, &nToIncr);
}
template <>
CVM_API void __copy<double> (int nSize, const double* pFrom, int nFromIncr, double* pTo, int nToIncr)
{
    CVM_ASSERT(pFrom, ((nFromIncr) * (nSize - 1) + 1) * sizeof(double))
    CVM_ASSERT(pTo,   ((nToIncr) * (nSize - 1) + 1) * sizeof(double))
    DCOPY (&nSize, pFrom, &nFromIncr, pTo, &nToIncr);
}
template <>
CVM_API void __copy<std::complex<float> > (int nSize, const std::complex<float>* pFrom, int nFromIncr, std::complex<float>* pTo, int nToIncr)
{
    CVM_ASSERT(pFrom, ((nFromIncr) * (nSize - 1) + 1) * sizeof(std::complex<float>))
    CVM_ASSERT(pTo,   ((nToIncr) * (nSize - 1) + 1) * sizeof(std::complex<float>))
    CCOPY (&nSize, pFrom, &nFromIncr, pTo, &nToIncr);
}
template <>
CVM_API void __copy<std::complex<double> > (int nSize, const std::complex<double>* pFrom, int nFromIncr, std::complex<double>* pTo, int nToIncr)
{
    CVM_ASSERT(pFrom, ((nFromIncr) * (nSize - 1) + 1) * sizeof(std::complex<double>))
    CVM_ASSERT(pTo,   ((nToIncr) * (nSize - 1) + 1) * sizeof(std::complex<double>))
    ZCOPY (&nSize, pFrom, &nFromIncr, pTo, &nToIncr);
}
template <>
CVM_API void __copy<int> (int nSize, const int* pFrom, int nFromIncr, int* pTo, int nToIncr)
{
    CVM_ASSERT(pFrom, ((nFromIncr) * (nSize - 1) + 1) * sizeof(int))
    CVM_ASSERT(pTo,   ((nToIncr) * (nSize - 1) + 1) * sizeof(int))
    for (int i = 0; i < nSize; ++i)
    {
        pTo[i * nToIncr] = pFrom[i * nFromIncr];
    }
}

template <>
CVM_API void __swap<float> (int nSize, float* p1, int n1Incr, float* p2, int n2Incr)
{
    CVM_ASSERT(p1, ((n1Incr) * (nSize - 1) + 1) * sizeof(float))
    CVM_ASSERT(p2, ((n2Incr) * (nSize - 1) + 1) * sizeof(float))
    SSWAP (&nSize, p1, &n1Incr, p2, &n2Incr);
}
template <>
CVM_API void __swap<double> (int nSize, double* p1, int n1Incr, double* p2, int n2Incr)
{
    CVM_ASSERT(p1, ((n1Incr) * (nSize - 1) + 1) * sizeof(double))
    CVM_ASSERT(p2, ((n2Incr) * (nSize - 1) + 1) * sizeof(double))
    DSWAP (&nSize, p1, &n1Incr, p2, &n2Incr);
}
template <>
CVM_API void __swap<std::complex<float> > (int nSize, std::complex<float>* p1, int n1Incr, std::complex<float>* p2, int n2Incr)
{
    CVM_ASSERT(p1, ((n1Incr) * (nSize - 1) + 1) * sizeof(std::complex<float>))
    CVM_ASSERT(p2, ((n2Incr) * (nSize - 1) + 1) * sizeof(std::complex<float>))
    CSWAP (&nSize, p1, &n1Incr, p2, &n2Incr);
}
template <>
CVM_API void __swap<std::complex<double> > (int nSize, std::complex<double>* p1, int n1Incr, std::complex<double>* p2, int n2Incr)
{
    CVM_ASSERT(p1, ((n1Incr) * (nSize - 1) + 1) * sizeof(std::complex<double>))
    CVM_ASSERT(p2, ((n2Incr) * (nSize - 1) + 1) * sizeof(std::complex<double>))
    ZSWAP (&nSize, p1, &n1Incr, p2, &n2Incr);
}
template <>
CVM_API void __swap<int> (int nSize, int* p1, int n1Incr, int* p2, int n2Incr)
{
    int n;
    CVM_ASSERT(p1, (n1Incr * (nSize - 1) + 1) * sizeof(int))
    CVM_ASSERT(p2, (n2Incr * (nSize - 1) + 1) * sizeof(int))
    for (int i = 0; i < nSize; ++i)
    {
        n = p1[i * n1Incr];
        p1[i * n1Incr] = p2[i * n2Incr];
        p2[i * n2Incr] = n;
    }
}

template<>
CVM_API void
__low_up<basic_srmatrix<float> > (basic_srmatrix<float>& m, int* nPivots) throw (cvmexception)
{
    int nOutInfo = 0;
    SGETRF (m._pm(), m._pn(), m, m._pld(), nPivots, &nOutInfo);
    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
    if (nOutInfo > 0) throw cvmexception (CVM_SINGULARMATRIX, nOutInfo);
}

template<>
CVM_API void
__low_up<basic_srmatrix<double> > (basic_srmatrix<double>& m, int* nPivots) throw (cvmexception)
{
    int nOutInfo = 0;
    DGETRF (m._pm(), m._pn(), m, m._pld(), nPivots, &nOutInfo);
    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
    if (nOutInfo > 0) throw cvmexception (CVM_SINGULARMATRIX, nOutInfo);
}

template<>
CVM_API void
__low_up<basic_scmatrix<float, std::complex<float> > >
    (basic_scmatrix<float, std::complex<float> >& m, int* nPivots) throw (cvmexception)
{
    int nOutInfo = 0;
    CGETRF (m._pm(), m._pn(), m, m._pld(), nPivots, &nOutInfo);
    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
    if (nOutInfo > 0) throw cvmexception (CVM_SINGULARMATRIX, nOutInfo);
}

template<>
CVM_API void
__low_up<basic_scmatrix<double, std::complex<double> > >
    (basic_scmatrix<double, std::complex<double> >& m, int* nPivots) throw (cvmexception)
{
    int nOutInfo = 0;
    ZGETRF (m._pm(), m._pn(), m, m._pld(), nPivots, &nOutInfo);
    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
    if (nOutInfo > 0) throw cvmexception (CVM_SINGULARMATRIX, nOutInfo);
}

template<>
CVM_API void
__low_up<basic_srbmatrix<float> >
    (basic_srbmatrix<float>& m, int* nPivots) throw (cvmexception)
{
    int nOutInfo = 0;
    const int nKL = m.lsize();
    const int nKU = m.usize();
    m.resize_lu (nKL, nKL + nKU);
    SGBTRF (m._pm(), m._pn(), &nKL, &nKU, m, m._pld(), nPivots, &nOutInfo);
    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
    if (nOutInfo > 0) throw cvmexception (CVM_SINGULARMATRIX, nOutInfo);
}

template<>
CVM_API void
__low_up<basic_srbmatrix<double> >
    (basic_srbmatrix<double>& m, int* nPivots) throw (cvmexception)
{
    int nOutInfo = 0;
    const int nKL = m.lsize();
    const int nKU = m.usize();
    m.resize_lu (nKL, nKL + nKU);
    DGBTRF (m._pm(), m._pn(), &nKL, &nKU, m, m._pld(), nPivots, &nOutInfo);
    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
    if (nOutInfo > 0) throw cvmexception (CVM_SINGULARMATRIX, nOutInfo);
}

template<>
CVM_API void
__low_up<basic_scbmatrix<float, std::complex<float> > >
    (basic_scbmatrix<float, std::complex<float> >& m, int* nPivots) throw (cvmexception)
{
    int nOutInfo = 0;
    const int nKL = m.lsize();
    const int nKU = m.usize();
    m.resize_lu (nKL, nKL + nKU);
    CGBTRF (m._pm(), m._pn(), &nKL, &nKU, m, m._pld(), nPivots, &nOutInfo);
    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
    if (nOutInfo > 0) throw cvmexception (CVM_SINGULARMATRIX, nOutInfo);
}

template<>
CVM_API void
__low_up<basic_scbmatrix<double, std::complex<double> > >
    (basic_scbmatrix<double, std::complex<double> >& m, int* nPivots) throw (cvmexception)
{
    int nOutInfo = 0;
    const int nKL = m.lsize();
    const int nKU = m.usize();
    m.resize_lu (nKL, nKL + nKU);
    ZGBTRF (m._pm(), m._pn(), &nKL, &nKU, m, m._pld(), nPivots, &nOutInfo);
    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
    if (nOutInfo > 0) throw cvmexception (CVM_SINGULARMATRIX, nOutInfo);
}

template<>
CVM_API int
__cholesky<basic_srmatrix<float> >
    (basic_srmatrix<float>& m)                              // input is symmetric, output is triangular
{
    int nOutInfo = 0;
    SPOTRF (Chars::pU(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
            1,
#endif
            m._pm(), m, m._pld(), &nOutInfo);

    return nOutInfo;
}

template<>
CVM_API int
__cholesky<basic_srmatrix<double> >
    (basic_srmatrix<double>& m)                           // input is symmetric, output is triangular
{
    int nOutInfo = 0;
    DPOTRF (Chars::pU(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
            1,
#endif
            m._pm(), m, m._pld(), &nOutInfo);

    return nOutInfo;
}

template<>
CVM_API int
__cholesky<basic_scmatrix<float, std::complex<float> > >
    (basic_scmatrix<float, std::complex<float> >& m)                            // input is hermitian, output is triangular
{
    int nOutInfo = 0;
    CPOTRF (Chars::pU(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
            1,
#endif
            m._pm(), m, m._pld(), &nOutInfo);

    return nOutInfo;
}

template<>
CVM_API int
__cholesky<basic_scmatrix<double, std::complex<double> > >
    (basic_scmatrix<double, std::complex<double> >& m)                         // input is hermitian, output is triangular
{
    int nOutInfo = 0;
    ZPOTRF (Chars::pU(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
            1,
#endif
            m._pm(), m, m._pld(), &nOutInfo);

    return nOutInfo;
}

template<>
CVM_API void
__bunch_kaufman<basic_srmatrix<float> >
    (basic_srmatrix<float>& m, int* nPivots) throw (cvmexception)        // input is symmetric, output is square
{
    int nOutInfo = 0;
    const int lwork = m.msize() * 64;
    basic_rvector<float> work (lwork);
    SSYTRF (Chars::pU(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
            1,
#endif
            m._pm(), m, m._pld(), nPivots, work, &lwork, &nOutInfo);

    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
    if (nOutInfo > 0) throw cvmexception (CVM_SINGULARMATRIX, nOutInfo);
}

template<>
CVM_API void
__bunch_kaufman<basic_srmatrix<double> >
    (basic_srmatrix<double>& m, int* nPivots) throw (cvmexception)        // input is symmetric, output is square
{
    int nOutInfo = 0;
    const int lwork = m.msize() * 64;
    basic_rvector<double> work (lwork);
    DSYTRF (Chars::pU(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
            1,
#endif
            m._pm(), m, m._pld(), nPivots, work, &lwork, &nOutInfo);

    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
    if (nOutInfo > 0) throw cvmexception (CVM_SINGULARMATRIX, nOutInfo);
}

template<>
CVM_API void
__bunch_kaufman<basic_scmatrix<float, std::complex<float> > >
    (basic_scmatrix<float, std::complex<float> >& m, int* nPivots) throw (cvmexception)       // input is hermitian, output is square
{
    int nOutInfo = 0;
    const int lwork = m.msize() * 64;
    basic_cvector<float, std::complex<float> > work (lwork);
    CHETRF (Chars::pU(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
            1,
#endif
            m._pm(), m, m._pld(), nPivots, work, &lwork, &nOutInfo);

    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
    if (nOutInfo > 0) throw cvmexception (CVM_SINGULARMATRIX, nOutInfo);
}

template<>
CVM_API void
__bunch_kaufman<basic_scmatrix<double, std::complex<double> > >
    (basic_scmatrix<double, std::complex<double> >& m, int* nPivots) throw (cvmexception)       // input is hermitian, output is square
{
    int nOutInfo = 0;
    const int lwork = m.msize() * 64;
    basic_cvector<double, std::complex<double> > work (lwork);
    ZHETRF (Chars::pU(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
            1,
#endif
            m._pm(), m, m._pld(), nPivots, work, &lwork, &nOutInfo);

    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
    if (nOutInfo > 0) throw cvmexception (CVM_SINGULARMATRIX, nOutInfo);
}

template<>
CVM_API void
__ger<float, basic_rmatrix<float>, basic_rvector<float> >
    (basic_rmatrix<float>& m, 
    const basic_rvector<float>& vCol, 
    const basic_rvector<float>& vRow, 
    float dAlpha)
{
    CVM_ASSERT(m.get(), vCol.size() * vRow.size() * sizeof(float))
    SGER (vCol._psize(), vRow._psize(), &dAlpha, vCol, vCol._pincr(), vRow, vRow._pincr(), m, m._pld());
}

template<>
CVM_API void
__ger<double, basic_rmatrix<double>, basic_rvector<double> >
    (basic_rmatrix<double>& m,
    const basic_rvector<double>& vCol,
    const basic_rvector<double>& vRow,
    double dAlpha)
{
    CVM_ASSERT(m.get(), vCol.size() * vRow.size() * sizeof(double))
    DGER (vCol._psize(), vRow._psize(), &dAlpha, vCol, vCol._pincr(), vRow, vRow._pincr(), m, m._pld());
}

template<>
CVM_API void 
__geru<std::complex<float>, basic_cmatrix<float, std::complex<float> >, basic_cvector<float, std::complex<float> > >
    (basic_cmatrix<float, std::complex<float> >& m, 
    const basic_cvector<float, std::complex<float> >& vCol, 
    const basic_cvector<float, std::complex<float> >& vRow, 
    std::complex<float> cAlpha)
{
    CVM_ASSERT(m, vCol.size() * vRow.size() * sizeof(std::complex<float>))
    CGERU (vCol._psize(), vRow._psize(), &cAlpha, vCol, vCol._pincr(), vRow, vRow._pincr(), m, m._pld());
}

template<>
CVM_API void
__geru<std::complex<double>, basic_cmatrix<double, std::complex<double> >, basic_cvector<double, std::complex<double> > >
    (basic_cmatrix<double, std::complex<double> >& m, 
    const basic_cvector<double, std::complex<double> >& vCol, 
    const basic_cvector<double, std::complex<double> >& vRow, 
    std::complex<double> cAlpha)
{
    CVM_ASSERT(m, vCol.size() * vRow.size() * sizeof(std::complex<double>))
    ZGERU (vCol._psize(), vRow._psize(), &cAlpha, vCol, vCol._pincr(), vRow, vRow._pincr(), m, m._pld());
}

template<>
CVM_API void
__gerc<std::complex<float>, basic_cmatrix<float, std::complex<float> >, basic_cvector<float, std::complex<float> > >
    (basic_cmatrix<float, std::complex<float> >& m, 
    const basic_cvector<float, std::complex<float> >& vCol, 
    const basic_cvector<float, std::complex<float> >& vRow,
    std::complex<float> cAlpha)
{
    CVM_ASSERT(m, vCol.size() * vRow.size() * sizeof(std::complex<float>))
    CGERC (vCol._psize(), vRow._psize(), &cAlpha, vCol, vCol._pincr(), vRow, vRow._pincr(), m, m._pld());
}

template<>
CVM_API void __gerc<std::complex<double>, basic_cmatrix<double, std::complex<double> >, basic_cvector<double, std::complex<double> > >
    (basic_cmatrix<double, std::complex<double> >& m, 
    const basic_cvector<double, std::complex<double> >& vCol, 
    const basic_cvector<double, std::complex<double> >& vRow,
    std::complex<double> cAlpha)
{
    CVM_ASSERT(m, vCol.size() * vRow.size() * sizeof(std::complex<double>))
    ZGERC (vCol._psize(), vRow._psize(), &cAlpha, vCol, vCol._pincr(), vRow, vRow._pincr(), m, m._pld());
}

template<>
CVM_API void
__poequ<float, basic_srsmatrix<float>, basic_rvector<float> >
    (const basic_srsmatrix<float>& m,
     basic_rvector<float>& vScalings, 
     float& dCond,
     float& dMax)
{
    int nOutInfo = 0;
    SPOEQU (m._pm(), m, m._pld(), vScalings, &dCond, &dMax, &nOutInfo);
    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
    if (nOutInfo > 0) throw cvmexception (CVM_NOTPOSITIVEDIAG, nOutInfo);
}

template<>
CVM_API void
__poequ<double, basic_srsmatrix<double>, basic_rvector<double> >
    (const basic_srsmatrix<double>& m,
     basic_rvector<double>& vScalings, 
     double& dCond,
     double& dMax)
{
    int nOutInfo = 0;
    DPOEQU (m._pm(), m, m._pld(), vScalings, &dCond, &dMax, &nOutInfo);
    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
    if (nOutInfo > 0) throw cvmexception (CVM_NOTPOSITIVEDIAG, nOutInfo);
}

template<>
CVM_API void
__poequ<float, basic_schmatrix<float, std::complex<float> >, basic_rvector<float> >
    (const basic_schmatrix<float, std::complex<float> >& m,
     basic_rvector<float>& vScalings, 
     float& dCond,
     float& dMax)
{
    int nOutInfo = 0;
    CPOEQU (m._pm(), m, m._pld(), vScalings, &dCond, &dMax, &nOutInfo);
    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
    if (nOutInfo > 0) throw cvmexception (CVM_NOTPOSITIVEDIAG, nOutInfo);
}

template<>
CVM_API void
__poequ<double, basic_schmatrix<double, std::complex<double> >, basic_rvector<double> >
    (const basic_schmatrix<double, std::complex<double> >& m,
     basic_rvector<double>& vScalings, 
     double& dCond,
     double& dMax)
{
    int nOutInfo = 0;
    ZPOEQU (m._pm(), m, m._pld(), vScalings, &dCond, &dMax, &nOutInfo);
    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
    if (nOutInfo > 0) throw cvmexception (CVM_NOTPOSITIVEDIAG, nOutInfo);
}

CVM_NAMESPACE_END
