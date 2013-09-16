/* CVM Class Library */
/* Copyright (C), Sergei Nikolaev, 1992-2006, http://cvmlib.com */

#include "cvm.h"

#if defined (_MSC_VER)
#   pragma warning(disable:4786)
#endif

extern "C" {
    void __stdcall XERBLA (const char* szSubName,
    #ifdef CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES
                        const unsigned int,
    #endif
                        const int* pnParam) throw (cvm::cvmexception)
    {
        throw cvm::cvmexception (CVM_WRONGMKLARG2, *pnParam, szSubName);
    }
}

#if !defined (CVM_STATIC) && (defined (_MSC_VER) || defined (__WATCOMC__))
BOOL APIENTRY DllMain (HANDLE /*hModule*/,
                       DWORD  ul_reason_for_call,
                       LPVOID /*lpReserved*/)
{
    switch (ul_reason_for_call)
    {
        case DLL_PROCESS_ATTACH:
        case DLL_THREAD_ATTACH:
        case DLL_THREAD_DETACH:
        case DLL_PROCESS_DETACH:
            break;
    }
    return TRUE;
}
#endif


CVM_NAMESPACE_BEG

// multithreaded synchronizer
class CriticalSection {
#if defined (CVM_MT)
private:
    bool mbOK;

    #if defined (WIN32) || defined (_WIN32)
        ::CRITICAL_SECTION mCriticalSection; 
    #else                                                       // POSIX Threads library assumed
        pthread_mutex_t mMutex;
        pthread_mutexattr_t mMutexAttr;
    #endif
#endif

public:
    CriticalSection() 
#if defined (CVM_MT)
        : mbOK (false)
#endif
    {
#if defined (CVM_MT)
    #if defined (WIN32) || defined (_WIN32)
        if (!::InitializeCriticalSectionAndSpinCount (&mCriticalSection, 0x80000400))
        {
            ::InitializeCriticalSection (&mCriticalSection);
        }
        mbOK = true;
    #else
        if (pthread_mutexattr_init (&mMutexAttr) == 0 &&
            pthread_mutexattr_setpshared (&mMutexAttr, PTHREAD_PROCESS_PRIVATE) == 0 &&
            pthread_mutex_init (&mMutex, &mMutexAttr) == 0)
        {
            mbOK = true;
        }
    #endif
#endif
    }

    ~CriticalSection()
    {
#if defined (CVM_MT)
    #if defined (WIN32) || defined (_WIN32)
        if (mbOK)
        {
            ::DeleteCriticalSection (&mCriticalSection);
        }
    #else
        pthread_mutexattr_destroy (&mMutexAttr);
        pthread_mutex_destroy (&mMutex);
    #endif
#endif
    }

    void enter()
#if defined (CVM_MT)
    throw (cvmexception)
#endif
    {
#if defined (CVM_MT)
    #if defined (WIN32) || defined (_WIN32)
        if (!mbOK)
        {
            throw cvmexception (CVM_SEMAPHOREERROR);
        }
        ::EnterCriticalSection (&mCriticalSection); 
    #else
        if (!mbOK || pthread_mutex_lock (&mMutex) != 0)
        {
            throw cvmexception (CVM_SEMAPHOREERROR);
        }
    #endif
#endif
    }

    void leave()
#if defined (CVM_MT)
    throw (cvmexception)
#endif
    {
#if defined (CVM_MT)
    #if defined (WIN32) || defined (_WIN32)
        if (!mbOK)
        {
            throw cvmexception (CVM_SEMAPHOREERROR);
        }
        ::LeaveCriticalSection (&mCriticalSection);
    #else
        if (!mbOK || pthread_mutex_unlock (&mMutex) != 0)
        {
            throw cvmexception (CVM_SEMAPHOREERROR);
        }
    #endif
#endif
    }
};


// these species must be global allowing to avoid problems in multithreading environments
CriticalSection gCS;
MemoryPool gPool;


class Lock
{
public:
    Lock ();
    ~Lock ();
};

Lock::Lock ()
{
    gCS.enter();
}
Lock::~Lock ()
{
    gCS.leave();
}


// 5.5.2 - noved out of cvn.h
cvmexception::cvmexception (int nCause, ...)
    : mnCause (nCause)
{
    va_list argList;
    va_start (argList, nCause);
#if defined (CVM_VSNPRINTF_S_DEFINED)
    const int nLength = CVM_VSNPRINTF (mszMsg, sizeof(mszMsg), sizeof(mszMsg) - 1, _get_message(mnCause), argList);
#else
    const int nLength = CVM_VSNPRINTF (mszMsg, sizeof(mszMsg) - 1, _get_message(mnCause), argList);
#endif
    va_end (argList);
    if (nLength >= (int) sizeof(mszMsg))
    {
        mszMsg[sizeof(mszMsg) - 1] = '\0';
    }
}

cvmexception::cvmexception (const cvmexception& e)
    : std::exception(e), mnCause (e.mnCause)
{
#if defined (CVM_STRCPY_S_DEFINED)
    strcpy_s (mszMsg, sizeof(mszMsg), e.mszMsg);
#else
    strcpy (mszMsg, e.mszMsg);
#endif
}


#ifdef CVM_USE_POOL_MANAGER
CVM_API void _cvm_assert (const void* pvBlock, size_t nBytes)
{
    Lock l;
    gPool.Assert (pvBlock, nBytes);
}
#else
CVM_API void _cvm_assert (const void*, size_t)
{
}
#endif  // CVM_USE_POOL_MANAGER

CVM_API tbyte* _cvmMalloc (size_t nBytes) throw (cvmexception)
{
    return gPool.Malloc (nBytes);
}

CVM_API tbyte* _cvmAddRef (const tbyte* pD)
{
    return gPool.AddRef (pD);
}

CVM_API int _cvmFree (tbyte*& pD)
{
    return gPool.Free (pD);
}

CVM_API void cvmExit()
{
#ifdef CVM_USE_POOL_MANAGER
    gPool.Clear();
#endif
}

#ifdef CVM_USE_POOL_MANAGER
#define CVM_PAGE_SIZE (0x1000)
#define CVM_HEAP_SIZE ((int) 0x40000000)

size_t _up_value (size_t n)                                     // the least power of 2 multiplied by 2
{
    if (n < CVM_PAGE_SIZE)                                      // let small objects be in one page
    {
        n = CVM_PAGE_SIZE;
    }
    else //if (n < CVM_HEAP_SIZE)
    {
        int i = 0;
        while (n >> i) ++i;
        if (i && n & (1 << (i - 1)) - 1) ++i;                   // obey warning C4554 :)
        n = 1 << i;
    }
    return n;
}
#endif  // CVM_USE_POOL_MANAGER

MemoryPool::MemoryPool()
{
}

MemoryPool::~MemoryPool()
{
#ifdef CVM_USE_POOL_MANAGER
    Clear();
#endif
}

#ifdef CVM_USE_POOL_MANAGER
#ifdef __BORLANDC__
#    pragma warn -8091
#endif

void MemoryPool::Clear()
{
    std::for_each (mOutBlocks.rbegin(), mOutBlocks.rend(), MemoryPool::DeletePtr());
    mOutBlocks.clear();
}

#ifdef __BORLANDC__
#    pragma warn +8091
#endif
#endif  // CVM_USE_POOL_MANAGER


tbyte* MemoryPool::Malloc (size_t nBytes) throw (cvmexception)
{
    Lock lock;
#ifdef CVM_USE_POOL_MANAGER
    if (nBytes >= CVM_HEAP_SIZE) throw cvmexception (CVM_WRONGSIZE, nBytes);
    if (nBytes == 0) return static_cast<tbyte*>(NULL);

    tbyte* pB = mMemoryBlocks.GetFreeBlock (nBytes);

    if (pB == NULL)     // There is no suitable memory block. Let's create a new one.
    {
        const size_t nUpBytes = _up_value (nBytes);
        const size_t nRest    = nUpBytes - nBytes;

        try
        {
            pB = AllocatorInstance<tbyte>().allocate(nUpBytes, NULL);
        }
        catch (const std::bad_alloc&)
        {
        }
        if (pB == NULL)
        {
            throw (cvmexception (CVM_OUTOFMEMORY));
        }

        mOutBlocks.push_back (pB);
        mMemoryBlocks.AddPair (pB, nBytes, nRest);
    }
#else
    tbyte* pB = NULL;
    try
    {
        pB = AllocatorInstance<tbyte>().allocate(nBytes, NULL);
    }
    catch (const std::bad_alloc&)
    {
    }
    if (pB == NULL)
    {
        throw (cvmexception (CVM_OUTOFMEMORY));
    }

    mMemoryBlocks.AddNew (pB, nBytes);
#endif
    return pB;
}

tbyte* MemoryPool::AddRef (const tbyte* pD)
{
    Lock l;
    return mMemoryBlocks.AddRef (pD);
}

int MemoryPool::Free (tbyte*& pToFree) throw (cvmexception)
{
    Lock l;
    int nRefCounter = mMemoryBlocks.FreeBlock (pToFree);
    if (!nRefCounter)
    {
        pToFree = NULL;
    }
    return nRefCounter;
}

CVM_NAMESPACE_END
