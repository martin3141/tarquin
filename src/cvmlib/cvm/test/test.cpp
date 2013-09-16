/* CVM Class Library regression test utility*/
/* Copyright (C), Sergei Nikolaev, 1992-2005, http://cvmlib.com */

#include "StdAfx.h"
#ifdef HAVE_CONFIG_H
#   include "../../../src/cvm.h"  // for kdevelop
#else
#   include "../src/cvm.h"
#endif

#if defined (_MSC_VER)
#   pragma warning(disable:4305)
#   if _MSC_VER < 1300
#       pragma warning(disable:4018)
#       pragma warning(disable:4018)
#       pragma warning(disable:4786)
#   endif
#endif

#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>

static time_t last_check = 0;

// multithreaded synchronizer
class CriticalSection {
private:
    bool mbOK;

#if defined (CVM_MT)
    #if defined (WIN32) || defined (_WIN32)
        ::CRITICAL_SECTION mCriticalSection; 
    #else                                                                       // POSIX Threads library assumed
pthread_mutex_t mMutex;
pthread_mutexattr_t mMutexAttr;
    #endif
#endif

public:
    CriticalSection () : mbOK (false)
    {
#if defined (CVM_MT)
    #if defined (WIN32) || defined (_WIN32)
        if (::InitializeCriticalSectionAndSpinCount (&mCriticalSection, 0x80000400))
        {
            mbOK = true;
        }
        else
        {
            ::InitializeCriticalSection (&mCriticalSection);
            mbOK = true;
        }
    #else

        if (pthread_mutexattr_init (&mMutexAttr) != 0)
        {
            std::cout << "FAILED TO pthread_mutexattr_init" << std::endl;
        }

        if (pthread_mutexattr_setpshared (&mMutexAttr, PTHREAD_PROCESS_PRIVATE) != 0)
        {
            std::cout << "FAILED TO pthread_mutexattr_setpshared" << std::endl;
        }

        if (pthread_mutex_init (&mMutex, &mMutexAttr) != 0)
        {
            std::cout << "FAILED TO pthread_mutex_init" << std::endl;
        }
    #endif
#endif
    }

    ~CriticalSection ()
    {
#if defined (CVM_MT)
    #if defined (WIN32) || defined (_WIN32)
        if (mbOK)
        {
            ::DeleteCriticalSection (&mCriticalSection);
        }
    #else
        if (pthread_mutexattr_destroy (&mMutexAttr) != 0)
        {
            std::cout << "FAILED TO pthread_mutexattr_destroy" << std::endl;
        }
        if (pthread_mutex_destroy (&mMutex) != 0)
        {
            std::cout << "FAILED TO pthread_mutex_destroy" << std::endl;
        }

    #endif
#endif
        mbOK = false;
    }

    void enter ()
    {
#if defined (CVM_MT)
    #if defined (WIN32) || defined (_WIN32)
        if (mbOK)
        {
            ::EnterCriticalSection (&mCriticalSection); 
        }
    #else

        if (pthread_mutex_lock (&mMutex) != 0)
        {
            std::cout << "FAILED TO pthread_mutex_lock" << std::endl;
        }
    #endif
#endif
    }

    void leave ()
    {
#if defined (CVM_MT)
    #if defined (WIN32) || defined (_WIN32)
        if (mbOK)
        {
            ::LeaveCriticalSection (&mCriticalSection);
        }
    #else
        if (pthread_mutex_unlock (&mMutex) != 0)
        {
            std::cout << "FAILED TO pthread_mutex_unlock" << std::endl;
        }
    #endif
#endif
    }
};


CriticalSection gCS;

class LockIt
{
public:
    LockIt ();
    ~LockIt ();
};

LockIt::LockIt ()
{
    gCS.enter();
}

LockIt::~LockIt ()
{
    gCS.leave();
}



#if defined (_MSC_VER)
    #pragma warning(disable:4700)
#endif

#ifndef CVM_NO_NAMESPACE
using namespace cvm;
#endif

#ifdef CVM_FLOAT
    #if defined (_MSC_VER)
        #define FILE_OUT "testout_win32_float.txt"
    #else
        #if defined (__BORLANDC__)
            #define FILE_OUT "testout_win32_borland_float.txt"
        #else
            #define FILE_OUT "testout_linux_float.txt"
        #endif
    #endif
#else
    #if defined (_MSC_VER)
        #define FILE_OUT "testout_win32.txt"
    #else
        #if defined (__BORLANDC__)
            #define FILE_OUT "testout_win32_borland.txt"
        #else
            #define FILE_OUT "testout_linux.txt"
        #endif
    #endif
#endif

class test_exception : public std::exception
{
    std::string m_s;
public:
    explicit test_exception (const char* szWhat) : m_s ("TEST FAILED: ") {m_s += szWhat;}
    virtual ~test_exception () throw () {}
    virtual const char* what () const throw ()
    {
        return m_s.c_str();
    }
};

tcomplex conj (tcomplex c)
{
    return tcomplex (c.real(), - c.imag());
}

treal mod (tcomplex c)
{
    return (treal) sqrt ((treal) c.real() * (treal) c.real() + (treal) c.imag() * (treal) c.imag());
}

void Report (const char* szMsg, std::ostream& os, int line) throw (test_exception)
{
    std::ostringstream oss;
    oss << szMsg << std::endl << "EXCEPTION ON LINE " << line << std::ends;
    std::string sMsg = oss.str ();
    os << sMsg << std::endl;
    throw test_exception (sMsg.c_str());
}

void CheckBool (bool b, bool bPattern, const char* szMsg, std::ostream& os, int line) throw (test_exception)
{
    LockIt l;
    time_t curr_time = time(NULL);
    time_t duration = curr_time - last_check;
    last_check = curr_time;
    if (duration > 0) os << "*** " << duration << " ";
    os << szMsg << std::endl << b << std::endl;
    if (b != bPattern)
    {
        Report (szMsg, os, line);
    }
}

void CheckBoolNoLock (bool b, bool bPattern, const char* szMsg, std::ostream& os, int line) throw (test_exception)
{
    time_t curr_time = time(NULL);
    time_t duration = curr_time - last_check;
    last_check = curr_time;
    if (duration > 0) os << "*** " << duration << " ";
    os << szMsg << std::endl << b << std::endl;
    if (b != bPattern)
    {
        Report (szMsg, os, line);
    }
}

void CheckInt (int v, int vPattern, const char* szMsg, std::ostream& os, int line) throw (test_exception)
{
    LockIt l;
    time_t curr_time = time(NULL);
    time_t duration = curr_time - last_check;
    last_check = curr_time;
    if (duration > 0) os << "*** " << duration << " ";
    os << szMsg << std::endl << v << std::endl;
    if (v != vPattern)
    {
        os << "Expected: " << vPattern << std::endl << "Returned " << v << std::endl;
        Report (szMsg, os, line);
    }
}

void CheckReal (treal v, treal vPattern, const char* szMsg, std::ostream& os, int line, treal rSp = cvmMachSp ()) throw (test_exception)
{
    LockIt l;
    time_t curr_time = time(NULL);
    time_t duration = curr_time - last_check;
    last_check = curr_time;
    if (duration > 0) os << "*** " << duration << " ";
    os << szMsg << std::endl << v << std::endl;
    const treal mp = (treal) fabs ((treal) vPattern);
    treal vn = v - vPattern;
    if (mp > (treal) 1.) vn /= mp;
    if (fabs (vn) > rSp)
    {
        os << "Expected: " << vPattern << std::endl << "Returned " << v << std::endl;
        Report (szMsg, os, line);
    }
}

void CheckComplex (tcomplex v, tcomplex vPattern, const char* szMsg, std::ostream& os, int line, treal rSp = cvmMachSp ()) throw (test_exception)
{
    LockIt l;
    time_t curr_time = time(NULL);
    time_t duration = curr_time - last_check;
    last_check = curr_time;
    if (duration > 0) os << "*** " << duration << " ";
    os << szMsg << std::endl << v << std::endl;
    const treal mp = mod (vPattern);
    tcomplex vn = v - vPattern;
    if (mp > (treal) 1.) vn /= mp;
    if (mod (vn) > rSp)
    {
        os << "Expected: " << vPattern << std::endl << "Returned " << v << std::endl;
        Report (szMsg, os, line);
    }
}

/*
void test_param (const srmatrix& m)
{
    treal r = m(1,1);
}
*/

std::ofstream os (FILE_OUT);    // share it among multiple executions

void cprint (const std::complex<treal>* p, int size)
{
    for (int i = 0; i < size; ++i)
    {
        std::cout << p[i] << " ";
    }
    std::cout << std::endl;
}

void print_solution (const srmatrix& a, const rvector& b)
{
    std::cout << a.solve(b);
}

srmatrix invert (const srmatrix& a)
{
    return a.inv();
}

treal& ret (rmatrix& m, const unsigned row, const unsigned col)
{
    return m.operator()(row+1, col+1).get();
}


#if !defined (WIN32) && !defined (_WIN32)
void *
#else
unsigned int __stdcall
#endif
TestBody (void*)
{
    {
        LockIt l;
        std::cout << "TESTS STARTED" << std::endl;
        if (os.bad ())
        {
            std::cout << "Error while creating file \"" FILE_OUT "\"" << std::endl;
        }
    }

    treal dPessimisticSp =
#ifdef CVM_FLOAT
        (treal) 5.e-4;
#else
        (treal) 1.e-12;
#endif

    treal dVeryPessimisticSp =
#ifdef CVM_FLOAT
        (treal) 5.e-3;
#else
        (treal) 1.e-7;
#endif

#ifdef CVM_FLOAT
        os.precision (7);
        std::cout.precision (5);
#else
        os.precision (5);
        std::cout.precision (17);
#endif
        os.setf (std::ios::scientific | std::ios::showpoint | std::ios::left);
        std::cout.setf (std::ios::scientific | std::ios::showpoint | std::ios::left); 

    try
    {
       
        last_check = time(NULL);

        
        {
            // Intel MKL 8.1 crash test
            int i, j;
            const int n = 1000;
            const int p = 100;

            rmatrix A(n, p);
            for (j = 0 ; j < p; ++j) {
               for (i = 0 ; i < n; ++i) {
                   A(i+1, j+1) = (treal)(i + j * p);
               }
            }
            rvector v(_cvm_min(n,p)) ;
            srmatrix mU(n) ;
            srmatrix mVH(p) ;

            v.svd(A, mU, mVH) ;

            rmatrix mv(n,p);
            mv.diag(0) = v;

#ifndef CVM_FLOAT    // this would be too naive to expect this precision from float-based algorithm for 1000x100 matrix, wouldn't it?
            CheckReal    ((A * ~mVH - mU * mv).norm(),  (treal) 0.,  "srmatrix svd", os, __LINE__, dVeryPessimisticSp);
            CheckReal    ((~A * mU - ~(mv * mVH)).norm(),  (treal) 0.,  "srmatrix svd", os, __LINE__, dVeryPessimisticSp);
#endif
        }

/*
        {
            treal d = -7.62199774163029530e-001;

            std::cout.precision (17);
            std::cout.setf (std::ios::scientific | std::ios::showpoint | std::ios::left); 

            {
                std::ofstream os ("out.txt");
                os.precision (17);
                os.setf (std::ios::scientific | std::ios::showpoint | std::ios::left);
                os << d;
            }
            {
                treal dcopy;
                std::ifstream is ("out.txt");
                //is.precision (17);
                is >> dcopy;

                std::cout << d << std::endl << dcopy << std::endl;
                assert(d == dcopy);
            }
        }
*/

        {
            LockIt l;

            // matrix in/out
            rmatrix m1 (3,4);
            rmatrix m2 (3,4);
            m1.randomize((treal) -5., (treal) 5.);

            {
                std::ofstream os ("_" FILE_OUT);
                os.precision (17);
                os.setf (std::ios::scientific | std::ios::showpoint | std::ios::left);
                os << m1;
            }
            {
                std::ifstream is ("_" FILE_OUT);
                is.precision(17);
                is >> m2;
            }

            CheckBoolNoLock (m1 == m2, true, "Matrix in/out", os, __LINE__);
        }

        {
            treal pi = (treal)3.1415926535897932384626433832795;
            treal ci = (treal)1.;
            tcomplex phase = exp(2*pi*ci); 

            rmatrix tmp(2,2);
            tmp(2,2) = 3.;
            cmatrix H(10,10);
            H(2,2) += phase*tmp(2,2);
            CheckComplex (H(2,2), tcomplex ((treal) 1.606474966574294e+03, treal (0.)),  
                          "tcomplex * type_proxy<treal>",  os, __LINE__, dPessimisticSp);
        }

        {
            rmatrix m(2,3);
            ret(m,0,1) = (treal)7.77;
            CheckReal (m(1,2), (treal) 7.77, "type_proxy::get",  os, __LINE__);
        }

        int i, j, l;

        treal    a1[100], a2[100], a3[100], a4[100];
        tcomplex c1[100], c2[100];
        for (i = 0; i < 100; i++)
        {
            a1[i] = (treal) (i + 1);
            a2[i] = (treal) (i + 1) / (treal) 10.;
            a3[i] = (treal) (i + 1) * (treal) 10.;
            a4[i] = (treal) (i + 1) / (treal) 100.;
            c1[i] = tcomplex(a1[i], a2[i]);
            c2[i] = tcomplex(a2[i], a4[i]);
        }

        const treal cs[] = {3., 0., 2., 1., -1., 2., 2., -1., 3., 0.,
                             0., 3., -1., -2., 0., -3., 5., 0.};
        const treal as[] = {1., 2., 1., 2., 5., -1., 1., -1., 20.};

// constructors
        rvector  rv;
        rvector  rv0  (10);
        rvector  rv1  (a1, 10);             // note: this constructor shares memory
        rvector  rv2  (a1, 10, 3);          // note: this constructor shares memory
        rvector  rv3  (11, (treal) 17.77);
        rvector  rv4  (rv2);                // note: this constructor copies memory, does not share

        rmatrix BIG_MAMA(100, 100);
        cmatrix BIG_MAMAC(100, 100);
        BIG_MAMA.randomize(0., 2.);
        BIG_MAMAC.randomize_real(0., 2.);
        BIG_MAMAC.randomize_imag(0., 2.);

        rmatrix  rm;
        rmatrix  rm0  (BIG_MAMA, 21, 34, 5, 6);
        rmatrix  rm1  (a1, 2, 3);           // note: this constructor shares memory
        rmatrix  rm2  (rm1);
        rmatrix  rm3  (rv2, true);          // column
        rmatrix  rm4  (rv2, false);         // row
        srmatrix srm;
        srmatrix srm0 (BIG_MAMA, 43, 47, 4);
        srmatrix srm1 (a1, 3);              // note: this constructor shares memory
        srmatrix srm2 (srm1);
        srmatrix srm30 (srm0);

        cvector  cv;
        cvector  cv0  (10);
        cvector  cv1  (a1, a2, 10);         // note: this constructor copies memory, does not share
        cvector  cv2  (a1, a2, 10, 3);      // note: this constructor copies memory, does not share
        cvector  cv3  (c1, 10, 3);          // note: this constructor shares memory
        cvector  cv4  (11, tcomplex ((treal) 1.3, (treal) 2.4));
        cvector  cv5  (cv3);
        cvector  cv6  (rv1, rv2);           // note: this constructor copies memory, does not share
        cvector  cv7  (a4, 10, true,  2);   // note: this constructor copies memory, does not share
        cvector  cv8  (a4, 10, false, 3);   // note: this constructor copies memory, does not share
        cvector  cv9  (rv2, true);
        cvector  cv10 (rv2, false);

        cmatrix  cm;
        cmatrix  cm0  (BIG_MAMAC, 68, 17, 5, 6);
        cmatrix  cm1  (a1, a2, 2, 3);       // note: this constructor copies memory, does not share
        cmatrix  cm2  (cm1);
        scmatrix scm;
        scmatrix scm0 (4);
        scmatrix scm1 (a1, a2, 3);          // note: this constructor copies memory, does not share
        scmatrix scm2 (scm1);

        srbmatrix srbm;
        srbmatrix srbm1 (a1, 4, 1, 2);

        scbmatrix scbm;
        scbmatrix scbm1 (c1, 4, 1, 2);

        srbmatrix srbm5 (5, 1, 2);
        scbmatrix scbm5 (5, 1, 2);

        srsmatrix srs1 (3);
        schmatrix sch1 (3);


// Array<TR,TC> derived features.
        CheckInt (rv   .size()  , 0,  "rv.size()",    os, __LINE__);
        CheckInt (rv0  .size()  , 10, "rv0.size()",   os, __LINE__);
        CheckInt (rv1  .size()  , 10, "rv1.size()",   os, __LINE__);
        CheckInt (rv2  .size()  , 10, "rv2.size()",   os, __LINE__);
        CheckInt (rv3  .size()  , 11, "rv3.size()",   os, __LINE__);
        CheckInt (rv4  .size()  , 10, "rv4.size()",   os, __LINE__);
        CheckInt (rm   .size()  , 0,  "rm.size()",    os, __LINE__);
        CheckInt (rm0  .size()  , 30, "rm0.size()",   os, __LINE__);
        CheckInt (rm1  .size()  , 6,  "rm1.size()",   os, __LINE__);
        CheckInt (rm2  .size()  , 6,  "rm2.size()",   os, __LINE__);
        CheckInt (rm3  .size()  , 10, "rm3.size()",   os, __LINE__);
        CheckInt (rm4  .size()  , 10, "rm4.size()",   os, __LINE__);
        CheckInt (srm  .size()  , 0,  "srm.size()",   os, __LINE__);
        CheckInt (srm0 .size()  , 16, "srm0.size()",  os, __LINE__);
        CheckInt (srm1 .size()  , 9,  "srm1.size()",  os, __LINE__);
        CheckInt (srm2 .size()  , 9,  "srm2.size()",  os, __LINE__);
        CheckInt (cv   .size()  , 0,  "cv.size()",    os, __LINE__);
        CheckInt (cv0  .size()  , 10, "cv0.size()",   os, __LINE__);
        CheckInt (cv1  .size()  , 10, "cv1.size()",   os, __LINE__);
        CheckInt (cv2  .size()  , 10, "cv2.size()",   os, __LINE__);
        CheckInt (cv3  .size()  , 10, "cv3.size()",   os, __LINE__);
        CheckInt (cv4  .size()  , 11, "cv4.size()",   os, __LINE__);
        CheckInt (cv5  .size()  , 10, "cv5.size()",   os, __LINE__);
        CheckInt (cv6  .size()  , 10, "cv6.size()",   os, __LINE__);
        CheckInt (cv7  .size()  , 10, "cv7.size()",   os, __LINE__);
        CheckInt (cv8  .size()  , 10, "cv8.size()",   os, __LINE__);
        CheckInt (cv9  .size()  , 10, "cv9.size()",   os, __LINE__);
        CheckInt (cv10 .size()  , 10, "cv10.size()",  os, __LINE__);
        CheckInt (cm   .size()  , 0,  "cm.size()",    os, __LINE__);
        CheckInt (cm0  .size()  , 30, "cm0.size()",   os, __LINE__);
        CheckInt (cm1  .size()  , 6,  "cm1.size()",   os, __LINE__);
        CheckInt (scm  .size()  , 0,  "scm.size()",   os, __LINE__);
        CheckInt (scm0 .size()  , 16, "scm0.size()",  os, __LINE__);
        CheckInt (scm1 .size()  , 9,  "scm1.size()",  os, __LINE__);
        CheckInt (scm2 .size()  , 9,  "scm2.size()",  os, __LINE__);
        CheckInt (srbm .size()  , 0,  "srbm.size()",  os, __LINE__);
        CheckInt (srbm1.size()  , 16, "srbm1.size()", os, __LINE__);
        CheckInt (scbm .size()  , 0,  "scbm.size()",  os, __LINE__);
        CheckInt (scbm1.size()  , 16, "scbm1.size()", os, __LINE__);
        CheckInt (srbm5.size()  , 20, "srbm5.size()",  os, __LINE__);   // 5 * (1 + 1 + 2)
        CheckInt (scbm5.size()  , 20, "scbm5.size()",  os, __LINE__);
        CheckInt (srs1.size()   , 9, "srs1.size()", os, __LINE__);
        CheckInt (sch1.size()   , 9, "sch1.size()", os, __LINE__);

        CheckInt (rv   .incr()  , 0, "rv.incr()",    os, __LINE__);
        CheckInt (rv0  .incr()  , 1, "rv0.incr()",   os, __LINE__);
        CheckInt (rv1  .incr()  , 1, "rv1.incr()",   os, __LINE__);
        CheckInt (rv2  .incr()  , 3, "rv2.incr()",   os, __LINE__);
        CheckInt (rv3  .incr()  , 1, "rv3.incr()",   os, __LINE__);
        CheckInt (rv4  .incr()  , 1, "rv4.incr()",   os, __LINE__);
        CheckInt (rm   .incr()  , 0, "rm.incr()",    os, __LINE__);
        CheckInt (rm0  .incr()  , 1, "rm0.incr()",   os, __LINE__);
        CheckInt (rm1  .incr()  , 1, "rm1.incr()",   os, __LINE__);
        CheckInt (rm2  .incr()  , 1, "rm2.incr()",   os, __LINE__);
        CheckInt (rm3  .incr()  , 1, "rm3.incr()",   os, __LINE__);
        CheckInt (rm4  .incr()  , 1, "rm4.incr()",   os, __LINE__);
        CheckInt (srm  .incr()  , 0, "srm.incr()",   os, __LINE__);
        CheckInt (srm0 .incr()  , 1, "srm0.incr()",  os, __LINE__);
        CheckInt (srm1 .incr()  , 1, "srm1.incr()",  os, __LINE__);
        CheckInt (srm2 .incr()  , 1, "srm2.incr()",  os, __LINE__);
        CheckInt (cv   .incr()  , 0, "cv.incr()",    os, __LINE__);
        CheckInt (cv0  .incr()  , 1, "cv0.incr()",   os, __LINE__);
        CheckInt (cv1  .incr()  , 1, "cv1.incr()",   os, __LINE__);
        CheckInt (cv2  .incr()  , 1, "cv2.incr()",   os, __LINE__);
        CheckInt (cv3  .incr()  , 3, "cv3.incr()",   os, __LINE__);
        CheckInt (cv4  .incr()  , 1, "cv4.incr()",   os, __LINE__);
        CheckInt (cv5  .incr()  , 1, "cv5.incr()",   os, __LINE__);
        CheckInt (cv6  .incr()  , 1, "cv6.incr()",   os, __LINE__);
        CheckInt (cv7  .incr()  , 1, "cv7.incr()",   os, __LINE__);
        CheckInt (cv8  .incr()  , 1, "cv8.incr()",   os, __LINE__);
        CheckInt (cv9  .incr()  , 1, "cv9.incr()",   os, __LINE__);
        CheckInt (cv10 .incr()  , 1, "cv10.incr()",  os, __LINE__);
        CheckInt (cm   .incr()  , 0, "cm.incr()",    os, __LINE__);
        CheckInt (cm0  .incr()  , 1, "cm0.incr()",   os, __LINE__);
        CheckInt (cm1  .incr()  , 1, "cm1.incr()",   os, __LINE__);
        CheckInt (cm2  .incr()  , 1, "cm2.incr()",   os, __LINE__);
        CheckInt (scm  .incr()  , 0, "scm.incr()",   os, __LINE__);
        CheckInt (scm0 .incr()  , 1, "scm0.incr()",  os, __LINE__);
        CheckInt (scm1 .incr()  , 1, "scm1.incr()",  os, __LINE__);
        CheckInt (scm2 .incr()  , 1, "scm2.incr()",  os, __LINE__);
        CheckInt (srbm .incr()  , 0, "srbm.incr()",  os, __LINE__);
        CheckInt (srbm1.incr()  , 1, "srbm1.incr()", os, __LINE__);
        CheckInt (scbm .incr()  , 0, "scbm.incr()",  os, __LINE__);
        CheckInt (scbm1.incr()  , 1, "scbm1.incr()", os, __LINE__);
        CheckInt (srs1.incr()  , 1, "srs1.incr()", os, __LINE__);
        CheckInt (sch1.incr()  , 1, "sch1.incr()", os, __LINE__);

        cmatrix cm1r (rm1);
        cmatrix cm1i (rm1, false);
        CheckComplex (cm1r[1][2], tcomplex(rm1(1,2),(treal)0.), "cmatrix(rmatrix)",   os, __LINE__);
        CheckComplex (cm1i[1][2], tcomplex((treal)0.,rm1(1,2)), "cmatrix(rmatrix)",   os, __LINE__);

        rm2.resize(4, 4);
        cm2.resize(4, 4);
        srmatrix srm3 (rm2);
        scmatrix scm3 (cm2);

        CheckInt (rm2  .size()  , 16, "rm2.size()",   os, __LINE__);
        CheckInt (cm2  .size()  , 16, "cm2.size()",   os, __LINE__);
        CheckInt (rm2  .incr()  , 1,  "rm2.incr()",   os, __LINE__);
        CheckInt (cm2  .incr()  , 1,  "cm2.incr()",   os, __LINE__);
        CheckInt (srm3 .size()  , 16, "srm3.size()",  os, __LINE__);
        CheckInt (scm3 .size()  , 16, "scm3.size()",  os, __LINE__);
        CheckInt (srm3 .incr()  , 1,  "srm3.incr()",  os, __LINE__);
        CheckInt (scm3 .incr()  , 1,  "scm3.incr()",  os, __LINE__);

        scmatrix scm4 (srm3, true);
        scmatrix scm5 (srm3, false);

        CheckInt (scm4  .size() , 16, "scm4.size()",  os, __LINE__);
        CheckInt (scm5  .size() , 16, "scm5.size()",  os, __LINE__);
        CheckInt (scm4  .incr() , 1 , "scm4.incr()",  os, __LINE__);
        CheckInt (scm5  .incr() , 1 , "scm5.incr()",  os, __LINE__);

        scmatrix scm6 (srm3, srm0);
        CheckInt (scm6  .size() , 16, "scm6.size()",  os, __LINE__);
        CheckInt (scm6  .incr() , 1 , "scm6.incr()",  os, __LINE__);

        srmatrix srm4 (srbm1);
        CheckInt (srm4 .size()  , 16, "srm4.size()",  os, __LINE__);
        CheckInt (srm4 .incr()  , 1,  "srm4.incr()",  os, __LINE__);

        scmatrix scm8 (scbm1);
        CheckInt (scm8 .size()  , 16, "scm8.size()",  os, __LINE__);
        CheckInt (scm8 .incr()  , 1,  "scm8.incr()",  os, __LINE__);

        srmatrix srm5 (rv2);
        CheckInt (srm5 .size()  , 100,"srm5.size()",  os, __LINE__);
        CheckInt (srm5 .incr()  , 1,  "srm5.incr()",  os, __LINE__);

        scmatrix scm7 (cv2);
        CheckInt (scm7 .size()  , 100,"srm7.size()",  os, __LINE__);
        CheckInt (scm7 .incr()  , 1,  "srm7.incr()",  os, __LINE__);

// Indexing and assignments
        CheckReal (rv2[1],  1.,  "rv2[1]",   os, __LINE__);
        CheckReal (rv2[10], 28., "rv2[10]",  os, __LINE__);

        treal r1 = (treal) -1.92;
        a1[3] = r1;
        CheckReal (rv1[4],     r1,  "rv2[4]",      os, __LINE__);     // memoty sharing check
        CheckReal (rv2[2],     r1,  "rv2[2]",      os, __LINE__);
        CheckReal (rm1[2][2],  r1,  "rm1[2][2]",   os, __LINE__);
        CheckReal (rm1(2, 2),  r1,  "rm1(2, 2)",   os, __LINE__);
        CheckReal (srm1[1][2], r1,  "srm1[1][2]",  os, __LINE__);
        CheckReal (srm1(1, 2), r1,  "srm1(1, 2)",  os, __LINE__);

        CheckComplex (cm1(2, 2), tcomplex ((treal) 4., (treal) 0.4),  "cm1(2, 2)",  os, __LINE__);

        treal* pr1 = srm30;
        srm30(4,4) = r1;
        CheckReal (srm30(4,4),  r1,  "srm30(4,4)",      os, __LINE__);     // memoty sharing check
        CheckReal (pr1[15],    r1,  "pr1[15]",        os, __LINE__);     // memoty sharing check

        treal* pr2 = srbm1;
        srbm1(1,2) = r1;
        CheckReal (srbm1(1,2), r1,  "srbm1(1,2)",     os, __LINE__);     // memoty sharing check
        CheckReal (pr2[5],     r1,  "pr2[5]",         os, __LINE__);     // memoty sharing check

        tcomplex cr1 = tcomplex ((treal) 1.07, (treal) -0.179);
        std::complex<treal>* pc2 = scbm1;
        scbm1(1,2) = cr1;
        CheckComplex (scbm1(1,2), cr1,  "scbm1(1,2)",     os, __LINE__);     // memoty sharing check
        CheckComplex (pc2[5],     cr1,  "pc2[5]",         os, __LINE__);     // memoty sharing check

        srs1.assign(as);
        sch1.assign((tcomplex*)cs);

        CheckReal (srs1[1][2], as[3], "srs1[1][2]",   os, __LINE__);
        CheckReal (srs1(3,2),  as[5], "srs1(3,2)",   os, __LINE__);
        CheckReal (srs1(3)[2], as[7], "srs1(3)[2]",  os, __LINE__);

        CheckComplex (sch1[1][2], tcomplex(cs[6],cs[7]), "sch1[1][2]",   os, __LINE__);
        CheckComplex (sch1(3,2),  tcomplex(cs[10],cs[11]), "sch1(3,2)",   os, __LINE__);
        CheckComplex (sch1(3)[2], tcomplex(cs[14],cs[15]), "sch1(3)[2]",  os, __LINE__);


// Array<TR,TC> derived features -  continued
        rv << rv1.normalize();
        CheckReal (rv(7), rv1[7], "rvector << rvector",    os, __LINE__);

        treal r2 = (treal) 0.;
        for (i = 0; i < 10; i++)
        {
            r2 += a1[i] * a1[i];
        }

        CheckReal (r2,         (treal) 1.,  "normalize",       os, __LINE__, dPessimisticSp);
        CheckReal (rv.norm(),  (treal) 1.,  "rv.norm()",       os, __LINE__);

        CheckInt (rv.indofmax () ,    10,  "rv.indofmax ()",  os, __LINE__);
        CheckInt (rv.indofmin () ,    1 ,  "rv.indofmin ()",  os, __LINE__);

        r1 = rv[10];
        CheckReal (rv.norminf (), r1,  "rv.norminf ()",   os, __LINE__);

        rv1.sum (rv1, rv);
        CheckReal (rv1[10], r1 + r1,  "sum",   os, __LINE__);
        rv1.diff (rv1, rv);
        CheckReal (rv1[10], r1,       "diff",   os, __LINE__);

        rv1 += rv;
        CheckReal (rv1[10], r1 + r1,  "+=, rvector",   os, __LINE__);
        rv1 -= rv;
        CheckReal (rv1[10], r1,       "-=, rvector",   os, __LINE__);
        rv1 += rv1;
        CheckReal (rv1[10], r1 + r1,  "+=, rvector self",   os, __LINE__);

        cv << cv1;
        CheckComplex (cv(7), cv1[7], "cvector << cvector",    os, __LINE__);

        cr1 = cv1[10];
        cv1 += cv;
        CheckComplex (cv1[10], cr1 + cr1,  "+=, cvector",   os, __LINE__);
        cv1 -= cv;
        CheckComplex (cv1[10], cr1,        "-=, cvector",   os, __LINE__);
        cv1 += cv1;
        CheckComplex (cv1[10], cr1 + cr1,  "+=, cvector self",   os, __LINE__);

        rm << rm1;
        CheckReal (rm(2, 1), rm1[2][1], "rmatrix << rmatrix",    os, __LINE__);

        r1 = rm(2, 2);
        rm += rm1;
        CheckReal (rm(2, 2), r1 + r1, "+=, rmatrix",   os, __LINE__);
        rm -= rm1;
        CheckReal (rm(2, 2), r1,      "-=, rmatrix",   os, __LINE__);
        rm += rm;
        CheckReal (rm(2, 2), r1 + r1, "+=, rmatrix, self",   os, __LINE__);

        cm << cm1;
        CheckComplex (cm(2, 1), cm1[2][1], "cmatrix << cmatrix",    os, __LINE__);

        cr1 = cm(2, 2);
        cm += cm1;
        CheckComplex (cm(2, 2), cr1 + cr1, "+=, cmatrix",   os, __LINE__);
        cm -= cm1;
        CheckComplex (cm(2, 2), cr1,       "-=, cmatrix",   os, __LINE__);
        cm += cm;
        CheckComplex (cm(2, 2), cr1 + cr1, "+=, cmatrix, self",   os, __LINE__);

        srm << srm1;
        r1 = srm(2, 2);
        srm += srm1;
        CheckReal (srm(2, 2), r1 + r1, "+=, srmatrix",   os, __LINE__);
        srm -= srm1;
        CheckReal (srm(2, 2), r1,      "-=, srmatrix",   os, __LINE__);
        srm += srm;
        CheckReal (srm(2, 2), r1 + r1, "+=, srmatrix, self",   os, __LINE__);

        scm << scm1;
        CheckComplex (scm(2, 1), scm1[2][1], "scmatrix << scmatrix",    os, __LINE__);

        cr1 = scm(2, 2);
        scm += scm1;
        CheckComplex (scm(2, 2), cr1 + cr1, "+=, scmatrix",   os, __LINE__);
        scm -= scm1;
        CheckComplex (scm(2, 2), cr1,       "-=, scmatrix",   os, __LINE__);
        scm += scm;
        CheckComplex (scm(2, 2), cr1 + cr1, "+=, scmatrix, self",   os, __LINE__);

        srbm1.set((treal)1.14);
        srbm << srbm1;
        srbm.set( (treal)-.684);
        r1 = srbm(2, 1);
        r2 = srbm1(2, 1);
        srbm += srbm1;
        CheckReal (srbm(2, 1), r1 + r2, "+=, srbmatrix",   os, __LINE__);
        srbm -= srbm1;
        CheckReal (srbm(2, 1), r1, "-=, srbmatrix",   os, __LINE__);

        tcomplex cr2 = tcomplex ((treal) 1.03, (treal) -0.79);
        scbm1.set(cr2);
        scbm << scbm1;
        scbm.set(cr1);
        cr1 = scbm(2, 1);
        cr2 = scbm1(2, 1);
        scbm += scbm1;

        CheckComplex (scbm(2, 1), cr1 + cr2, "+=, scbmatrix",   os, __LINE__);
        scbm -= scbm1;
        CheckComplex (scbm(2, 1), cr1, "-=, scbmatrix",   os, __LINE__);


        srsmatrix srs2;
        schmatrix sch2;
        srs2 << srs1;
        sch2 << sch1;
        CheckReal ((srs2 - srs1).norm(),  (treal) 0. ,  "srsmatrix <<",  os, __LINE__);
        CheckReal ((sch2 - sch1).norm(),  (treal) 0. ,  "schmatrix <<",  os, __LINE__);

        srsmatrix srs2sub (srs2, 2, 2);
        CheckReal    (srs2sub(1,2), srs2(2,3), "srsmatrix submatrix ctr",   os, __LINE__);
        schmatrix sch2sub (sch2, 2, 2);
        CheckComplex (sch2sub(1,2), sch2(2,3), "schmatrix submatrix ctr",   os, __LINE__);


        r2 = (treal) 1.13;
        cr2 = tcomplex ((treal) 1.03, (treal) -0.79);

        treal rs1 = srs2(2,3);
        tcomplex cs1 = sch2(2,3);
        srs2 *= r2;
        sch2 *= r2;
        CheckReal (srs2(2,3),  rs1 * r2 ,  "srsmatrix *= TR",  os, __LINE__);
        CheckComplex (sch2(2,3),  cs1 * r2 ,  "schmatrix *= TR",  os, __LINE__);
        srs2 /= r2;
        sch2 /= r2;
        CheckReal (srs2(2,3),  rs1,  "srsmatrix /= TR",  os, __LINE__);
        CheckComplex (sch2(2,3),  cs1,  "schmatrix /= TR",  os, __LINE__, dPessimisticSp);
        CheckComplex ((sch2 * cr2)(2,3),  cs1 * cr2,  "schmatrix * TC",  os, __LINE__, dPessimisticSp);

        rvector vrs1(3);
        vrs1.randomize(3., 7.);
        CheckReal ((srs1 * vrs1 - srmatrix(srs1) * vrs1).norm(), (treal) 0.,  "srsmatrix * rvector",  os, __LINE__, dPessimisticSp);

        cvector vch1(3);
        vch1.randomize_real(3., 7.);
        vch1.randomize_imag(-3., 7.);
        CheckReal ((sch1 * vch1 - scmatrix(sch1) * vch1).norm(), (treal) 0.,  "schmatrix * cvector",  os, __LINE__, dPessimisticSp);

        r1 = rv1(9);
        rv1 *= r2;
        CheckReal (rv1(9),     r1 * r2,  "*=, rvector",     os, __LINE__);
        r1 = rm3(7,1);
        rm3 *= r2;
        CheckReal (rm3(7,1),   r1 * r2,  "*=, rmatrix",     os, __LINE__);
        r1 = srm4(1,2);
        srm4 *= r2;
        CheckReal (srm4(1,2),  r1 * r2,  "*=, srmatrix",    os, __LINE__);
        r1 = srbm1(1,2);
        srbm1 *= r2;
        CheckReal (srbm1(1,2), r1 * r2,  "*=, srbmatrix",   os, __LINE__);
        cr1 = scbm1(1,2);
        scbm1 *= cr2;
        CheckComplex (scbm1(1,2), cr1 * cr2,  "*=, scbmatrix",   os, __LINE__);
        r1 = rv1(9);
        rv1 /= r2;
        CheckReal (rv1(9),     r1 / r2,  "/=, rvector",     os, __LINE__);
        r1 = rm3(7,1);

        rm3 /= r2;
        CheckReal (rm3(7,1),   r1 / r2,  "/=, rmatrix",     os, __LINE__);
        r1 = srm4(1,2);
        srm4 /= r2;
        CheckReal (srm4(1,2),  r1 / r2,  "/=, srmatrix",    os, __LINE__);
        r1 = srbm1(1,2);
        srbm1 /= r2;
        CheckReal (srbm1(1,2), r1 / r2,  "/=, srbmatrix",   os, __LINE__);
        cr1 = scbm1(2,1);
        scbm1 /= cr2;
        CheckComplex (scbm1(2,1), cr1 / cr2,  "/=, scbmatrix",   os, __LINE__);


        cr1 = cv1(9);
        cv1 *= r2;
        CheckComplex (cv1(9),       cr1 * r2,  "cvector *= treal",    os, __LINE__);
        cr1 = cm1(2, 2);
        cm1 *= r2;
        CheckComplex (cm1(2, 2),    cr1 * r2,  "cmatrix *= treal",    os, __LINE__);
        cr1 = scm1(2, 2);
        scm1 *= r2;
        CheckComplex (scm1(2, 2),   cr1 * r2,  "scmatrix *= treal",   os, __LINE__);
        cr1 = scbm1(1, 2);
        scbm1 *= r2;
        CheckComplex (scbm1(1, 2),  cr1 * r2,  "scbmatrix *= treal",   os, __LINE__);
        cr1 = cv1(9);
        cv1 /= r2;
        CheckComplex (cv1(9),       cr1 / r2,  "cvector /= treal",    os, __LINE__);
        cr1 = cm1(2, 2);
        cm1 /= r2;
        CheckComplex (cm1(2, 2),    cr1 / r2,  "cmatrix /= treal",    os, __LINE__);
        cr1 = scm1(2, 2);
        scm1 /= r2;
        CheckComplex (scm1(2, 2),   cr1 / r2,  "scmatrix /= treal",   os, __LINE__);
        cr1 = scbm1(2, 1);
        scbm1 /= r2;
        CheckComplex (scbm1(2, 1),  cr1 / r2,  "scbmatrix /= treal",   os, __LINE__);


        cr2 = tcomplex ((treal) 1.03, (treal) -0.79);

        cr1 = cv1(9);
        cv1 *= cr2;
        CheckComplex (cv1(9),       cr1 * cr2, "cvector *= tcomplex",  os, __LINE__);
        cr1 = cm1(2, 2);
        cm1 *= cr2;
        CheckComplex (cm1(2, 2),    cr1 * cr2, "cmatrix *= tcomplex",  os, __LINE__);
        cr1 = scm1(2, 2);
        scm1 *= cr2;
        CheckComplex (scm1(2, 2),   cr1 * cr2, "scmatrix *= tcomplex", os, __LINE__);
        cr1 = scbm1(2, 1);
        scbm1 *= cr2;
        CheckComplex (scbm1(2, 1),  cr1 * cr2, "scbmatrix *= tcomplex", os, __LINE__);
        cr1 = cv1(9);
        cv1 /= cr2;
        CheckComplex (cv1(9),       cr1 / cr2, "cvector /= tcomplex",  os, __LINE__);
        cr1 = cm1(2, 2);
        cm1 /= cr2;
        CheckComplex (cm1(2, 2),    cr1 / cr2, "cmatrix /= tcomplex",  os, __LINE__);
        cr1 = scm1(2, 2);
        scm1 /= cr2;
        CheckComplex (scm1(2, 2),   cr1 / cr2, "scmatrix /= tcomplex", os, __LINE__);
        cr1 = scbm1(1, 2);
        scbm1 /= cr2;
        CheckComplex (scbm1(2, 2),  cr1 / cr2, "scbmatrix /= tcomplex", os, __LINE__);

        srbm << srbm1;
        CheckReal (srbm(2, 3), srbm1(2, 3), "srbmatrix << srbmatrix",    os, __LINE__);
        CheckReal (srbm(1, 4), srbm1(1, 4), "srbmatrix << srbmatrix",    os, __LINE__);
        scbm << scbm1;
        CheckComplex (scbm(2, 3), scbm1(2, 3), "scbmatrix << scbmatrix",    os, __LINE__);
        CheckComplex (scbm(1, 4), scbm1(1, 4), "scbmatrix << scbmatrix",    os, __LINE__);


        srs2.set ((treal) 2.3);
        CheckReal (srs2(1,3), (treal) 2.3,  "srsmatrix.set",  os, __LINE__, dPessimisticSp);
        CheckReal (srs2(3,2), (treal) 2.3,  "srsmatrix.set",  os, __LINE__, dPessimisticSp);

        sch2.set_real((treal) 2.3);
        CheckReal (sch2.real()(1,3), (treal) 2.3,  "schmatrix.set_real",  os, __LINE__, dPessimisticSp);
        CheckReal (sch2.real()(3,2), (treal) 2.3,  "schmatrix.set_real",  os, __LINE__, dPessimisticSp);



        r1 = (treal) -0.127;
        rv.set(r1);
        CheckReal (rv[1],     r1,  "rvector = treal",      os, __LINE__);
        rm.set(r1);
        CheckReal (rm[1][2],  r1,  "rmatrix = treal",      os, __LINE__);
        srm.set(r1);
        CheckReal (srm(2, 2), r1,  "srmatrix = treal",     os, __LINE__);
        srbm.set(r1);
        CheckReal (srbm(2, 3), r1, "srbmatrix = treal",    os, __LINE__);
        CheckReal (srbm(1, 4), 0,  "srbmatrix = treal, 0", os, __LINE__);

        cr2 = tcomplex ((treal) 1.3, (treal) -0.9);        
        cv.set(cr1);
        CheckComplex (cv[7],     cr1,  "cvector = tcomplex",   os, __LINE__);
        cm.set(cr1);
        CheckComplex (cm[2][3],  cr1,  "cmatrix = tcomplex",   os, __LINE__);
        scm.set(cr1);
        CheckComplex (scm(3, 2), cr1,  "scmatrix = tcomplex",  os, __LINE__);
        scbm.set(cr1);
        CheckComplex (scbm(3, 2),cr1,  "scbmatrix = tcomplex", os, __LINE__);

        rv.assign(a2);
        CheckReal (rv[3],      a2[2], "rvector = treal*",      os, __LINE__);
        rm.assign(a2);
        CheckReal (rm(1, 3),   a2[4], "rmatrix = treal*",      os, __LINE__);
        srm.assign(a2);
        CheckReal (srm(3, 3),  a2[8], "srmatrix = treal*",     os, __LINE__);
        srbm.assign(a2);
        CheckReal (srbm(1, 1), a2[2], "srbmatrix = treal*",    os, __LINE__);
        CheckReal (srbm(2, 3), a2[9], "srbmatrix = treal*",    os, __LINE__);

        cv.assign(c1);
        CheckComplex (cv[3],      c1[2], "cvector = tcomplex*",   os, __LINE__);
        cm.assign(c1);
        CheckComplex (cm(1,3),    c1[4], "cmatrix = tcomplex*",   os, __LINE__);
        scm.assign(c1);
        CheckComplex (scm(3, 3),  c1[8], "scmatrix = tcomplex*",  os, __LINE__);
        scbm.assign(c1);
        CheckComplex (scbm(1, 1), c1[2], "scbmatrix = tcomplex*",    os, __LINE__);
        CheckComplex (scbm(2, 3), c1[9], "scbmatrix = tcomplex*",    os, __LINE__);

        {   // sub-assignment
            rvector rv2(4);
            rv2.randomize((treal) -3., (treal) 2.);
            rv.assign(3, rv2);
            CheckReal (rv[3], rv2[1], "rvector subvector assignment",    os, __LINE__);
            CheckReal (rv[6], rv2[4], "rvector subvector assignment",    os, __LINE__);

            cvector cv2(4);
            cv2.randomize_real((treal) -3., (treal) 2.);
            cv2.randomize_imag((treal) -3., (treal) 2.);
            cv.assign(3, cv2);
            CheckComplex (cv[3], cv2[1], "cvector subvector assignment",    os, __LINE__);
            CheckComplex (cv[6], cv2[4], "cvector subvector assignment",    os, __LINE__);

            rm.randomize((treal) -3., (treal) 2.);
            rm2.assign(2, 2, rm);
            CheckReal (rm2(2,2), rm(1,1), "rmatrix submatrix assignment",    os, __LINE__);
            CheckReal (rm2(3,4), rm(2,3), "rmatrix submatrix assignment",    os, __LINE__);

            srmatrix srm (5);
            srm.randomize((treal) -3., (treal) 2.);
            srm.assign(2, 2, rm);
            CheckReal (srm(2,2), rm(1,1), "srmatrix submatrix assignment",    os, __LINE__);
            CheckReal (srm(3,4), rm(2,3), "srmatrix submatrix assignment",    os, __LINE__);

            cm.randomize_real((treal) -3., (treal) 2.);
            cm.randomize_imag((treal) -3., (treal) 2.);
            cm2.assign(2, 2, cm);
            CheckComplex (cm2(2,2), cm(1,1), "cmatrix submatrix assignment",    os, __LINE__);
            CheckComplex (cm2(3,4), cm(2,3), "cmatrix submatrix assignment",    os, __LINE__);

            scmatrix scm (5);
            scm.randomize_real((treal) -3., (treal) 2.);
            scm.randomize_imag((treal) -3., (treal) 2.);
            scm.assign(2, 2, cm);
            CheckComplex (scm(2,2), cm(1,1), "scmatrix submatrix assignment",    os, __LINE__);
            CheckComplex (scm(3,4), cm(2,3), "scmatrix submatrix assignment",    os, __LINE__);

            int ns = srs1.msize();
            srs1.resize(5);
            srs2.randomize((treal) -3., (treal) 2.);
            srs1.assign(3,srs2);
            CheckReal (srs1(3,3), srs2(1,1), "srsmatrix submatrix assignment",    os, __LINE__);
            CheckReal (srs1(4,5), srs2(2,3), "srsmatrix submatrix assignment",    os, __LINE__);
            srs1.resize(ns);

            ns = sch1.msize();
            sch1.resize(5);
            sch2.randomize_real((treal) -3., (treal) 2.);
            sch2.randomize_imag((treal) -3., (treal) 2.);
            sch1.assign(3,sch2);
            CheckComplex (sch1(3,3), sch2(1,1), "schmatrix submatrix assignment",    os, __LINE__);
            CheckComplex (sch1(4,5), sch2(2,3), "schmatrix submatrix assignment",    os, __LINE__);
            sch1.resize(ns);
        }

        srs2 = srs1;
        CheckBool (srs1 == srs2, true,            "srsmatrix ==",   os, __LINE__);
        srs2.set(2,3,srs2(2,3) + (treal)0.000001);
        CheckBool (srs1 == srs2, false,            "srsmatrix ==",   os, __LINE__);
        CheckBool (srs1 != srs2, true,            "srsmatrix !=",   os, __LINE__);

        sch2 = sch1;
        CheckBool (sch1 == sch2, true,            "schmatrix ==",   os, __LINE__);
        sch2.set(2,3,sch2(2,3) + tcomplex ((treal)0.000001, (treal)0.00001));
        CheckBool (sch1 == sch2, false,            "schmatrix ==",   os, __LINE__);
        CheckBool (sch1 != sch2, true,            "schmatrix !=",   os, __LINE__);



        rv1 = rv;
        CheckReal (rv1[3],       rv(3),        "rvector = rvector",                             os, __LINE__);
        CheckBool (rv1 == rv, true,            "rvector ==",                                    os, __LINE__);
        CheckBool (rv1 != rv, false,           "rvector !=",                                    os, __LINE__);
        rm1 = rm;
        CheckReal (rm1[2][3],    rm(2, 3),     "rmatrix = rmatrix",                             os, __LINE__);
        CheckBool (rm1 == rm, true,            "rmatrix ==",                                    os, __LINE__);
        CheckBool (rm1 != rm, false,           "rmatrix !=",                                    os, __LINE__);
        CheckBool (rm1[1] == rm[1], true,      "rmatrix = rmatrix, rm1[1] == rm[1]",            os, __LINE__);
        CheckBool (rm1(2) != rm(2), false,     "rmatrix = rmatrix, rm1(2) != rm(2)",            os, __LINE__);
        srm1 = srm;
        CheckReal (srm1[2][3],   srm(2, 3),    "srmatrix = srmatrix",                           os, __LINE__);
        CheckBool (srm1 == srm, true,          "srmatrix ==",                                   os, __LINE__);
        CheckBool (srm1 != srm, false,         "srmatrix !=",                                   os, __LINE__);
        CheckBool (srm1[1] == srm[1], true,    "srmatrix = srmatrix, srm1[1] == srm[1]",        os, __LINE__);
        CheckBool (srm1(2) != srm(2), false,   "srmatrix = srmatrix, srm1(2) != srm(2)",        os, __LINE__);
        srbm1 = srbm;
        CheckReal (srbm1[2][1],  srbm(2, 1),   "srbmatrix = srbmatrix",                         os, __LINE__);
        CheckBool (srbm1 == srbm, true,        "srbmatrix ==",                                  os, __LINE__);
        CheckBool (srbm1 != srbm, false,       "srbmatrix !=",                                  os, __LINE__);
        CheckBool (srbm1[1] == srbm[1], true,  "srbmatrix = srbmatrix, srbm1[1] == srbm[1]",    os, __LINE__);
        CheckBool (srbm1(2) != srbm(2), false, "srbmatrix = srbmatrix, srbm1(2) != srbm(2)",    os, __LINE__);

        cv1 = cv;
        CheckComplex  (cv1[4],   cv(4),        "cvector = cvector",                             os, __LINE__);
        CheckBool (cv1 == cv, true,            "cvector ==",                                    os, __LINE__);
        CheckBool (cv1 != cv, false,           "cvector !=",                                    os, __LINE__);
        cm1 = cm;
        CheckComplex  (cm1[2][1],   cm(2,1),   "cmatrix = cmatrix",                             os, __LINE__);
        CheckBool (cm1 == cm, true,            "cmatrix ==",                                    os, __LINE__);
        CheckBool (cm1 != cm, false,           "cmatrix !=",                                    os, __LINE__);
        CheckBool (cm1[1] == cm[1], true,      "cmatrix = cmatrix, cm1[1] == cm[1]",            os, __LINE__);
        CheckBool (cm1(2) != cm(2), false,     "cmatrix = cmatrix, cm1(2) != cm(2)",            os, __LINE__);
        scm1 = scm;
        CheckComplex  (scm1[2][1],   scm(2,1), "scmatrix = scmatrix",                           os, __LINE__);
        CheckBool (scm1 == scm, true,          "scmatrix ==",                                   os, __LINE__);
        CheckBool (scm1 != scm, false,         "scmatrix !=",                                   os, __LINE__);
        CheckBool (scm1[1] == scm[1], true,    "scmatrix = scmatrix, scm1[1] == scm[1]",        os, __LINE__);
        CheckBool (scm1(2) != scm(2), false,   "scmatrix = scmatrix, scm1(2) != scm(2)",        os, __LINE__);
        scbm1 = scbm;
        CheckComplex (scbm1[2][1],  scbm(2, 1),   "scbmatrix = scbmatrix",                      os, __LINE__);
        CheckBool (scbm1 == scbm, true,        "scbmatrix ==",                                  os, __LINE__);
        CheckBool (scbm1 != scbm, false,       "scbmatrix !=",                                  os, __LINE__);
        CheckBool (scbm1[1] == scbm[1], true,  "scbmatrix = scbmatrix, scbm1[1] == scbm[1]",    os, __LINE__);
        CheckBool (scbm1(2) != scbm(2), false, "scbmatrix = scbmatrix, scbm1(2) != scbm(2)",    os, __LINE__);

//        rv2 = rv + rv1;   // wouldn't work because rv and rv1 share the same array!
        rv3.resize(10);
        rv3 = rv + rv1;
        CheckReal (rv3[1],   rv(1) + rv1[1],                     "rvector + rvector",           os, __LINE__);
        CheckReal (rv3[10],  rv(10) + rv1[10],                   "rvector + rvector",           os, __LINE__);
        rv3 = rv - rv1;
        CheckReal (rv3[1],   rv(1) - rv1[1],                     "rvector - rvector",           os, __LINE__);
        CheckReal (rv3[10],  rv(10) - rv1[10],                   "rvector - rvector",           os, __LINE__);
        cv3 = cv + cv1;
        CheckComplex (cv3[1],   cv(1) + cv1[1],                  "cvector + cvector",           os, __LINE__);
        CheckComplex (cv3[10],  cv(10) + cv1[10],                "cvector + cvector",           os, __LINE__);
        cv3 = cv - cv1;
        CheckComplex (cv3[1],   cv(1) - cv1[1],                  "cvector - cvector",           os, __LINE__);
        CheckComplex (cv3[10],  cv(10) - cv1[10],                "cvector - cvector",           os, __LINE__);
        rm2.resize(2,3);
        rm = rm1 + rm2;
        CheckReal (rm[1][1],   rm1(1,1) + rm2(1,1),              "rmatrix + rmatrix",           os, __LINE__);
        CheckReal (rm[2].norm(), (rm1[2] + rm2[2]).norm(),       "rmatrix + rmatrix",           os, __LINE__);
        rm = rm1 - rm2;
        CheckReal (rm[1][1],   rm1(1,1) - rm2(1,1),              "rmatrix - rmatrix",           os, __LINE__);
        CheckReal (rm(3).norm(), (rm1(3) - rm2(3)).norm(),       "rmatrix - rmatrix",           os, __LINE__);
        cm2.resize(2,3);
        cm = cm1 + cm2;
        CheckComplex (cm[1][1], cm1(1,1) + cm2(1,1),             "cmatrix + cmatrix",           os, __LINE__);
        CheckComplex (cm[2].norm(), (cm1[2] + cm2[2]).norm(),    "cmatrix + cmatrix",           os, __LINE__);
        cm = cm1 - cm2;
        CheckComplex (cm[1][1],   cm1(1,1) - cm2(1,1),           "cmatrix - cmatrix",           os, __LINE__);
        CheckComplex (cm(3).norm(), (cm1(3) - cm2(3)).norm(),    "cmatrix - cmatrix",           os, __LINE__);
        srm = srm1 + srm2;
        CheckReal (srm[1][1],   srm1(1,1) + srm2(1,1),           "srmatrix + srmatrix",         os, __LINE__);
        CheckReal (srm[3].norm(), (srm1[3] + srm2[3]).norm(),    "srmatrix + srmatrix",         os, __LINE__);
        srm = srm1 - srm2;
        CheckReal (srm[1][1],   srm1(1,1) - srm2(1,1),           "srmatrix - srmatrix",         os, __LINE__);
        CheckReal (srm(3).norm(), (srm1(3) - srm2(3)).norm(),    "srmatrix - srmatrix",         os, __LINE__);
        scm = scm1 + scm2;
        CheckComplex (scm[1][1],   scm1(1,1) + scm2(1,1),        "scmatrix + scmatrix",         os, __LINE__);
        CheckComplex (scm(3).norm(), (scm1(3) + scm2(3)).norm(), "scmatrix + scmatrix",         os, __LINE__);
        scm = scm1 - scm2;
        CheckComplex (scm[1][1],   scm1(1,1) - scm2(1,1),        "scmatrix - scmatrix",         os, __LINE__);
        CheckComplex (scm(3).norm(), (scm1(3) - scm2(3)).norm(), "scmatrix - scmatrix",         os, __LINE__);

        srbmatrix srbm2 (a2, 4, 1, 2);
        srbm = srbm1 + srbm2;
        CheckReal (srbm[1][1],   srbm1(1,1) + srbm2(1,1),        "srbmatrix + srbmatrix",       os, __LINE__);
        CheckReal (srbm(3).norm(), (srbm1(3) + srbm2(3)).norm(), "srbmatrix + srbmatrix",       os, __LINE__);
        srbm = srbm1 - srbm2;
        CheckReal (srbm[1][1],   srbm1(1,1) - srbm2(1,1),        "srbmatrix - srbmatrix",       os, __LINE__);
        CheckReal (srbm(3).norm(), (srbm1(3) - srbm2(3)).norm(), "srbmatrix - srbmatrix",       os, __LINE__);

        scbmatrix scbm2 (c1, 4, 1, 2);
        scbm = scbm1 + scbm2;
        CheckComplex (scbm[1][1],   scbm1(1,1) + scbm2(1,1),        "scbmatrix + scbmatrix",       os, __LINE__);
        CheckReal    (scbm(3).norm(), (scbm1(3) + scbm2(3)).norm(), "scbmatrix + scbmatrix",       os, __LINE__);
        scbm = scbm1 - scbm2;
        CheckComplex (scbm[1][1],   scbm1(1,1) - scbm2(1,1),        "scbmatrix - scbmatrix",       os, __LINE__);
        CheckReal    (scbm(3).norm(), (scbm1(3) - scbm2(3)).norm(), "scbmatrix - scbmatrix",       os, __LINE__);

        srs2 = srs1;
        rs1 = srs1(1,2);
        CheckReal    ((srs1 + srs2)(1,2), rs1 + rs1, "srsmatrix + srsmatrix",       os, __LINE__);
        CheckReal    ((srs1 - srs2).norm(), (treal)0., "srsmatrix - srsmatrix",       os, __LINE__);

        sch2 = sch1;
        cs1 = sch1(1,2);
        CheckComplex ((sch1 + sch2)(1,2), cs1 + cs1, "schmatrix + schmatrix",       os, __LINE__);
        CheckReal    ((sch1 - sch2).norm(), (treal)0., "schmatrix - schmatrix",       os, __LINE__);


        int n1 = -2;
        r1     = -2.;
        rv1 = rv * r1;
        rv3 = n1 * rv;
        CheckReal    (rv3[3],   rv1[3],         "rvector * number",                             os, __LINE__);
        rv3 = r1 * rv;
        CheckReal    (rv3[3],   rv1[3],         "rvector * number",                             os, __LINE__);
        cv1 = cv * r1;
        cv3 = n1 * cv;
        CheckComplex (cv3[3],   cv1[3],         "cvector * number",                             os, __LINE__);
        cv3 = r1 * cv;
        CheckComplex (cv3[3],   cv1[3],         "cvector * number",                             os, __LINE__);
        rm1 = rm * r1;
        rm2 = n1 * rm;
        CheckReal    (rm1(2,3), rm2(2,3),       "rmatrix * number",                             os, __LINE__);
        rm2 = r1 * rm;
        CheckReal    (rm1(2,3), rm2(2,3),       "rmatrix * number",                             os, __LINE__);
        cm1 = cm * r1;
        cm2 = n1 * cm;
        CheckComplex (cm2(2,3), cm1(2,3),       "cmatrix * number",                             os, __LINE__);
        cm2 = r1 * cm;
        CheckComplex (cm2(2,3), cm1(2,3),       "cmatrix * number",                             os, __LINE__);
        srm1 = srm * r1;
        srm2 = n1 * srm;
        CheckReal    (srm1(2,3), srm2(2,3),     "srmatrix * number",                            os, __LINE__);
        srm2 = r1 * srm;
        CheckReal    (srm1(2,3), srm2(2,3),     "srmatrix * number",                            os, __LINE__);
        scm.assign(c1);
        scm1 = scm * r1;
        scm2 = n1 * scm;
        CheckComplex (scm2(2,3), scm1(2,3),     "scmatrix * number",                            os, __LINE__);
        scm2 = r1 * scm;
        CheckComplex (scm2(2,3), scm1(2,3),     "scmatrix * number",                            os, __LINE__);
        srbm1 = srbm * r1;
        srbm2 = n1 * srbm;
        CheckReal    (srbm1(2,3), srbm2(2,3),   "srbmatrix * number",                           os, __LINE__);
        srbm2 = r1 * srbm;
        CheckReal    (srbm1(2,3), srbm2(2,3),   "srbmatrix * number",                           os, __LINE__);
        scbm1 = scbm * r1;
        scbm2 = n1 * scbm;
        CheckComplex (scbm1(2,3), scbm2(2,3),   "scbmatrix * number",                           os, __LINE__);
        scbm2 = r1 * scbm;
        CheckComplex (scbm1(2,3), scbm2(2,3),   "scbmatrix * number",                           os, __LINE__);
        cr1 = r1;
        scbm2 = cr1 * scbm;
        CheckComplex (scbm1(2,3), scbm2(2,3),   "scbmatrix * number",                           os, __LINE__);

        rv1 = rv / r1;
        CheckReal    (rv1[10],   rv[10] / r1,       "rvector / number",                         os, __LINE__);
        cv1 = cv / r1;
        CheckComplex (cv1[10],   cv[10] / r1,       "cvector / number",                         os, __LINE__);
        rm1 = rm / r1;
        CheckReal    (rm1(2,3), rm(2,3) / r1,       "rmatrix / number",                         os, __LINE__);
        cm1 = cm / r1;
        CheckComplex (cm1(2,3), cm(2,3) / r1,       "cmatrix / number",                         os, __LINE__);

        srm1 = srm / r1;
        CheckReal    (srm1(2,3), srm(2,3) / r1,     "srmatrix / number",                        os, __LINE__);
        scm1 = scm / r1;
        CheckComplex (scm1(2,3), scm(2,3) / r1,     "scmatrix / number",                        os, __LINE__);
        srbm1 = srbm / r1;
        CheckReal    (srbm1(2,3), srbm(2,3) / r1,   "srbmatrix / number",                       os, __LINE__);
        scbm1 = scbm / r1;
        CheckComplex (scbm1(2,3), scbm(2,3) / r1,   "scbmatrix / number",                       os, __LINE__);
        scbm1 = scbm / cr1;
        CheckComplex (scbm1(2,3), scbm(2,3) / cr1,  "scbmatrix / number",                       os, __LINE__);

        cv1 = cv  * cr2;
        cv3 = cr2 * cv;
        CheckComplex (cv3[3],   cv1[3],         "cvector * cmplx number",                       os, __LINE__);
        cm1 = cm * cr2;
        cm2 = cr2 * cm;
        CheckComplex (cm2(2,3), cm1(2,3),       "cmatrix * cmplx number",                       os, __LINE__);
        scm1 = scm * cr2;
        scm2 = cr2 * scm;
        CheckComplex (scm2(2,3), scm1(2,3),     "scmatrix * cmplx number",                      os, __LINE__);
        scbm1 = scbm * cr2;
        scbm2 = cr2 * scbm;
        CheckComplex (scbm2(2,3), scbm1(2,3),   "scbmatrix * cmplx number",                     os, __LINE__);

        rv1 = - rv;
        CheckReal    (rv1[10],   - rv[10],      "- rvector",                                    os, __LINE__);
        cv1 = - cv;
        CheckComplex (cv1[10],   - cv[10],      "- cvector",                                    os, __LINE__);
        rm1 = - rm;
        CheckReal    (rm1(2,3), - rm(2,3),      "- rmatrix",                                    os, __LINE__);
        cm1 = - cm;
        CheckComplex (cm1(2,3), - cm(2,3),      "- cmatrix",                                    os, __LINE__);
        srm1 = - srm;
        CheckReal    (srm1(2,3), - srm(2,3),    "- srmatrix",                                   os, __LINE__);
        scm1 = - scm;
        CheckComplex (scm1(2,3), - scm(2,3),    "- scmatrix",                                   os, __LINE__);
        srbm.assign(a2);
        srbm1 = - srbm;
        CheckReal    (srbm1(2,3), - srbm(2,3),  "- srbmatrix",                                  os, __LINE__);
        scbm.assign(c1);
        scbm1 = - scbm;
        CheckComplex (scbm1(2,3), - scbm(2,3),  "- scbmatrix",                                  os, __LINE__);


        rv1.set((treal)1.17);
        rv2.set((treal)-0.31);
        rm2.set((treal)9.01);
        srbm1.set((treal)13.1);
        srbm2.set((treal)5.51);
        cv1.set(tcomplex (2,1));
        cv2.set(tcomplex (-1,3));
        cm2.set(tcomplex (-4,3));
        rv1.resize (2);
        rv2.resize (3);
        rv2.mult (rv1, rm2);

        CheckReal    (rv2[1], rv1 * rm2(1),  "mult",                                            os, __LINE__, dPessimisticSp);
        rv1.mult (rm2, rv2);
        CheckReal    (rv1[1], rv2 * rm2[1],  "mult",                                            os, __LINE__, dPessimisticSp);

        cv1.resize (2);
        cv2.resize (3);
        cv2.mult (cv1, cm2);

        CheckComplex (cv2[1], cv1 * cm2(1),  "mult",                                            os, __LINE__, dPessimisticSp);
        cv1.mult (cm2, cv2);
        CheckComplex (cv1[1], cv2 * cm2[1],  "mult",                                            os, __LINE__, dPessimisticSp);

        rv1.resize (3);
        rv1.mult (srm2, rv2);
        CheckReal    (rv1[1], rv2 * srm2[1],  "mult",                                           os, __LINE__, dPessimisticSp);
        rv2.mult (rv1, srm2);

        CheckReal    (rv2[1], rv1 * srm2(1),  "mult",                                           os, __LINE__, dPessimisticSp);

        cv1.resize (3);
        cv1.mult (scm2, cv2);
        CheckComplex (cv1[1], cv2 * scm2[1],  "mult",                                           os, __LINE__, dPessimisticSp);
        cv2.mult (cv1, scm2);
        CheckComplex (cv2[1], cv1 * scm2(1),  "mult",                                           os, __LINE__, dPessimisticSp);

        rv1.resize (4);
        rv2.resize (4);
        rv2.mult (rv1, srbm2);
        CheckReal    (rv2[1], rv1 * srbm2(1),  "mult",                                          os, __LINE__, dPessimisticSp);
        rv1.mult (srbm2, rv2);
        CheckReal    (rv1[1], rv2 * srbm2[1],  "mult",                                          os, __LINE__, dPessimisticSp);

        cv1.resize (4);
        cv2.resize (4);
        cv2.mult (cv1, scbm2);
        CheckComplex (cv2[1], cv1 * scbm2(1),  "mult",                                          os, __LINE__, dPessimisticSp);
        cv1.mult (scbm2, cv2);
        CheckComplex (cv1[1], cv2 * scbm2[1],  "mult",                                          os, __LINE__, dPessimisticSp);

        rm1.resize (3, 2);
        rm1[3].assign(a1);
        rm3.resize (2, 2);
        rm4.resize (3, 3);
        rm3.mult (rm2, rm1);
        CheckReal    (rm3(2,2), rm2[2] * rm1(2),  "mult",                                       os, __LINE__, dPessimisticSp);
        rm4.mult (rm1, rm2);
        CheckReal    (rm4(3,3), rm1[3] * rm2(3),  "mult",                                       os, __LINE__, dPessimisticSp);
        srm4.resize(3);
        srm4.mult (rm1, rm2);

        CheckReal    (srm4(3,3), rm1[3] * rm2(3), "mult",                                       os, __LINE__, dPessimisticSp);
        rm4.resize (3, 2);
        rm1.mult (srm4, rm4);
        CheckReal    (rm1(3,2),  srm4[3] * rm4(2), "mult",                                      os, __LINE__, dPessimisticSp);
        srbm1.resize(3);
        rm1.mult (srbm1, rm4);
        CheckReal    (rm1(3,2),  srbm1[3] * rm4(2), "mult",                                     os, __LINE__, dPessimisticSp);
        rm1.mult (~srbm1, rm4);
        CheckReal    (rm1(3,2),  srbm1(3) * rm4(2), "mult",                                     os, __LINE__, dPessimisticSp);
        srbm1.mult (rm1, rm2);
        CheckReal    (srbm1(2,2), rm1[2] * rm2(2),  "mult",                                     os, __LINE__, dPessimisticSp);

        r1 = (treal) -0.031;
        r2 = (treal) 0.319;
        rm1.randomize(1., 2.);
        rm2.randomize(0., 1.);
        rm3.randomize(0., 1.);
        rmatrix rm3_dub = rm3;

        rm3.gemm (rm2, false, rm1, false, r1, r2);
        CheckReal    ((rm3 - (rm2 * rm1 * r1 + rm3_dub * r2)).norm2(), (treal) 0.,  "gemm",      os, __LINE__, dVeryPessimisticSp);
        rm3_dub = rm3;
        rm3 << ~rm3;
        rm3.gemm (rm1, true, rm2, true, r1, r2);
        CheckReal    ((~rm3 - (rm2 * rm1 * r1 + rm3_dub * r2)).norm2(), (treal) 0.,  "gemm",      os, __LINE__, dVeryPessimisticSp);

        srbm1.randomize((treal) -1., (treal) 3.);
        rmatrix rm1_dub = rm1;
        rm1.gemm (srbm1, false, rm4, false, r1, r2);
        CheckReal    ((rm1 - (srbm1 * rm4 * r1 + rm1_dub * r2)).norm2(), (treal) 0.,  "gemm",      os, __LINE__, dVeryPessimisticSp);




        cm1.resize (3, 2);
        cm1[3].assign(c1);
        cmatrix cm3 (2, 2), cm4 (3,3);
        cm3.assign(c2);
        cm3.mult (cm2, cm1);
        CheckComplex (cm3(2,2), cm2[2] * cm1(2),  "mult",                                       os, __LINE__, dPessimisticSp);
        cm4.mult (cm1, cm2);
        CheckComplex (cm4(3,3), cm1[3] * cm2(3),  "mult",                                       os, __LINE__, dPessimisticSp);
        scm4.resize(3);
        scm4.mult (cm1, cm2);
        CheckComplex (scm4(3,3), cm1[3] * cm2(3), "mult",                                       os, __LINE__, dPessimisticSp);
        cm4.resize (3, 2);
        cm1.mult (scm4, cm4);
        CheckComplex (cm1(3,2),  scm4[3] * cm4(2), "mult",                                      os, __LINE__, dPessimisticSp);
        scbm.resize(3);
        scbm.set(tcomplex ((treal) 1.23, (treal) -0.912));
        cm1.mult (scbm, cm4);
        CheckComplex (cm1(3,2),  scbm[3] * cm4(2), "mult",                                     os, __LINE__, dPessimisticSp);
        cm1.mult (~scbm, cm4);
        CheckComplex (cm1(3,2),  ~(scbm(3)) * cm4(2), "mult",                                     os, __LINE__, dPessimisticSp);
        scbm1.resize(3);
        scbm1.mult (cm1, cm2);
        CheckComplex (scbm1(2,2), cm1[2] * cm2(2),  "mult",                                     os, __LINE__, dPessimisticSp);


        cm1.randomize_real(0., 1.);
        cm2.randomize_real(0., 1.);
        cm3.randomize_real(0., 1.);
        scbm.randomize_real(0., 1.);
        cmatrix cm3_dub = cm3;
        cm3.gemm (cm2, false, cm1, false, cr1, cr2);
        CheckReal    ((cm3 - (cm2 * cm1 * cr1 + cm3_dub * cr2)).norm(), (treal) 0.,  "gemm",      os, __LINE__, dPessimisticSp);
        cmatrix cm1_dub = cm1;
        cm1.gemm (scbm, false, cm4, false, cr1, cr2);
        CheckReal    ((cm1 - (scbm * cm4 * cr1 + cm1_dub * cr2)).norm(), (treal) 0.,  "gemm",      os, __LINE__, dPessimisticSp);

        cr1 = tcomplex ((treal)-1.14,(treal)3.22);
        cr2 = tcomplex ((treal)2.04,(treal)-4.2);
        cm1_dub << cm1;
        cm1.conj();
        cm1.gemm (cm4, true, scbm, true, cr1, cr2);
        CheckReal    ((~cm1 - (scbm * cm4 * conj(cr1) + cm1_dub * conj(cr2))).norm2(), (treal) 0.,  "gemm", os, __LINE__, dVeryPessimisticSp);
        cm1.conj();

        rv1.randomize((treal)0., (treal)1.);
        rv2.randomize((treal)0., (treal)1.);
        CheckReal    (rv1 * rv2, rv1[1]*rv2[1]+rv1[2]*rv2[2]+rv1[3]*rv2[3]+rv1[4]*rv2[4],  "scalar product", os, __LINE__, dPessimisticSp);
        cv1.randomize_real((treal)0., (treal)1.);

        cv1.randomize_imag((treal)0., (treal)1.);
        cv2.randomize_real((treal)0., (treal)1.);
        cv2.randomize_imag((treal)0., (treal)1.);
        CheckComplex (cv1 * cv2, cv1[1]*cv2[1]+cv1[2]*cv2[2]+cv1[3]*cv2[3]+cv1[4]*cv2[4],  "scalar product",  os, __LINE__, dPessimisticSp);
        CheckComplex (cv1 % cv2, conj(cv1[1])*cv2[1]+conj(cv1[2])*cv2[2]+conj(cv1[3])*cv2[3]+conj(cv1[4])*cv2[4],  "scalar product, conj",  os, __LINE__, dPessimisticSp);

        CheckReal    ((rm1[2] - (~rm1)(2)).norm(), (treal) 0.,  "~",                            os, __LINE__);

        cvector cm1_2_conj (cm1[2].size());
        cm1_2_conj = cm1[2];
        cm1_2_conj.conj();
        CheckReal    ((cm1_2_conj - (~cm1)(2)).norm(), (treal) 0.,  "~",                     os, __LINE__);
        CheckReal    ((srbm1[2] - (~srbm1)(2)).norm(), (treal) 0.,  "~",                        os, __LINE__);
        CheckReal    ((~(scbm1[2]) - (~scbm1)(2)).norm(), (treal) 0.,  "~",                     os, __LINE__);

        rv1.resize (3);
        rv2.resize (2);
        rv1 = rm1 * rv2;
        CheckReal    (rv1[3], rv2 * rm1[3],  "rmatrix * rvector",                               os, __LINE__, dPessimisticSp);
        rv2 = rv1 * rm1;
        CheckReal    (rv2[2], rv1 * rm1(2),  "rvector * rmatrix",                               os, __LINE__, dPessimisticSp);
        cv1.resize (3);
        cv2.resize (2);
        cv1 = cm1 * cv2;
        CheckComplex (cv1[3], cv2 * cm1[3],  "cmatrix * cvector",                               os, __LINE__, dPessimisticSp);
        cv2 = cv1 * cm1;
        CheckComplex (cv2[2], cv1 * cm1(2),  "cvector * cmatrix",                               os, __LINE__, dPessimisticSp);

        rv2.resize (3);
        rv2 = srm4 * rv1;
        CheckReal    (rv2[3], rv1 * srm4[3],  "srmatrix * rvector",                             os, __LINE__, dPessimisticSp);
        rv2 = rv1 * srm4;
        CheckReal    (rv2[3], rv1 * srm4(3),  "rvector * srmatrix",                             os, __LINE__, dPessimisticSp);
        cv2.resize (3);
        cv2 = scm4 * cv1;
        CheckComplex (cv2[3], cv1 * scm4[3],  "scmatrix * cvector",                             os, __LINE__, dPessimisticSp);
        cv2 = cv1 * scm4;
        CheckComplex (cv2[3], cv1 * scm4(3),  "cvector * scmatrix",                             os, __LINE__, dPessimisticSp);

        srbm1.normalize();
        rv1.normalize();
        rv2 = srbm1 * rv1;
        CheckReal    (rv2[3], rv1 * srbm1[3],  "srbmatrix * rvector",                           os, __LINE__);
        rv2 = rv1 * srbm1;
        CheckReal    (rv2[3], rv1 * srbm1(3),  "rvector * srbmatrix",                           os, __LINE__);

        scbm1.normalize();
        cv1.normalize();
        cv2 = scbm1 * cv1;
        CheckComplex (cv2[3], cv1 * scbm1[3],  "scbmatrix * cvector",                           os, __LINE__);
        cv2 = cv1 * scbm1;
        CheckComplex (cv2[3], cv1 * scbm1(3),  "cvector * scbmatrix",                           os, __LINE__);

        rv2.resize (2);
        rm1 = rv1.rank1update (rv2);
        CheckReal    (rm1(3,1), rv1[3] * rv2[1],  "rank1update",                                os, __LINE__);
        rm1.rank1update (rv1, rv2);
        CheckReal    (rm1(3,2), rv1[3] * rv2[2],  "rank1update",                                os, __LINE__);

        cv2.resize (2);
        cv1.normalize();
        cv2.normalize();
        cm1 = cv1.rank1update_u (cv2);
        CheckComplex (cm1(3,1), cv1[3] * cv2[1],  "rank1update_u",                              os, __LINE__);
        cm1.rank1update_u (cv1, cv2);
        CheckComplex (cm1(3,2), cv1[3] * cv2[2],  "rank1update_u",                              os, __LINE__);
        cm1 = cv1.rank1update_c (cv2);
        CheckComplex (cm1(3,1), cv1[3] * conj (cv2[1]),  "rank1update_c",                       os, __LINE__);
        cm1.rank1update_c (cv1, cv2);
        CheckComplex (cm1(3,2), cv1[3] * conj (cv2[2]),  "rank1update_c",                       os, __LINE__);

        srm4.assign(a3);
        srm4(3,3) = -(treal) 1.;
        srm4.normalize();
        CheckReal    (srm4.cond(), (treal) 1. / (srm4.norminf() * srm4.inv().norminf()),  "cond", os, __LINE__, dVeryPessimisticSp);
        CheckReal    (srm4.det(), srm4(1,1) * srm4(2,2) * srm4(3,3) -
                                  srm4(1,1) * srm4(2,3) * srm4(3,2) -
                                  srm4(1,2) * srm4(2,1) * srm4(3,3) +
                                  srm4(1,2) * srm4(2,3) * srm4(3,1) +
                                  srm4(1,3) * srm4(2,1) * srm4(3,2) -
                                  srm4(1,3) * srm4(2,2) * srm4(3,1),  "det", os, __LINE__);

        scm4.assign(c2);
        scm4.normalize();
        CheckReal    (scm4.cond(), (treal) 1. / (scm4.norminf() * scm4.inv().norminf()),  "cond", os, __LINE__, dVeryPessimisticSp);
        CheckComplex (scm4.det(), scm4(1,1) * scm4(2,2) * scm4(3,3) -
                                  scm4(1,1) * scm4(2,3) * scm4(3,2) -
                                  scm4(1,2) * scm4(2,1) * scm4(3,3) +
                                  scm4(1,2) * scm4(2,3) * scm4(3,1) +
                                  scm4(1,3) * scm4(2,1) * scm4(3,2) -
                                  scm4(1,3) * scm4(2,2) * scm4(3,1),  "det", os, __LINE__);


        r1 = (treal) 2.;
        rv1.resize (4);
        rm1.resize (4, 4);
        srbm1.resize (4);
        rv1.set(1.);
        rm1 << eye_real(4);
        srm4 << eye_real(4);
        srbm1 << srbmatrix (eye_real(4), 0, 0);

        CheckReal    (rv1.norm(), r1,  "rvector norm", os, __LINE__);
        CheckReal    (rm1.norm(), r1,  "rmatrix norm", os, __LINE__);
        CheckReal    (srm4.norm(), r1,  "srmatrix norm", os, __LINE__);
        CheckReal    (srbm1.norm(), r1,  "srbmatrix norm", os, __LINE__);

        r1 = (treal) 2. * (treal) sqrt ((treal) 2.);
        cv1.resize (4);


        cm1.resize (4, 4);
        scm1.resize (4);
        cv1.set(tcomplex (1, 1));
        cm1 << scmatrix (cv1);
        scm1 = cm1;
        scbm2 = scbmatrix(cm1, scbm2.lsize(), scbm2.usize());
        CheckReal    (cv1.norm(), r1,  "cvector norm", os, __LINE__);
        CheckReal    (cm1.norm(), r1,  "cmatrix norm", os, __LINE__);
        CheckReal    (scm1.norm(), r1,  "scmatrix norm", os, __LINE__);
        CheckReal    (scbm2.norm(), r1,  "scbmatrix norm", os, __LINE__);

        // mix
        scbm2.set(tcomplex((treal)1.23, (treal)-0.977));

        cm1 = scbm2;
        CheckComplex (cm1(2,3), scbm2(2,3),  "mix cmatrix  scbm", os, __LINE__);
        CheckComplex (cm1(4,1), scbm2(4,1),  "mix cmatrix  scbm", os, __LINE__);

        cm1 = cm1 + scbm2;
        cm1 += scbm2;
        CheckComplex (cm1(2,3), scbm2(2,3) * 3,  "mix cmatrix  scbm", os, __LINE__);
        CheckComplex (cm1(4,1), scbm2(4,1) * (treal)3.,  "mix cmatrix  scbm", os, __LINE__);
        CheckComplex (cm1(1,2), 3 * scbm2(2,3),  "mix cmatrix  scbm", os, __LINE__);
        CheckComplex (cm1(2,1), 3. * scbm2(2,1),  "mix cmatrix  scbm", os, __LINE__);

        rm1 = srbm2;
        CheckReal    (rm1(2,3), srbm2(2,3),  "mix rmatrix  srbm", os, __LINE__);
        CheckReal    (rm1(4,1), srbm2(4,1),  "mix rmatrix  srbm", os, __LINE__);

        rm1 = rm1 + srbm2;
        rm1 += srbm2;
        CheckReal    (rm1(2,3), srbm2(2,3) * 3.,  "mix rmatrix  srbm", os, __LINE__);
        CheckReal    (rm1(4,1), 3. * srbm2(4,1),  "mix matrix  srbm", os, __LINE__);
        CheckReal    (rm1(4,1), 3 * srbm2(4,1),  "mix rmatrix  srbm", os, __LINE__);
        CheckReal    (rm1(2,3), srbm2(2,3) * 3,  "mix rmatrix  srbm", os, __LINE__);

        scbm1.resize(4);
        for (j = 1; j <= 4; j++)
        {
            for (i = 1; i <= 4; i++)
            {
                rm1(i,j)  = - (treal) ((j - 1) * 4 + i);
                srm4(i,j) = - (treal) ((j - 1) * 4 + i);
                cm1(i,j)  = - (treal) ((j - 1) * 4 + i);
                scm1(i,j) = - (treal) ((j - 1) * 4 + i);
            }
            srbm1(j,j) = (treal) j;
            scbm1(j,j) = (tcomplex) (treal)j;
        }

        CheckReal    (rm1.norm1(),   (treal) (13 + 14 + 15 + 16),  "rmatrix norm1", os, __LINE__);
        CheckReal    (srm4.norm1(),  (treal) (13 + 14 + 15 + 16),  "srmatrix norm1", os, __LINE__);
        CheckReal    (srbm1.norm1(), (treal) 4,                    "srbmatrix norm1", os, __LINE__);
        CheckReal    (cm1.norm1(),   (treal) (13 + 14 + 15 + 16),  "cmatrix norm1", os, __LINE__);
        CheckReal    (scm1.norm1(),  (treal) (13 + 14 + 15 + 16),  "scmatrix norm1", os, __LINE__);
        CheckReal    (scbm1.norm1(), (treal) 4,                    "scbmatrix norm1", os, __LINE__);

        CheckReal    (rm1.norminf(),   (treal) (4 + 8 + 12 + 16),  "rmatrix norminf", os, __LINE__);
        CheckReal    (srm4.norminf(),  (treal) (4 + 8 + 12 + 16),  "srmatrix norminf", os, __LINE__);
        CheckReal    (srbm1.norminf(), (treal) 4,                  "srbmatrix norminf", os, __LINE__);
        CheckReal    (cm1.norminf(),   (treal) (4 + 8 + 12 + 16),  "cmatrix norminf", os, __LINE__);
        CheckReal    (scm1.norminf(),  (treal) (4 + 8 + 12 + 16),  "scmatrix norminf", os, __LINE__);
        CheckReal    (scbm1.norminf(), (treal) 4,                  "scbmatrix norminf", os, __LINE__);

        CheckReal    (eye_real(6)(6,6),   (treal) 1.,  "eye_real", os, __LINE__);
        CheckComplex (eye_complex(6)(6,6),   tcomplex (1, 0),  "eye_complex", os, __LINE__);

        rv2.resize (4);
        srmatrix rmU(4), rmVH(4);
        rv1 = srm4.svd (rmU, rmVH);
        rv2.svd (srm4, rmU, rmVH);
        CheckBool    (rv1 == rv2,  true,  "srmatrix svd", os, __LINE__);
        srm1 << srmatrix (rv1);
        CheckReal    ((srm4 * ~rmVH - rmU * srm1).norm(),  (treal) 0.,  "srmatrix svd", os, __LINE__, dPessimisticSp);
        CheckReal    ((~srm4 * rmU - ~(srm1 * rmVH)).norm(),  (treal) 0.,  "srmatrix svd", os, __LINE__, dPessimisticSp);

        rv1 = srbm2.svd (rmU, rmVH);
        rv2.svd (srbm2);

        CheckReal    ((rv1 - rv2).norm(),  (treal) 0.,  "srbmatrix svd", os, __LINE__, dVeryPessimisticSp);
        rv2.svd (srbm2, rmU, rmVH);
        srm1 << srmatrix (rv1);
        CheckReal    ((srbm2 * ~rmVH - rmU * srm1).norm(),  (treal) 0.,  "srbmatrix svd", os, __LINE__, dPessimisticSp);
        CheckReal    ((~srbm2 * rmU - ~(srm1 * rmVH)).norm(),  (treal) 0.,  "srbmatrix svd", os, __LINE__, dPessimisticSp);


        // test case from Martin
        // http://www.vni.com/products/jmsl/v25/api/com/imsl/math/SVDEx1.html
        rmatrix A(6,4);
        A(1,1) = 1;
        A(1,2) = 2;
        A(1,3) = 1;
        A(1,4) = 4;
        A(2,1) = 3;
        A(2,2) = 2;
        A(2,3) = 1;
        A(2,4) = 3;
        A(3,1) = 4;
        A(3,2) = 3;
        A(3,3) = 1;
        A(3,4) = 4;
        A(4,1) = 2;
        A(4,2) = 1;
        A(4,3) = 3;
        A(4,4) = 1;
        A(5,1) = 1;
        A(5,2) = 5;
        A(5,3) = 2;
        A(5,4) = 2;
        A(6,1) = 1;
        A(6,2) = 2;
        A(6,3) = 2;
        A(6,4) = 3;
        srmatrix U(6), V(4);
        const rvector singVal = A.svd(U,V);

        rmatrix singValM (A);
        singValM.set(0.);
        singValM(1,1) = singVal(1);
        singValM(2,2) = singVal(2);
        singValM(3,3) = singVal(3);
        singValM(4,4) = singVal(4);

        CheckReal    ((A * ~V - U * singValM).norm(),  (treal) 0.,  "rmatrix svd", os, __LINE__, dPessimisticSp);
        CheckReal    ((~A * U - ~(singValM * V)).norm(),  (treal) 0.,  "rmatrix svd", os, __LINE__, dPessimisticSp);

        CheckReal (singVal[1], (treal) 1.148501791155974e+001, "rmatrix svd", os, __LINE__, dPessimisticSp);
        CheckReal (singVal[2], (treal) 3.269751214412497e+000, "rmatrix svd", os, __LINE__, dPessimisticSp);
        CheckReal (singVal[3], (treal) 2.653356162007834e+000, "rmatrix svd", os, __LINE__, dPessimisticSp);
        CheckReal (singVal[4], (treal) 2.088729672440923e+000, "rmatrix svd", os, __LINE__, dPessimisticSp);
        CheckReal (cvm::_abs(U(1, 1)), cvm::_abs((treal)-0.38047558632), "rmatrix svd", os, __LINE__, dVeryPessimisticSp);
        CheckReal (cvm::_abs(U(1, 2)), cvm::_abs((treal)-0.11967099264), "rmatrix svd", os, __LINE__, dVeryPessimisticSp);
        CheckReal (cvm::_abs(U(1, 3)), cvm::_abs((treal)-0.43908282438), "rmatrix svd", os, __LINE__, dVeryPessimisticSp);
        CheckReal (cvm::_abs(U(1, 4)), cvm::_abs((treal)0.56539958591), "rmatrix svd", os, __LINE__, dVeryPessimisticSp);
        CheckReal (cvm::_abs(U(1, 5)), cvm::_abs((treal)0.024311516146),"rmatrix svd", os, __LINE__, dVeryPessimisticSp);
        CheckReal (cvm::_abs(U(1, 6)), cvm::_abs((treal)-0.5725868611), "rmatrix svd", os, __LINE__, dVeryPessimisticSp);
        CheckReal (cvm::_abs(U(2, 1)), cvm::_abs((treal)-0.40375371317), "rmatrix svd", os, __LINE__, dVeryPessimisticSp);
        CheckReal (cvm::_abs(U(2, 2)), cvm::_abs((treal)-0.34511083711), "rmatrix svd", os, __LINE__, dVeryPessimisticSp);
        CheckReal (cvm::_abs(U(2, 3)), cvm::_abs((treal)0.05657618529),"rmatrix svd", os, __LINE__, dVeryPessimisticSp);
        CheckReal (cvm::_abs(U(2, 4)), cvm::_abs((treal)-0.21477557652), "rmatrix svd", os, __LINE__, dVeryPessimisticSp);
        CheckReal (cvm::_abs(U(2, 5)), cvm::_abs((treal)0.80890058873),"rmatrix svd", os, __LINE__, dVeryPessimisticSp);
        CheckReal (cvm::_abs(U(2, 6)), cvm::_abs((treal)0.11929741721),"rmatrix svd", os, __LINE__, dVeryPessimisticSp);
        CheckReal (cvm::_abs(U(3, 1)), cvm::_abs((treal)-0.54512048625), "rmatrix svd", os, __LINE__, dVeryPessimisticSp);
        CheckReal (cvm::_abs(U(3, 2)), cvm::_abs((treal)-0.42926489349), "rmatrix svd", os, __LINE__, dVeryPessimisticSp);
        CheckReal (cvm::_abs(U(3, 3)), cvm::_abs((treal)-0.051392692809), "rmatrix svd", os, __LINE__, dVeryPessimisticSp);
        CheckReal (cvm::_abs(U(3, 4)), cvm::_abs((treal)-0.43214416281), "rmatrix svd", os, __LINE__, dVeryPessimisticSp);
        CheckReal (cvm::_abs(U(3, 5)), cvm::_abs((treal)-0.57232764817), "rmatrix svd", os, __LINE__, dVeryPessimisticSp);
        CheckReal (cvm::_abs(U(3, 6)), cvm::_abs((treal)0.040330924871),"rmatrix svd", os, __LINE__, dVeryPessimisticSp);
        CheckReal (cvm::_abs(U(4, 1)), cvm::_abs((treal)-0.264784294), "rmatrix svd", os, __LINE__, dVeryPessimisticSp);
        CheckReal (cvm::_abs(U(4, 2)), cvm::_abs((treal)0.068319525327),"rmatrix svd", os, __LINE__, dVeryPessimisticSp);
        CheckReal (cvm::_abs(U(4, 3)), cvm::_abs((treal)0.88386086743),"rmatrix svd", os, __LINE__, dVeryPessimisticSp);
        CheckReal (cvm::_abs(U(4, 4)), cvm::_abs((treal)0.21525369818),"rmatrix svd", os, __LINE__, dVeryPessimisticSp);
        CheckReal (cvm::_abs(U(4, 5)), cvm::_abs((treal)-0.06252092259), "rmatrix svd", os, __LINE__, dVeryPessimisticSp);
        CheckReal (cvm::_abs(U(4, 6)), cvm::_abs((treal)-0.30621669907), "rmatrix svd", os, __LINE__, dVeryPessimisticSp);
        CheckReal (cvm::_abs(U(5, 1)), cvm::_abs((treal)-0.4463101123), "rmatrix svd", os, __LINE__, dVeryPessimisticSp);
        CheckReal (cvm::_abs(U(5, 2)), cvm::_abs((treal)0.81682762328),"rmatrix svd", os, __LINE__, dVeryPessimisticSp);
        CheckReal (cvm::_abs(U(5, 3)), cvm::_abs((treal)-0.14189967506), "rmatrix svd", os, __LINE__, dVeryPessimisticSp);
        CheckReal (cvm::_abs(U(5, 4)), cvm::_abs((treal)-0.32126958427), "rmatrix svd", os, __LINE__, dVeryPessimisticSp);
        CheckReal (cvm::_abs(U(5, 5)), cvm::_abs((treal)0.062133782096),"rmatrix svd", os, __LINE__, dVeryPessimisticSp);
        CheckReal (cvm::_abs(U(5, 6)), cvm::_abs((treal)-0.079935268), "rmatrix svd", os, __LINE__, dVeryPessimisticSp);
        CheckReal (cvm::_abs(U(6, 1)), cvm::_abs((treal)-0.35462865661), "rmatrix svd", os, __LINE__, dVeryPessimisticSp);
        CheckReal (cvm::_abs(U(6, 2)), cvm::_abs((treal)0.10214739916),"rmatrix svd", os, __LINE__, dVeryPessimisticSp);
        CheckReal (cvm::_abs(U(6, 3)), cvm::_abs((treal)0.0043184439799),"rmatrix svd", os, __LINE__, dVeryPessimisticSp);
        CheckReal (cvm::_abs(U(6, 4)), cvm::_abs((treal)0.54580022185),"rmatrix svd", os, __LINE__, dVeryPessimisticSp);
        CheckReal (cvm::_abs(U(6, 5)), cvm::_abs((treal)-0.098794626562), "rmatrix svd", os, __LINE__, dVeryPessimisticSp);
        CheckReal (cvm::_abs(U(6, 6)), cvm::_abs((treal)0.74573957611),"rmatrix svd", os, __LINE__, dVeryPessimisticSp);

        CheckReal (cvm::_abs((~V)(1, 1)), cvm::_abs((treal)-4.442941288423535e-001), "rmatrix svd", os, __LINE__, dPessimisticSp);
        CheckReal (cvm::_abs((~V)(2, 1)), cvm::_abs((treal)-5.580672381903871e-001), "rmatrix svd", os, __LINE__, dPessimisticSp);
        CheckReal (cvm::_abs((~V)(3, 1)), cvm::_abs((treal)-3.243861032062802e-001), "rmatrix svd", os, __LINE__, dPessimisticSp);
        CheckReal (cvm::_abs((~V)(4, 1)), cvm::_abs((treal)-6.212385538433783e-001), "rmatrix svd", os, __LINE__, dPessimisticSp);
        CheckReal (cvm::_abs((~V)(1, 2)), cvm::_abs((treal)5.555312577999473e-001), "rmatrix svd", os, __LINE__, dVeryPessimisticSp);
        CheckReal (cvm::_abs((~V)(2, 2)), cvm::_abs((treal)-6.542987401123238e-001), "rmatrix svd", os, __LINE__, dPessimisticSp);
        CheckReal (cvm::_abs((~V)(3, 2)), cvm::_abs((treal)-3.513606455925113e-001), "rmatrix svd", os, __LINE__, dPessimisticSp);
        CheckReal (cvm::_abs((~V)(4, 2)), cvm::_abs((treal)3.739303103834293e-001),"rmatrix svd", os, __LINE__, dPessimisticSp);
        CheckReal (cvm::_abs((~V)(1, 3)), cvm::_abs((treal)-4.353789666739416e-001), "rmatrix svd", os, __LINE__, dPessimisticSp);
        CheckReal (cvm::_abs((~V)(2, 3)), cvm::_abs((treal)2.774569004588126e-001),"rmatrix svd", os, __LINE__, dPessimisticSp);
        CheckReal (cvm::_abs((~V)(3, 3)), cvm::_abs((treal)-7.320995334295977e-001), "rmatrix svd", os, __LINE__, dPessimisticSp);
        CheckReal (cvm::_abs((~V)(4, 3)), cvm::_abs((treal)4.444019542237462e-001),"rmatrix svd", os, __LINE__, dPessimisticSp);
        CheckReal (cvm::_abs((~V)(1, 4)), cvm::_abs((treal)-5.517543874418699e-001), "rmatrix svd", os, __LINE__, dPessimisticSp);
        CheckReal (cvm::_abs((~V)(2, 4)), cvm::_abs((treal)-4.283360651798634e-001), "rmatrix svd", os, __LINE__, dPessimisticSp);
        CheckReal (cvm::_abs((~V)(3, 4)), cvm::_abs((treal)4.851284633245337e-001),"rmatrix svd", os, __LINE__, dPessimisticSp);
        CheckReal (cvm::_abs((~V)(4, 4)), cvm::_abs((treal)5.260662365874236e-001),"rmatrix svd", os, __LINE__, dPessimisticSp);

        rmatrix rm6(3,4);
        for (j = 1; j <= 4; j++)
        {
            for (i = 1; i <= 3; i++)
            {
                rm6(i,j)  = - (treal) ((j - 1) * 4 + i);
            }
        }

        rv1.resize (3);
        rv2.resize (3);
        rmU.resize(3);
        rmVH.resize(4);
        rv1 = rm6.svd (rmU, rmVH);
        rv2.svd (rm6, rmU, rmVH);
        CheckBool    (rv1 == rv2,  true,  "srmatrix svd", os, __LINE__);

        singValM << rm6;
        singValM.set(0.);
        singValM(1,1) = rv1(1);
        singValM(2,2) = rv1(2);
        singValM(3,3) = rv1(3);
        CheckReal    ((rm6 * ~rmVH - rmU * singValM).norm(),  (treal) 0.,  "rmatrix svd", os, __LINE__, dPessimisticSp);
        CheckReal    ((~rm6 * rmU - ~(singValM * rmVH)).norm(),  (treal) 0.,  "rmatrix svd", os, __LINE__, dPessimisticSp);


        rv1.resize (4);
        rv2.resize (4);
        scmatrix cmU(4), cmVH(4);
        rv1 = scm1.svd (cmU, cmVH);
        rv2.svd (scm1, cmU, cmVH);
        CheckBool    (rv1 == rv2,  true,  "scmatrix svd", os, __LINE__);
        cv1 << cvector (rv1);
        scm << scmatrix (cv1);
        CheckReal    ((scm1 * ~cmVH - cmU * scm).norm(),  (treal) 0.,  "scmatrix svd", os, __LINE__, dPessimisticSp);
        CheckReal    ((~scm1 * cmU - ~(scm * cmVH)).norm(),   (treal) 0.,  "scmatrix svd", os, __LINE__, dPessimisticSp);

        scbm1(4,3)=-cr1;
        rv1 = scbm1.svd (cmU, cmVH);
        rv2.svd (scbm1);
        CheckReal    ((rv1 - rv2).norm(),  (treal) 0.,  "scbmatrix svd", os, __LINE__, dVeryPessimisticSp);
        rv2.svd (scbm1, cmU, cmVH);
        scm1 << scmatrix (srmatrix(rv1));
        CheckReal    ((scbm1 * ~cmVH - cmU * scm1).norm(),  (treal) 0.,  "scbmatrix svd", os, __LINE__, dPessimisticSp);
        CheckReal    ((~scbm1 * cmU - ~(scm1 * cmVH)).norm(),  (treal) 0.,  "scbmatrix svd", os, __LINE__, dPessimisticSp);

        srm4 *=     (treal) -1.;

        srm4(3,3) = (treal) 1.;

        srm4(4,4) = (treal) 1.;
        rv.resize (4);
        srm4.normalize();
        rv.solve (srm4, rv1);
        CheckReal    ((srm4 * rv - rv1).norm(),   0,  "srmatrix solve", os, __LINE__, dPessimisticSp);
        rv = srm4.solve (rv1);
        CheckReal    ((srm4 * rv - rv1).norm(),   0,  "srmatrix solve", os, __LINE__, dPessimisticSp);
        rv.solve (srbm2, rv1);
        CheckReal    ((srbm2 * rv - rv1).norm(),   0,  "srbmatrix solve", os, __LINE__, dPessimisticSp);
        rv = srbm2.solve (rv1);
        CheckReal    ((srbm2 * rv - rv1).norm(),   0,  "srbmatrix solve", os, __LINE__, dPessimisticSp);
        scm1.assign(c2);
        scm1(3,3) = tcomplex (1, -1);
        scm1(4,4) = tcomplex (1, -1);
        cv.resize (4);
        cv1.resize (4);
        cv.solve (scm1, cv1);
        CheckReal    ((scm1 * cv - cv1).norm(),   0,  "scmatrix solve", os, __LINE__, dPessimisticSp);
        cv.solve (scbm1, cv1);
        CheckReal    ((scbm1 * cv - cv1).norm(),   0,  "scbmatrix solve", os, __LINE__, dPessimisticSp);
        cv = scbm1.solve (cv1);
        CheckReal    ((scbm1 * cv - cv1).norm(),   0,  "scbmatrix solve", os, __LINE__, dPessimisticSp);



        srm.resize  (3);
        scm.resize  (3);
        cv .resize  (3);
        cv1.resize  (3);

        srm(1,1) = (treal) 0.1;  srm(1,2) = (treal) 0.2;  srm(1,3) = (treal) 0.1;
        srm(2,1) = (treal) 0.11; srm(2,2) = (treal) -2.9; srm(2,3) = (treal) -8.4;
        srm(3,1) = (treal) 0.;   srm(3,2) = (treal) 2.91; srm(3,3) = (treal) 8.2;
        cv.eig (srm, scm);
        cv1 = srm.eig (scm);
        CheckBool    (cv == cv1,   true,  "srmatrix eig", os, __LINE__);
        CheckReal    ((scmatrix (srm) * scm(1) - scm(1) * cv1(1)).norm(),   0,  "srmatrix eig", os, __LINE__, dPessimisticSp);
        CheckReal    ((scmatrix (srm) * scm(2) - scm(2) * cv1(2)).norm(),   0,  "srmatrix eig", os, __LINE__, dPessimisticSp);
        CheckReal    ((scmatrix (srm) * scm(3) - scm(3) * cv1(3)).norm(),   0,  "srmatrix eig", os, __LINE__, dPessimisticSp);

        cv.eig (srm, scm, false);
        cv1 = srm.eig (scm, false);
        CheckBool    (cv == cv1,   true,  "srmatrix eig, left", os, __LINE__);

        CheckReal    ((~scm(1) * scmatrix (srm) - ~scm(1) * cv1(1)).norm(),   0,  "srmatrix eig, left", os, __LINE__, dPessimisticSp);
        CheckReal    ((~scm(2) * scmatrix (srm) - ~scm(2) * cv1(2)).norm(),   0,  "srmatrix eig, left", os, __LINE__, dPessimisticSp);
        CheckReal    ((~scm(3) * scmatrix (srm) - ~scm(3) * cv1(3)).norm(),   0,  "srmatrix eig, left", os, __LINE__, dPessimisticSp);

        srm(2,2) = (treal) 2.9;
        cv.eig (srm, scm);
        cv1 = srm.eig (scm);
        CheckBool    (cv == cv1,   true,  "srmatrix eig", os, __LINE__);
        CheckReal    ((scmatrix (srm) * scm(1) - scm(1) * cv1(1)).norm(),   0,  "srmatrix eig", os, __LINE__, dPessimisticSp);
        CheckReal    ((scmatrix (srm) * scm(2) - scm(2) * cv1(2)).norm(),   0,  "srmatrix eig", os, __LINE__, dPessimisticSp);
        CheckReal    ((scmatrix (srm) * scm(3) - scm(3) * cv1(3)).norm(),   0,  "srmatrix eig", os, __LINE__, dPessimisticSp);
        cv.eig (srm, scm, false);
        cv1 = srm.eig (scm, false);



        CheckBool    (cv == cv1,   true,  "srmatrix eig, left", os, __LINE__);
        CheckReal    ((~scm(1) * scmatrix (srm) - ~scm(1) * cv1(1)).norm(),   0,  "srmatrix eig, left", os, __LINE__, dPessimisticSp);
        CheckReal    ((~scm(2) * scmatrix (srm) - ~scm(2) * cv1(2)).norm(),   0,  "srmatrix eig, left", os, __LINE__, dPessimisticSp);
        CheckReal    ((~scm(3) * scmatrix (srm) - ~scm(3) * cv1(3)).norm(),   0,  "srmatrix eig, left", os, __LINE__, dPessimisticSp);

        scm1.resize  (3);
        cv.eig (scm1, scm);
        cv1 = scm1.eig (scm);
        CheckBool    (cv == cv1,   true,  "scmatrix eig", os, __LINE__);
        CheckReal    ((scm1 * scm(1) - scm(1) * cv1(1)).norm(),   0,  "scmatrix eig", os, __LINE__, dPessimisticSp);
        CheckReal    ((scm1 * scm(2) - scm(2) * cv1(2)).norm(),   0,  "scmatrix eig", os, __LINE__, dPessimisticSp);
        CheckReal    ((scm1 * scm(3) - scm(3) * cv1(3)).norm(),   0,  "scmatrix eig", os, __LINE__, dPessimisticSp);

        cv.eig (scm1, scm, false);
        cv1 = scm1.eig (scm , false);

        CheckReal    ((cv - cv1).norm(),  0.,  "scmatrix eig, left", os, __LINE__, dPessimisticSp);
        CheckReal    ((~scm(1) * scm1 - ~scm(1) * cv1(1)).norm(),   0,  "scmatrix eig, left", os, __LINE__, dPessimisticSp);
        CheckReal    ((~scm(2) * scm1 - ~scm(2) * cv1(2)).norm(),   0,  "scmatrix eig, left", os, __LINE__, dPessimisticSp);
        CheckReal    ((~scm(3) * scm1 - ~scm(3) * cv1(3)).norm(),   0,  "scmatrix eig, left", os, __LINE__, dPessimisticSp);

        scm.resize  (4);
        cv .resize  (4);
        cv1.resize  (4);
        cv.eig (srbm2, scm);
        cv1 = srbm2.eig (scm);
        CheckBool    (cv == cv1,   true,  "srbmatrix eig", os, __LINE__);
        CheckReal    ((scmatrix (srbm2) * scm(1) - scm(1) * cv1(1)).norm(),   0,  "srbmatrix eig", os, __LINE__, dPessimisticSp);
        CheckReal    ((scmatrix (srbm2) * scm(2) - scm(2) * cv1(2)).norm(),   0,  "srbmatrix eig", os, __LINE__, dPessimisticSp);
        CheckReal    ((scmatrix (srbm2) * scm(3) - scm(3) * cv1(3)).norm(),   0,  "srbmatrix eig", os, __LINE__, dPessimisticSp);
        CheckReal    ((scmatrix (srbm2) * scm(4) - scm(4) * cv1(4)).norm(),   0,  "srbmatrix eig", os, __LINE__, dPessimisticSp);

        cv.eig (scbm2, scm);
        cv1 = scbm2.eig (scm);
        CheckBool    (cv == cv1,   true,  "scbmatrix eig", os, __LINE__);
        CheckReal    ((scbm2 * scm(1) - scm(1) * cv1(1)).norm(),   0,  "scbmatrix eig", os, __LINE__, dPessimisticSp);
        CheckReal    ((scbm2 * scm(2) - scm(2) * cv1(2)).norm(),   0,  "scbmatrix eig", os, __LINE__, dPessimisticSp);
        CheckReal    ((scbm2 * scm(3) - scm(3) * cv1(3)).norm(),   0,  "scbmatrix eig", os, __LINE__, dPessimisticSp);
        CheckReal    ((scbm2 * scm(4) - scm(4) * cv1(4)).norm(),   0,  "scbmatrix eig", os, __LINE__, dPessimisticSp);


        rvector b(4), x(4);
        srsmatrix B(4);
        srmatrix EV(4);

        B.set(1,1,(treal)1.00000000000000e+000);
        B.set(2,1,(treal)5.55244534996568e-001); B.set(2,2,(treal)2.00000000000000e+000);
        B.set(3,1,(treal)1.00000000000000e+003); B.set(3,2,(treal)1.38811133749142e+000); B.set(3,3,(treal)3.00000000000000e+000);
        B.set(4,1,(treal)1.94335587248799e+000); B.set(4,2,(treal)2.22097813998627e+000); B.set(4,3,(treal)2.49860040748456e+000); B.set(4,4,(treal)4.00000000000000e+000);

        b.eig (B, EV);
        x = B.eig ();

        CheckReal    ((x - b).norm(),   0.,  "srsmatrix eig", os, __LINE__, dPessimisticSp);
        CheckReal    ((B * EV(1) - EV(1) * x(1)).norm(),   0,  "srsmatrix eig", os, __LINE__, dVeryPessimisticSp);
        CheckReal    ((B * EV(2) - EV(2) * x(2)).norm(),   0,  "srsmatrix eig", os, __LINE__, dVeryPessimisticSp);
        CheckReal    ((B * EV(3) - EV(3) * x(3)).norm(),   0,  "srsmatrix eig", os, __LINE__, dVeryPessimisticSp);
        // put schmatrix here

        rv1 = cv1.real();
        rv2 = cv1.imag();
        CheckReal    (rv1[4] - cv1(4).real(),   0,  "cvector::real", os, __LINE__);
        CheckReal    (rv2[3] - cv1(3).imag(),   0,  "cvector::imag", os, __LINE__);

        rm1 = cm1.real();
        rm2 << cm1.imag();
        CheckReal    (rm1[3][2] - cm1(3,2).real(),   0,  "cmatrix::real", os, __LINE__);
        CheckReal    (rm2(2, 3) - cm1(2,3).imag(),   0,  "cmatrix::imag", os, __LINE__);

        srm = scm1.real();
        srm1 << scm1.imag();
        CheckReal    (srm[3][2]  - scm1(3,2).real(),   0,  "scmatrix::real", os, __LINE__);
        CheckReal    (srm1(2, 3) - scm1(2,3).imag(),   0,  "scmatrix::imag", os, __LINE__);

        cv1.set(tcomplex ((treal)-13.45, (treal)1.778));
        cr1 = cv1[2];
        cv << cv1.conj();
        CheckComplex (cv1(2),   conj(cr1),  "cvector::conj", os, __LINE__);
        cv.conj (cv1);
        CheckComplex (cv(2),   cr1,  "cvector::conj", os, __LINE__);
        cv1 = ~cv;
        CheckComplex (cv1(2),   conj(cr1),  "cvector ~", os, __LINE__);

        cr1 = cm1(2,3);
        cm << cm1.conj();
        CheckComplex (cm1(3,2),   conj(cr1),  "cmatrix::conj", os, __LINE__);

        cm.conj (cm1);
        CheckComplex (cm(2,3),   cr1,  "cmatrix::conj", os, __LINE__);
        cm.resize(2,3);
        cm1.resize(3,2);
        cm1 = ~cm;
        CheckComplex (cm1(3,2),   conj(cr1),  "cmatrix ~", os, __LINE__);

        cr1 = scm1(2,3);
        scm << scm1.conj();
        CheckComplex (scm1(3,2),   conj(cr1),  "scmatrix::conj", os, __LINE__);
        scm.conj (scm1);
        CheckComplex (scm(2,3),   cr1,  "scmatrix::conj", os, __LINE__);
        scm1 = ~scm;
        CheckComplex (scm1(3,2),   conj(cr1),  "scmatrix ~", os, __LINE__);

        cr1 = scbm1(2,3);
        scbm << scbm1.conj();
        CheckComplex (scbm1(3,2),   conj(cr1),  "scbmatrix::conj", os, __LINE__);
        scbm.conj (scbm1);
        CheckComplex (scbm(2,3),   cr1,  "scbmatrix::conj", os, __LINE__);
        scbm1 = ~scbm;
        CheckComplex (scbm1(3,2),   conj(cr1),  "scbmatrix ~", os, __LINE__);


        r1 = (treal) 1.389;
        rv1.set(r1);
        cr1 = cv1[3];
        cv1.assign_real(rv1);
        CheckReal    (cv1(3).real(),   r1,  "cvector::assign_real", os, __LINE__);
        CheckReal    (cv1(3).imag(),   cr1.imag(),  "cvector::assign_real", os, __LINE__);
        cv1.assign_imag(rv1);
        CheckReal    (cv1(2).real(),   r1,  "cvector::assign_imag", os, __LINE__);
        CheckReal    (cv1(2).imag(),   r1,  "cvector::assign_imag", os, __LINE__);

        rm1.resize(3,2);
        rm1.set(r1);
        cr1 = cm1(3,2);
        cm1.assign_real(rm1);
        CheckReal    (cm1(3,2).real(),   r1,  "cmatrix::assign_real", os, __LINE__);
        CheckReal    (cm1(3,2).imag(),   cr1.imag(),  "cmatrix::assign_real", os, __LINE__);
        cm1.assign_imag(rm1);
        CheckReal    (cm1(2,2).real(),   r1,  "cmatrix::assign_imag", os, __LINE__);
        CheckReal    (cm1(2,2).imag(),   r1,  "cmatrix::assign_imag", os, __LINE__);


        srm.set(r1);
        cr1 = scm1(3,2);
        scm1.assign_real(srm);
        CheckReal    (scm1(3,2).real(),   r1,  "scmatrix::assign_real", os, __LINE__);
        CheckReal    (scm1(3,2).imag(),   cr1.imag(),  "scmatrix::assign_real", os, __LINE__);
        scm1.assign_imag(srm);
        CheckReal    (scm1(2,2).real(),   r1,  "scmatrix::assign_imag", os, __LINE__);
        CheckReal    (scm1(2,2).imag(),   r1,  "scmatrix::assign_imag", os, __LINE__);

        srbm.resize(4);
        srbm.set(r1);
        scbm1.resize_lu (1, 2);
        cr1 = scbm1(1, 2);
        scbm1.assign_real(srbm);
        CheckReal    (scbm1(1,2).real(),   r1,  "scbmatrix::assign_real", os, __LINE__);
        CheckReal    (scbm1(1,2).imag(),   cr1.imag(),  "scbmatrix::assign_real", os, __LINE__);
        scbm1.assign_imag(srbm);
        CheckReal    (scbm1(1,2).real(),   r1,  "scbmatrix::assign_imag", os, __LINE__);

        CheckReal    (scbm1(1,2).imag(),   r1,  "scbmatrix::assign_imag", os, __LINE__);

        rm.set(1.);
        rm.normalize();
        rm(2,3) = (treal) -1.1;
        rm(1,2) = (treal) 1.e-9;
        CheckInt     (rm.rowofmax (),   2,  "rmatrix::rowofmax", os, __LINE__);

        CheckInt     (rm.colofmax (),   3,  "rmatrix::colofmax", os, __LINE__);
        CheckInt     (rm.rowofmin (),   1,  "rmatrix::rowofmin", os, __LINE__);
        CheckInt     (rm.colofmin (),   2,  "rmatrix::colofmin", os, __LINE__);
        CheckInt     (rm.msize (),      2,  "rmatrix::msize", os, __LINE__);
        CheckInt     (rm.nsize (),      3,  "rmatrix::nsize", os, __LINE__);
        rm1 << rm.swap_rows (1,2);
        CheckReal    (rm1(1,3),   (treal) -1.1,  "rmatrix::swap_rows", os, __LINE__);
        rm1.swap_cols (1,2);
        CheckReal    (rm1(2,1),   (treal) 1.e-9,  "rmatrix::swap_cols", os, __LINE__);

        srm.set(1.);
        srm.normalize();
        srm(2,3) = (treal) -1.1;
        srm(1,2) = (treal) 1.e-9;
        CheckInt     (srm.rowofmax (),   2,  "srmatrix::rowofmax", os, __LINE__);
        CheckInt     (srm.colofmax (),   3,  "srmatrix::colofmax", os, __LINE__);
        CheckInt     (srm.rowofmin (),   1,  "srmatrix::rowofmin", os, __LINE__);
        CheckInt     (srm.colofmin (),   2,  "srmatrix::colofmin", os, __LINE__);
        CheckInt     (srm.msize (),      3,  "srmatrix::msize", os, __LINE__);
        CheckInt     (srm.nsize (),      3,  "srmatrix::nsize", os, __LINE__);
        srm1 << srm.swap_rows (1,2);
        CheckReal    (srm1(1,3),   (treal) -1.1,  "srmatrix::swap_rows", os, __LINE__);
        srm1.swap_cols (1,2);
        CheckReal    (srm1(2,1),   (treal) 1.e-9,  "srmatrix::swap_cols", os, __LINE__);

        cm.set(tcomplex ((treal) 1., treal (1.)));
        cm.normalize();
        cm(2,3) = tcomplex ((treal) 1.1, treal (1.1));
        cm(1,2) = tcomplex ((treal) 1.e-9, (treal) 0.);
        CheckInt     (cm.rowofmax (),   2,  "cmatrix::rowofmax", os, __LINE__);
        CheckInt     (cm.colofmax (),   3,  "cmatrix::colofmax", os, __LINE__);
        CheckInt     (cm.rowofmin (),   1,  "cmatrix::rowofmin", os, __LINE__);
        CheckInt     (cm.colofmin (),   2,  "cmatrix::colofmin", os, __LINE__);
        CheckInt     (cm.msize (),      2,  "cmatrix::msize", os, __LINE__);
        CheckInt     (cm.nsize (),      3,  "cmatrix::nsize", os, __LINE__);
        cm1 << cm.swap_rows (1,2);
        CheckComplex (cm1(1,3),   tcomplex ((treal) 1.1, treal (1.1)),  "cmatrix::swap_rows", os, __LINE__);
        cm1.swap_cols (1,2);
        CheckComplex (cm1(2,1),   tcomplex ((treal) 1.e-9, (treal) 0.),  "cmatrix::swap_cols", os, __LINE__);

        scm.set(tcomplex ((treal) 1., treal (1.)));
        scm.normalize();
        scm(2,3) = tcomplex ((treal) 1.1, treal (1.1));
        scm(1,2) = tcomplex ((treal) 1.e-9, (treal) 0.);
        CheckInt     (scm.rowofmax (),   2,  "scmatrix::rowofmax", os, __LINE__);
        CheckInt     (scm.colofmax (),   3,  "scmatrix::colofmax", os, __LINE__);
        CheckInt     (scm.rowofmin (),   1,  "scmatrix::rowofmin", os, __LINE__);
        CheckInt     (scm.colofmin (),   2,  "scmatrix::colofmin", os, __LINE__);
        CheckInt     (scm.msize (),      3,  "scmatrix::msize", os, __LINE__);
        CheckInt     (scm.nsize (),      3,  "scmatrix::nsize", os, __LINE__);
        scm1 << scm.swap_rows (1,2);
        CheckComplex (scm1(1,3),   tcomplex ((treal) 1.1, treal (1.1)),  "scmatrix::swap_rows", os, __LINE__);
        scm1.swap_cols (1,2);
        CheckComplex (scm1(2,1),   tcomplex ((treal) 1.e-9, (treal) 0.),  "scmatrix::swap_cols", os, __LINE__);

        srbm.diag(0) = rv;
        srbm.normalize();
        srbm(2,3) = (treal) -1.1;
        srbm(1,2) = (treal) 1.e-7;

        CheckInt     (srbm.rowofmax (),   2,  "srbmatrix::rowofmax", os, __LINE__);
        CheckInt     (srbm.colofmax (),   3,  "srbmatrix::colofmax", os, __LINE__);
        CheckInt     (srbm.rowofmin (),   3,  "srbmatrix::rowofmin", os, __LINE__);
        CheckInt     (srbm.colofmin (),   1,  "srbmatrix::colofmin", os, __LINE__);
        CheckInt     (srbm.msize (),      4,  "srbmatrix::msize", os, __LINE__);
        CheckInt     (srbm.nsize (),      4,  "srbmatrix::nsize", os, __LINE__);

        scbm.diag(0).set(tcomplex ((treal) 1., treal (1.)));
        scbm.normalize();
        scbm(2,3) = tcomplex ((treal) 2., treal (1.));
        scbm(1,2) = tcomplex ((treal) -1.e-10, treal (-1.e-10));
        CheckInt     (scbm.rowofmax (),   2,  "scbmatrix::rowofmax", os, __LINE__);
        CheckInt     (scbm.colofmax (),   3,  "scbmatrix::colofmax", os, __LINE__);
        CheckInt     (scbm.rowofmin (),   4,  "scbmatrix::rowofmin", os, __LINE__);
        CheckInt     (scbm.colofmin (),   1,  "scbmatrix::colofmin", os, __LINE__);
        CheckInt     (scbm.msize (),      4,  "scbmatrix::msize", os, __LINE__);
        CheckInt     (scbm.nsize (),      4,  "scbmatrix::nsize", os, __LINE__);


        for (i = 0; i < 100; i++)
        {
            a1[i] = (treal) sqrt ((treal) (i + 1));
            a2[i] = (treal) (i + 1) / (treal) 10.;
            c1[i] = tcomplex(a1[i], a2[i]);
        }

        rm2.set((treal)-0.34);
        rm2(2,3) = (treal) 0.;

        CheckInt     (rm2.rank (),      2,  "rmatrix::rank", os, __LINE__);
        rm2.assign(a1);
        CheckInt     (rm2.rank (),      4,  "rmatrix::rank", os, __LINE__);

        srm2.assign(a1);
        srm2[2].set((treal)0.);
        CheckInt     (srm2.rank (),     2,  "srmatrix::rank", os, __LINE__);
        srm2.diag(0).set((treal)0.);
        CheckInt     (srm2.rank (),     2,  "srmatrix::rank", os, __LINE__);

        cm2.resize (3, 4);
        cm2.assign(c1);
        CheckInt     (cm2.rank (),      3,  "cmatrix::rank", os, __LINE__);

        scm2.assign(c1);
        CheckInt     (scm2.rank (),     3,  "scmatrix::rank", os, __LINE__);
        scm2.diag(0).set((treal)0.);
        CheckInt     (scm2.rank (),     3,  "scmatrix::rank", os, __LINE__);

        srbm.assign(a1);
        CheckInt     (srbm.rank (),     4,  "srbmatrix::rank", os, __LINE__);

        scbm.assign(c1);
        CheckInt     (scbm.rank (),     4,  "scbmatrix::rank", os, __LINE__);


        r1 = (treal) -8.76;
        srm2.set(r1);
        srm2.diag(1).set(0.);
        CheckReal    (srm2(1,1),   r1,  "srmatrix::diag", os, __LINE__);
        CheckReal    (srm2(1,2),   (treal) 0.,  "srmatrix::diag", os, __LINE__);

        cr1 = tcomplex ((treal) -8.76, (treal) -3.6);
        scm2.set(cr1);

        scm2.diag(1).set(0L);
        CheckComplex (scm2(1,1),   cr1,  "scmatrix::diag", os, __LINE__);
        CheckComplex (scm2(1,2),   0,    "scmatrix::diag", os, __LINE__);

        srbm.set(r1);
        srbm.diag(1).set(0.);
        CheckReal    (srbm(1,1),   r1,  "srbmatrix::diag", os, __LINE__);
        CheckReal    (srbm(1,2),   (treal) 0.,  "srbmatrix::diag", os, __LINE__);

        scbm.set(cr1);
        scbm.diag(1).set(0.);
        CheckComplex (scbm(1,1),   cr1,  "scbmatrix::diag", os, __LINE__);
        CheckComplex (scbm(1,2),   tcomplex ((treal) 0., (treal) 0.),  "srbmatrix::diag", os, __LINE__);


        srm2.set(r1);
        srm2++;
        CheckReal    (srm2(1,1),   r1 + 1,  "srmatrix++", os, __LINE__);
        ++srm2;
        CheckReal    (srm2(2,2),   r1 + 2,  "++srmatrix", os, __LINE__);
        srm2--;
        CheckReal    (srm2(1,1),   r1 + 1,  "srmatrix--", os, __LINE__);
        --srm2;
        CheckReal    (srm2(2,2),   r1,      "--srmatrix", os, __LINE__);

        scm2.set(cr1);
        scm2++;
        CheckComplex (scm2(1,1),   cr1 + tcomplex (1), "scmatrix++", os, __LINE__);
        ++scm2;
        CheckComplex (scm2(2,2),   cr1 + tcomplex (2), "++scmatrix", os, __LINE__);
        scm2--;
        CheckComplex (scm2(1,1),   cr1 + tcomplex (1), "scmatrix--", os, __LINE__);
        --scm2;
        CheckComplex (scm2(2,2),   cr1,                "--scmatrix", os, __LINE__);

        srbm.set(r1);
        srbm++;
        CheckReal    (srbm(1,1),   r1 + 1,  "srbmatrix++", os, __LINE__);
        ++srbm;
        CheckReal    (srbm(2,2),   r1 + 2,  "++srbmatrix", os, __LINE__);
        srbm--;
        CheckReal    (srbm(1,1),   r1 + 1,  "srbmatrix--", os, __LINE__);
        --srbm;
        CheckReal    (srbm(2,2),   r1,      "--srbmatrix", os, __LINE__);

        scbm.set(cr1);
        scbm++;
        CheckComplex (scbm(1,1),   cr1 + tcomplex (1), "scbmatrix++", os, __LINE__);
        ++scbm;
        CheckComplex (scbm(2,2),   cr1 + tcomplex (2), "++scbmatrix", os, __LINE__);
        scbm--;
        CheckComplex (scbm(1,1),   cr1 + tcomplex (1), "scbmatrix--", os, __LINE__);
        --scbm;
        CheckComplex (scbm(2,2),   cr1,                "--scbmatrix", os, __LINE__);


        srm << srm2.identity();
        CheckReal    (srm(1,1),    1,  "srmatrix::identity", os, __LINE__);
        CheckReal    (srm2(1,2),   0,  "srmatrix::identity", os, __LINE__);

        scm << scm2.identity();
        CheckComplex (scm(1,1),    1,  "scmatrix::identity", os, __LINE__);
        CheckComplex (scm2(1,2),   0,  "scmatrix::identity", os, __LINE__);


        srbm << srbm1.identity();
        CheckReal    (srbm(1,1),    1,  "srbmatrix::identity", os, __LINE__);
        CheckReal    (srbm1(1,2),   0,  "srbmatrix::identity", os, __LINE__);

        scbm << scbm1.identity();
        CheckComplex (scbm(1,1),    1,  "scbmatrix::identity", os, __LINE__);
        CheckComplex (scbm1(1,2),   0,  "scbmatrix::identity", os, __LINE__);


        rm.assign(a2);
        rm1 << rm.transpose();
        CheckReal    (rm(1,2),    a2[1],  "rmatrix::transpose", os, __LINE__);
        rm1.resize (rm.nsize(), rm.msize());
        rm1.transpose (rm);
        CheckReal    (rm1(1,2),   a2[2],  "rmatrix::transpose", os, __LINE__);

        srm.assign(a2);
        srm1 << srm.transpose();
        CheckReal    (srm(1,2),   a2[1],  "srmatrix::transpose", os, __LINE__);
        srm1.transpose (srm);
        CheckReal    (srm1(1,2),  a2[3],  "srmatrix::transpose", os, __LINE__);

        srbm.resize_lu (2, 1);
        srbm.assign(a2);
        srbm1 << srbm.transpose();
        CheckReal    (srbm(1,2),   a2[2],  "srbmatrix::transpose", os, __LINE__);
        srbm.transpose();
        srbm.assign(a2);
        srbm1.transpose (srbm);
        CheckReal    (srbm1(1,2),  a2[2],  "srbmatrix::transpose", os, __LINE__);

        scbm1 = scbm.assign(c1);
        scbm1.conj();
        CheckComplex (scbm1(1,2),  conj(scbm(2,1)),  "scbmatrix::conj", os, __LINE__);
        scbm1.assign(c1);
        scbm.conj(scbm1);
        CheckComplex (scbm1(1,3),  conj(scbm(3,1)),  "scbmatrix::conj", os, __LINE__);
        CheckComplex (scbm1(1,4),  conj(scbm(4,1)),  "scbmatrix::conj", os, __LINE__);


        srm.resize(3);
        srm1.resize(3);
        srm(1,1) = (treal) 2;    srm(1,2) = (treal) -0.8; srm(1,3) = (treal) -0.7;
        srm(2,1) = (treal) -0.4; srm(2,2) = (treal) -1;   srm(2,3) = (treal) -0.8;
        srm(3,1) = (treal) -0.6; srm(3,2) = (treal) -1.2; srm(3,3) = (treal) -0.9;

        srm1 = srm.exp();
        CheckReal    (srm1(1,1), (treal)  8.484495096274699e+000, "srmatrix::exp", os, __LINE__, dPessimisticSp);
        CheckReal    (srm1(1,2), (treal) -1.555963610758445e+000, "srmatrix::exp", os, __LINE__, dPessimisticSp);
        CheckReal    (srm1(1,3), (treal) -1.484978300761370e+000, "srmatrix::exp", os, __LINE__, dPessimisticSp);
        CheckReal    (srm1(2,1), (treal) -7.330690267073194e-001, "srmatrix::exp", os, __LINE__, dPessimisticSp);
        CheckReal    (srm1(2,2), (treal)  6.959256837027834e-001, "srmatrix::exp", os, __LINE__, dPessimisticSp);
        CheckReal    (srm1(2,3), (treal) -2.385221030493092e-001, "srmatrix::exp", os, __LINE__, dPessimisticSp);
        CheckReal    (srm1(3,1), (treal) -1.324167433420492e+000, "srmatrix::exp", os, __LINE__, dPessimisticSp);
        CheckReal    (srm1(3,2), (treal) -3.128703759020610e-001, "srmatrix::exp", os, __LINE__, dPessimisticSp);
        CheckReal    (srm1(3,3), (treal)  8.267946985957282e-001, "srmatrix::exp", os, __LINE__, dPessimisticSp);

        scm.resize(3);
        scm1.resize(3);
        scm(1,1) = tcomplex ((treal) 1.e-01, (treal) 2.e-01); 
        scm(1,2) = tcomplex ((treal) 3.e-01, (treal) 4.e-01);
        scm(1,3) = tcomplex ((treal) 5.e-01, (treal) 6.e-01);
        scm(2,1) = tcomplex ((treal) 1.e+00, (treal) 1.e+00);
        scm(2,2) = tcomplex ((treal) 0.e+00, (treal) 0.e+00);
        scm(2,3) = tcomplex ((treal) 0.e+00, (treal) 0.e+00);
        scm(3,1) = tcomplex ((treal) -1.e-01,(treal) -1.e-01);
        scm(3,2) = tcomplex ((treal) -3.e-01,(treal) -3.e-01);
        scm(3,3) = tcomplex ((treal) -5.e-01,(treal) -5.e-01);

        scm(2,2)=-cr1;
        scm1 = scm.exp();
        CheckComplex (scm1(1,1),  tcomplex ((treal) -4.816680814596321e+000, (treal) -4.855745768190474e+001), "scmatrix::exp", os, __LINE__, dPessimisticSp);
        CheckComplex (scm1(2,1),  tcomplex((treal) -5.878045515841980e+002, (treal) -7.809398663068483e+002), "scmatrix::exp", os, __LINE__, dPessimisticSp);
        CheckComplex (scm1(3,1),  tcomplex((treal) 1.111774160999558e+001, (treal)  3.979363145886382e+001), "scmatrix::exp", os, __LINE__, dPessimisticSp);
        CheckComplex (scm1(1,2),  tcomplex((treal) -1.623884970745376e+002, (treal)  -2.805917519984524e+002), "scmatrix::exp", os, __LINE__, dPessimisticSp);
        CheckComplex (scm1(2,2),  tcomplex((treal) -5.604582009475869e+003, (treal)   -3.219074690441815e+003), "scmatrix::exp", os, __LINE__, dPessimisticSp);
        CheckComplex (scm1(3,2),  tcomplex((treal) 1.715815440786858e+002, (treal)  2.129974004265882e+002), "scmatrix::exp", os, __LINE__, dPessimisticSp);
        CheckComplex (scm1(1,3),  tcomplex((treal) 1.710263520348249e+000, (treal)  -3.149555947204208e+000), "scmatrix::exp", os, __LINE__, dPessimisticSp);
        CheckComplex (scm1(2,3),  tcomplex((treal) -1.432034221529735e+001, (treal)  -7.375809596051487e+001), "scmatrix::exp", os, __LINE__, dPessimisticSp);
        CheckComplex (scm1(3,3),  tcomplex((treal) -4.639004479804901e-002, (treal)  2.814422951041492e+000), "scmatrix::exp", os, __LINE__, dPessimisticSp);


        srbm1.resize(2);
        srbm1.resize_lu(0, 1);
        srbm1(1,1) = (treal) 1.3;
        srbm1(1,2) = (treal) -11.2;
        srbm1(2,2) = (treal) 4.1;
        os << srbm1;


        srm1 << srbm1.exp();
        CheckReal    (srm1(1,1), (treal) 3.669296667619233e+000, "srbmatrix::exp", os, __LINE__, dPessimisticSp);
        CheckReal    (srm1(1,2), (treal) -2.266839637189685e+002, "srbmatrix::exp", os, __LINE__, dPessimisticSp);
        CheckReal    (srm1(2,1), (treal) 0, "srbmatrix::exp", os, __LINE__, dPessimisticSp);
        CheckReal    (srm1(2,2), (treal) 6.034028759736115e+001, "srbmatrix::exp", os, __LINE__, dPessimisticSp);

        iarray aPivots(3);
        srmatrix mLU (3), mLU2 (3), mLo(3), mUp(3);

        mLU = srm.low_up (aPivots);
        mLU2.low_up (srm, aPivots);
        CheckBool    (mLU == mLU2, true, "srmatrix::low_up", os, __LINE__);

        mLo.identity ();
        mLo(2,1) = mLU(2,1);
        mLo(3,1) = mLU(3,1);
        mLo(3,2) = mLU(3,2);

        mUp(1,1) = mLU(1,1);
        mUp(1,2) = mLU(1,2);
        mUp(1,3) = mLU(1,3);
        mUp(2,2) = mLU(2,2);
        mUp(2,3) = mLU(2,3);
        mUp(3,3) = mLU(3,3);

        mLU = mLo * mUp;
        for (l = 3; l >= 1; l--) {
            mLU.swap_rows (l, aPivots[l]);
        }
        CheckReal    ((srm - mLU).norminf(), (treal) 0, "srmatrix::low_up", os, __LINE__, dPessimisticSp);


        scmatrix cmLU (3), cmLU2 (3), cmLo(3), cmUp(3);
        cmLU = scm.low_up (aPivots);
        cmLU2.low_up (scm, aPivots);
        CheckBool    (cmLU == cmLU2, true, "scmatrix::low_up", os, __LINE__);

        cmLo.identity ();
        cmLo(2,1) = cmLU(2,1);
        cmLo(3,1) = cmLU(3,1);
        cmLo(3,2) = cmLU(3,2);

        cmUp(1,1) = cmLU(1,1);
        cmUp(1,2) = cmLU(1,2);
        cmUp(1,3) = cmLU(1,3);
        cmUp(2,2) = cmLU(2,2);

        cmUp(2,3) = cmLU(2,3);
        cmUp(3,3) = cmLU(3,3);

        cmLU = cmLo * cmUp;
        for (l = 3; l >= 1; l--) {
                cmLU.swap_rows (l, aPivots[l]);
        }
        CheckReal    ((scm - cmLU).norminf(), (treal) 0, "scmatrix::low_up", os, __LINE__, dPessimisticSp);


        srm1 << srm.inv();
        CheckReal    (((srm1 * srm)--).norm(), (treal) 0, "srmatrix::inv", os, __LINE__, dPessimisticSp);
        scm1 << scm.inv();
        CheckReal    (((scm1 * scm)--).norm(), (treal) 0, "scmatrix::inv", os, __LINE__, dPessimisticSp);


        rv.resize(11);
        rv(1)  = (treal) 2.2;
        rv(2)  = (treal) 1.3;
        rv(3)  = (treal) 1.1;
        rv(4)  = (treal) - 0.9;
        rv(5)  = (treal) 0.2;
        rv(6)  = (treal) - 0.45;
        rv(7)  = (treal) 45;
        rv(8)  = (treal) - 30;
        rv(9)  = (treal) 10;
        rv(10) = (treal) 3;
        rv(11) = (treal) 3.2;

        srm1.polynom (srm, rv);
        srm2 << srm.polynom (rv);
        CheckBool    (srm1 == srm2, true, "srmatrix::polynom", os, __LINE__);
        CheckReal    (srm1(1,1), (treal)  1.415106245372072e+004, "srmatrix::polynom", os, __LINE__, dPessimisticSp);
        CheckReal    (srm1(1,2), (treal)  8.018578436580816e+002, "srmatrix::polynom", os, __LINE__, dPessimisticSp);
        CheckReal    (srm1(1,3), (treal)  1.516628273102821e+002, "srmatrix::polynom", os, __LINE__, dPessimisticSp);
        CheckReal    (srm1(2,1), (treal)  6.009153894255998e+002, "srmatrix::polynom", os, __LINE__, dPessimisticSp);
        CheckReal    (srm1(2,2), (treal)  8.458618026988163e+003, "srmatrix::polynom", os, __LINE__, dPessimisticSp);
        CheckReal    (srm1(2,3), (treal)  6.668127559823842e+003, "srmatrix::polynom", os, __LINE__, dPessimisticSp);
        CheckReal    (srm1(3,1), (treal) -9.855925384439991e+001, "srmatrix::polynom", os, __LINE__, dPessimisticSp);
        CheckReal    (srm1(3,2), (treal)  1.020217780733232e+004, "srmatrix::polynom", os, __LINE__, dPessimisticSp);
        CheckReal    (srm1(3,3), (treal)  8.075071634102441e+003, "srmatrix::polynom", os, __LINE__, dPessimisticSp);

        rv.resize(3);
        srm1.polynom (srm, rv);
        CheckReal    (srm1(1,1), (treal)  1.001400000000000e+001, "srmatrix::polynom", os, __LINE__, dPessimisticSp);
        CheckReal    (srm1(1,2), (treal) -9.960000000000001e-001, "srmatrix::polynom", os, __LINE__, dPessimisticSp);
        CheckReal    (srm1(1,3), (treal) -1.053000000000000e+000, "srmatrix::polynom", os, __LINE__, dPessimisticSp);
        CheckReal    (srm1(2,1), (treal) -4.320000000000001e-001, "srmatrix::polynom", os, __LINE__, dPessimisticSp);
        CheckReal    (srm1(2,2), (treal)  3.408000000000000e+000, "srmatrix::polynom", os, __LINE__, dPessimisticSp);
        CheckReal    (srm1(2,3), (treal)  9.400000000000004e-001, "srmatrix::polynom", os, __LINE__, dPessimisticSp);
        CheckReal    (srm1(3,1), (treal) -9.780000000000000e-001, "srmatrix::polynom", os, __LINE__, dPessimisticSp);
        CheckReal    (srm1(3,2), (treal)  1.476000000000000e+000, "srmatrix::polynom", os, __LINE__, dPessimisticSp);
        CheckReal    (srm1(3,3), (treal)  3.439000000000000e+000, "srmatrix::polynom", os, __LINE__, dPessimisticSp);


        cv.resize(11);
        cv(1)  = tcomplex ((treal) 2.2, (treal) -1);
        cv(2)  = tcomplex ((treal) 1.3, (treal) -0.6);
        cv(3)  = tcomplex ((treal) 1.1, (treal) 2.3);
        cv(4)  = tcomplex ((treal) -0.9);
        cv(5)  = tcomplex ((treal) 0.2, (treal) 1);
        cv(6)  = tcomplex ((treal) -0.45, (treal) 2);
        cv(7)  = tcomplex ((treal) 45, (treal) -17.3);
        cv(8)  = tcomplex ((treal) -30);
        cv(9)  = tcomplex ((treal) 10, (treal) 1.5);
        cv(10) = tcomplex ((treal) 3);
        cv(11) = tcomplex ((treal) 3.2, (treal) -18.9);


        scm(1,1) = tcomplex ((treal) 2., (treal) 0.1); 
        scm(1,2) = tcomplex ((treal) -8.e-001, (treal) -1);
        scm(1,3) = tcomplex ((treal) -7.e-001, (treal) -2.1);
        scm(2,1) = tcomplex ((treal) -4.e-001, (treal) -0.1);
        scm(2,2) = tcomplex ((treal) -1.e+000);
        scm(2,3) = tcomplex ((treal) -8.e-001, (treal) 3.1);
        scm(3,1) = tcomplex ((treal) -6.e-001,(treal) 1);
        scm(3,2) = tcomplex ((treal) -1.2e+000);
        scm(3,3) = tcomplex ((treal) -9.e-001,(treal) -5.4);

        scm1.polynom (scm, cv);

        CheckComplex (scm1(1,1),  tcomplex ((treal)  5.264016832618990e+006, (treal) -1.051212804982833e+007), "scmatrix::polynom", os, __LINE__, dPessimisticSp);
        CheckComplex (scm1(1,2),  tcomplex ((treal)  9.386518437203571e+006, (treal) -9.534002545240149e+006), "scmatrix::polynom", os, __LINE__, dPessimisticSp);
        CheckComplex (scm1(1,3),  tcomplex ((treal)  2.313187132312614e+007, (treal)  4.742508767071142e+007), "scmatrix::polynom", os, __LINE__, dPessimisticSp);
        CheckComplex (scm1(2,1),  tcomplex ((treal) -1.143556158726668e+007, (treal)  2.626370923270145e+007), "scmatrix::polynom", os, __LINE__, dPessimisticSp);
        CheckComplex (scm1(2,2),  tcomplex ((treal) -2.183671220461629e+007, (treal)  2.471364343201455e+007), "scmatrix::polynom", os, __LINE__, dPessimisticSp);
        CheckComplex (scm1(2,3),  tcomplex ((treal) -6.325599106881835e+007, (treal) -1.133746860502928e+008), "scmatrix::polynom", os, __LINE__, dPessimisticSp);
        CheckComplex (scm1(3,1),  tcomplex ((treal)  1.143469364494270e+007, (treal) -4.448575764049879e+007), "scmatrix::polynom", os, __LINE__, dPessimisticSp);
        CheckComplex (scm1(3,2),  tcomplex ((treal)  2.832544852276585e+007, (treal) -4.473797233313387e+007), "scmatrix::polynom", os, __LINE__, dPessimisticSp);
        CheckComplex (scm1(3,3),  tcomplex ((treal)  1.291773725514465e+008, (treal)  1.634454865648127e+008), "scmatrix::polynom", os, __LINE__, dPessimisticSp);

        scm2 << scm.polynom (cv);
        CheckReal    ((scm1.normalize() - scm2.normalize()).norm(), (treal) 0., "scmatrix::polynom", os, __LINE__, dVeryPessimisticSp);

        srbm.resize_lu(1,0);
        bool bThrew = false;
        try
        {
            srbm(1,2) = (treal) 1.;
        }
        catch (cvmexception e)
        {
            if (e.cause() == CVM_READ_ONLY_ACCESS) bThrew = true;
        }
        CheckBool (bThrew, true, "srbmatrix read only exception", os, __LINE__);

        srbm.diag(0).set(1.);
        srbm.diag(-1).set(1.);

        srm << srbm.exp();

        CheckReal    (srm(1,1), (treal)  2.718281828459041e+000, "srbmatrix::exp", os, __LINE__, dPessimisticSp);
        CheckReal    (srm(1,2), (treal)  0., "srbmatrix::exp", os, __LINE__, dPessimisticSp);
        CheckReal    (srm(1,3), (treal)  0., "srbmatrix::exp", os, __LINE__, dPessimisticSp);
        CheckReal    (srm(1,4), (treal)  0., "srbmatrix::exp", os, __LINE__, dPessimisticSp);

        CheckReal    (srm(2,1), (treal)  2.718281828459041e+000, "srbmatrix::exp", os, __LINE__, dPessimisticSp);
        CheckReal    (srm(2,2), (treal)  2.718281828459041e+000, "srbmatrix::exp", os, __LINE__, dPessimisticSp);
        CheckReal    (srm(2,3), (treal)  0., "srbmatrix::exp", os, __LINE__, dPessimisticSp);
        CheckReal    (srm(2,4), (treal)  0., "srbmatrix::exp", os, __LINE__, dPessimisticSp);

        CheckReal    (srm(3,1), (treal)  1.359140914229521e+000, "srbmatrix::exp", os, __LINE__, dPessimisticSp);
        CheckReal    (srm(3,2), (treal)  2.718281828459041e+000, "srbmatrix::exp", os, __LINE__, dPessimisticSp);
        CheckReal    (srm(3,3), (treal)  2.718281828459041e+000, "srbmatrix::exp", os, __LINE__, dPessimisticSp);
        CheckReal    (srm(3,4), (treal)  0., "srbmatrix::exp", os, __LINE__, dPessimisticSp);

        CheckReal    (srm(4,1), (treal)  4.530469714098402e-001, "srbmatrix::exp", os, __LINE__, dPessimisticSp);
        CheckReal    (srm(4,2), (treal)  1.359140914229521e+000, "srbmatrix::exp", os, __LINE__, dPessimisticSp);
        CheckReal    (srm(4,3), (treal)  2.718281828459041e+000, "srbmatrix::exp", os, __LINE__, dPessimisticSp);
        CheckReal    (srm(4,4), (treal)  2.718281828459041e+000, "srbmatrix::exp", os, __LINE__, dPessimisticSp);

        srm -= (srmatrix) srbm;

        scbm.resize_lu(1,0);
        bThrew = false;
        try
        {
            scbm(1,2) = (treal) 1.;
        }
        catch (cvmexception e)
        {
            if (e.cause() == CVM_READ_ONLY_ACCESS) bThrew = true;
        }
        CheckBool (bThrew, true, "scbmatrix read only exception", os, __LINE__);

        scbm.diag(0).set(tcomplex ((treal) 1., (treal) 1.));
        scbm.diag(-1).set(tcomplex ((treal) 1., (treal) 1.));

        scm << scbm.exp();

        CheckComplex (scm(1,1),  tcomplex ((treal)  1.468693939915887e+000, (treal) 2.287355287178844e+000), "scbmatrix::exp", os, __LINE__, dPessimisticSp);
        CheckComplex (scm(2,1),  tcomplex ((treal)  -8.186613472629570e-001, (treal) 3.756049227094730e+000), "scbmatrix::exp", os, __LINE__, dPessimisticSp);
        CheckComplex (scm(3,1),  tcomplex ((treal)  -2.287355287178843e+000, (treal) 1.468693939915886e+000), "scbmatrix::exp", os, __LINE__, dPessimisticSp);
        CheckComplex (scm(4,1),  tcomplex ((treal)  -1.252016409031576e+000, (treal) -2.728871157543187e-001), "scbmatrix::exp", os, __LINE__, dPessimisticSp);

        CheckComplex (scm(1,2),  tcomplex ((treal)  0., (treal) 0.), "scbmatrix::exp", os, __LINE__, dPessimisticSp);
        CheckComplex (scm(2,2),  tcomplex ((treal)  1.468693939915887e+000, (treal) 2.287355287178844e+000), "scbmatrix::exp", os, __LINE__, dPessimisticSp);
        CheckComplex (scm(3,2),  tcomplex ((treal)  -8.186613472629570e-001, (treal) 3.756049227094730e+000), "scbmatrix::exp", os, __LINE__, dPessimisticSp);
        CheckComplex (scm(4,2),  tcomplex ((treal)  -2.287355287178843e+000, (treal) 1.468693939915886e+000), "scbmatrix::exp", os, __LINE__, dPessimisticSp);

        CheckComplex (scm(1,3),  tcomplex ((treal)  0., (treal) 0.), "scbmatrix::exp", os, __LINE__, dPessimisticSp);
        CheckComplex (scm(2,3),  tcomplex ((treal)  0., (treal) 0.), "scbmatrix::exp", os, __LINE__, dPessimisticSp);
        CheckComplex (scm(3,3),  tcomplex ((treal)  1.468693939915887e+000, (treal) 2.287355287178844e+000), "scbmatrix::exp", os, __LINE__, dPessimisticSp);
        CheckComplex (scm(4,3),  tcomplex ((treal)  -8.186613472629570e-001, (treal) 3.756049227094730e+000), "scbmatrix::exp", os, __LINE__, dPessimisticSp);

        CheckComplex (scm(1,4),  tcomplex ((treal)  0., (treal) 0.), "scbmatrix::exp", os, __LINE__, dPessimisticSp);
        CheckComplex (scm(2,4),  tcomplex ((treal)  0., (treal) 0.), "scbmatrix::exp", os, __LINE__, dPessimisticSp);
        CheckComplex (scm(3,4),  tcomplex ((treal)  0., (treal) 0.), "scbmatrix::exp", os, __LINE__, dPessimisticSp);
        CheckComplex (scm(4,4),  tcomplex ((treal)  1.468693939915887e+000, (treal) 2.287355287178844e+000), "scbmatrix::exp", os, __LINE__, dPessimisticSp);

        r1 = (treal) 1.217;
        r2 = (treal) -.179;
        rv << rv2 << rv1;
        rv2.normalize();

        rv1.gemv (true, srm, r1, rv2, r2);
        rv = rv2 * srm * r1 + r2 * rv;
        CheckReal    ((rv - rv1).norm(), (treal) 0., "rvector::gemv", os, __LINE__, dPessimisticSp);

        rv = rv1;
        rv1.gemv (false, srm, r1, rv2, r2);
        rv = srm * rv2 * r1 + r2 * rv;
        CheckReal    ((rv - rv1).norm(), (treal) 0., "rvector::gemv", os, __LINE__, dPessimisticSp);

        rv = rv1;
        rv1.gbmv (true, srbm, r1, rv2, r2);
        rv = rv2 * srbm * r1 + r2 * rv;
        CheckReal    ((rv - rv1).norm(), (treal) 0., "rvector::gbmv", os, __LINE__, dPessimisticSp);

        rv = rv1;
        rv1.gbmv (false, srbm, r1, rv2, r2);
        rv = srbm * rv2 * r1 + r2 * rv;
        CheckReal    ((rv - rv1).norm(), (treal) 0., "rvector::gbmv", os, __LINE__, dPessimisticSp);

        cv2 << cv1;
        cv << cv2;
        cv2.normalize();

        cv1.gemv (true, scm, cr1, cv2, cr2);
        cv = cv2 * scm * cr1 + cr2 * cv;
        CheckReal    ((cv - cv1).norm(), (treal) 0., "cvector::gemv", os, __LINE__, dPessimisticSp);

        cv = cv1;
        cv1.gemv (false, scm, cr1, cv2, cr2);
        cv = scm * cv2 * cr1 + cr2 * cv;
        CheckReal    ((cv - cv1).norm(), (treal) 0., "cvector::gemv", os, __LINE__, dPessimisticSp);

        cv = cv1;
        cv1.gbmv (true, scbm, cr1, cv2, cr2);
        cv = cv2 * scbm * cr1 + cr2 * cv;
        CheckReal    ((cv - cv1).norm(), (treal) 0., "cvector::gemv", os, __LINE__, dPessimisticSp);

        cv = cv1;
        cv1.gbmv (false, scbm, cr1, cv2, cr2);

        cv = scbm * cv2 * cr1 + cr2 * cv;
        CheckReal    ((cv.normalize() - cv1.normalize()).norm(), (treal) 0., "cvector::gemv", os, __LINE__, dPessimisticSp);

        {
            treal a[] = {1., 2., 3., -4., 5., -6.};
            const rvector vr(a, 6);
            const cvector vc((std::complex<treal>*) a, 3);
            CheckReal    (vr.norm1(), (treal) 21., "rvector::norm1", os, __LINE__);
            CheckReal    (vc.norm1(), (treal) 15.04631765340644, "cvector::norm1", os, __LINE__, dPessimisticSp);
        }
        {
            treal a[] = {1., 2., 3., 4., 5., 6.};
            const rmatrix m1(a, 2, 3);
            rmatrix m2(2, 3);
            rmatrix m(2, 3);
            m2.set(1.);

            CheckReal    (m.sum(m1, m2)(2,2), (treal) 5., "rmatrix::sum", os, __LINE__);
            CheckReal    (m.sum(m, m2)(1,3), (treal) 7., "rmatrix::sum", os, __LINE__);
        }
        {
            treal a[] = {1., 2., 3., 4., 5., 6.};
            const rmatrix m1(a, 2, 3);
            rmatrix m2(2, 3);
            rmatrix m(2, 3);
            m2.set(1.);

            CheckReal    (m.diff(m1, m2)(2,2), (treal) 3., "rmatrix::sum", os, __LINE__);
            CheckReal    (m.diff(m, m2)(1,3), (treal) 3., "rmatrix::sum", os, __LINE__);
        }
        {
            treal a[] = {1., 2., 3., 4., 5., 6.,
                        7., 8., 9., 10., 11., 12.};
            const cmatrix ma ((std::complex<treal>*) a, 2, 3);
            cmatrix mb (2, 3);
            cmatrix m (2, 3);
            mb.set (std::complex<treal>(1.,1.));

            CheckComplex (m.sum(ma, mb)(2,2), tcomplex (8., 9.), "cmatrix::sum" , os, __LINE__);
            CheckComplex (m.sum(m, mb)(1,3), tcomplex (11., 12.), "cmatrix::sum" , os, __LINE__);
        }
        {
            treal a[] = {1., 2., 3., 4., 5., 6.,
                        7., 8., 9., 10., 11., 12.};
            const cmatrix ma ((std::complex<treal>*) a, 2, 3);
            cmatrix mb (2, 3);
            cmatrix m (2, 3);
            mb.set (std::complex<treal>(1.,1.));

            CheckComplex (m.diff(ma, mb)(2,2), tcomplex (6., 7.), "cmatrix::diff" , os, __LINE__);
            CheckComplex (m.diff(m, mb)(1,3), tcomplex (7., 8.), "cmatrix::diff" , os, __LINE__);
        }
        {
            treal a[] = {1., 2., 3., 4., 5., 6.};
            const srbmatrix m1(a,3,1,0);
            srbmatrix m2(3,1,0);
            srbmatrix m(3,1,0);
            m2.set(1.);
            CheckReal    (m.sum(m1, m2)(2,2), (treal) 4., "srbmatrix::sum", os, __LINE__);
            CheckReal    (m(2,3), (treal) 0., "srbmatrix::sum", os, __LINE__);
            CheckReal    (m(2,1), (treal) 3., "srbmatrix::sum", os, __LINE__);
            CheckReal    (m.sum(m, m2)(2,1), (treal) 4., "srbmatrix::sum", os, __LINE__);
            CheckReal    (m(3,3), (treal) 7., "srbmatrix::sum", os, __LINE__);
        }
        {
            treal a[] = {1., 2., 3., 4., 5., 6., 7., 8.,
                        9., 10., 11., 12.};
            const scbmatrix m1((std::complex<treal>*)a,3,1,0);
            scbmatrix m2(3,1,0);
            scbmatrix m(3,1,0);
            m2.set(std::complex<treal>(1.,1.));
            CheckComplex (m.sum(m1, m2)(2,1), tcomplex (4., 5.), "scbmatrix::sum" , os, __LINE__);
            CheckComplex (m(3,3), tcomplex (10., 11.), "scbmatrix::sum" , os, __LINE__);
            CheckComplex (m(2,3), tcomplex (0., 0.), "scbmatrix::sum" , os, __LINE__);

            CheckComplex (m.sum(m, m2)(2,1), tcomplex (5., 6.), "scbmatrix::sum" , os, __LINE__);
            CheckComplex (m(3,3), tcomplex (11., 12.), "scbmatrix::sum" , os, __LINE__);
            CheckComplex (m(2,3), tcomplex (0., 0.), "scbmatrix::sum" , os, __LINE__);
        }
        {
            treal a[] = {1., 2., 3., 4., 5., 6., 7., 8.,
                        9., 10., 11., 12.};
            const scbmatrix m1((std::complex<treal>*)a,3,1,0);
            scbmatrix m2(3,1,0);
            scbmatrix m(3,1,0);
            m2.set(std::complex<treal>(1.,1.));
            CheckComplex (m.diff(m1, m2)(2,1), tcomplex (2., 3.), "scbmatrix::diff" , os, __LINE__);
            CheckComplex (m(3,3), tcomplex (8., 9.), "scbmatrix::diff" , os, __LINE__);
            CheckComplex (m(2,3), tcomplex (0., 0.), "scbmatrix::diff" , os, __LINE__);
            CheckComplex (m.diff(m, m2)(2,1), tcomplex (1., 2.), "scbmatrix::diff" , os, __LINE__);
            CheckComplex (m(3,3), tcomplex (7., 8.), "scbmatrix::diff" , os, __LINE__);
            CheckComplex (m(2,3), tcomplex (0., 0.), "scbmatrix::diff" , os, __LINE__);
        }
        {
            treal a[] = {1., 2., 3., 2., 5., 6., 3., 6., 9.};
            const srsmatrix m1(a, 3);
            srsmatrix m2(3);
            srsmatrix m(3);
            m2.set(1.);

            CheckReal    (m.sum(m1, m2)(2,2), (treal) 6., "srsmatrix::sum", os, __LINE__);
            CheckReal    (m(1,3), (treal) 4., "srsmatrix::sum", os, __LINE__);
            CheckReal    (m(3,1), (treal) 4., "srsmatrix::sum", os, __LINE__);
            CheckReal    (m.sum(m, m2)(2,2), (treal) 7., "srsmatrix::sum", os, __LINE__);
            CheckReal    (m(1,3), (treal) 5., "srsmatrix::sum", os, __LINE__);
            CheckReal    (m(3,1), (treal) 5., "srsmatrix::sum", os, __LINE__);
        }
        {
            treal a[] = {1., 2., 3., 2., 5., 6., 3., 6., 9.};
            const srsmatrix m1(a, 3);
            srsmatrix m2(3);
            srsmatrix m(3);
            m2.set(1.);

            CheckReal    (m.diff(m1, m2)(2,2), (treal) 4., "srsmatrix::diff", os, __LINE__);
            CheckReal    (m(1,3), (treal) 2., "srsmatrix::diff", os, __LINE__);
            CheckReal    (m(3,1), (treal) 2., "srsmatrix::diff", os, __LINE__);
            CheckReal    (m.diff(m, m2)(2,2), (treal) 3., "srsmatrix::diff", os, __LINE__);
            CheckReal    (m(1,3), (treal) 1., "srsmatrix::diff", os, __LINE__);
            CheckReal    (m(3,1), (treal) 1., "srsmatrix::diff", os, __LINE__);
        }
        {
            treal a[] = {1., 0., 2., 1., -1., 2., 2., -1., 2., 0.,
                        0., 3., -1., -2., 0., -3., 3., 0.};
            treal b[] = {1., 0., 1., 1., 1., 1., 1., -1., 1., 0.,
                        1., 1., 1., -1., 1., -1., 1., 0.};
            schmatrix m1((std::complex<treal>*)a,3);
            schmatrix m2((std::complex<treal>*)b,3);
            schmatrix m(3);
            CheckComplex (m.sum(m1, m2)(2,1), tcomplex (3., 2.), "schmatrix::sum" , os, __LINE__);
            CheckComplex (m(3,3), tcomplex (4., 0.), "schmatrix::sum" , os, __LINE__);
            CheckComplex (m(2,3), tcomplex (1., -4.), "schmatrix::sum" , os, __LINE__);
            CheckComplex (m(3,2), tcomplex (1., 4.), "schmatrix::sum" , os, __LINE__);
        }
        {
            treal a[] = {1., 0., 2., 1., -1., 2., 2., -1., 2., 0.,
                        0., 3., -1., -2., 0., -3., 3., 0.};
            treal b[] = {1., 0., 1., 1., 1., 1., 1., -1., 1., 0.,
                        1., 1., 1., -1., 1., -1., 1., 0.};
            schmatrix m1((std::complex<treal>*)a,3);
            schmatrix m2((std::complex<treal>*)b,3);
            schmatrix m(3);
            CheckComplex (m.diff(m1, m2)(2,1), tcomplex (1., 0.), "schmatrix::diff" , os, __LINE__);
            CheckComplex (m(3,3), tcomplex (2., 0.), "schmatrix::diff" , os, __LINE__);
            CheckComplex (m(2,3), tcomplex (-1., -2.), "schmatrix::diff" , os, __LINE__);
            CheckComplex (m(3,2), tcomplex (-1., 2.), "schmatrix::diff" , os, __LINE__);
        }

        {
            treal a[] = {3., 0., 2., 1., -1., 2., 2., -1., 3., 0.,
                        0., 3., -1., -2., 0., -3., 5., 0.};
            const schmatrix m((std::complex<treal>*)a,3);
            scmatrix h = m.cholesky();
            CheckReal ((~h * h - m).norm(), (treal) 0., "schmatrix::cholesky" , os, __LINE__, dPessimisticSp);
        }
        {
            treal a[] = {1., 0., 2., 1., -1., 2., 2., -1., 2., 0.,
                        0., 3., -1., -2., 0., -3., 3., 0.};
            schmatrix m((std::complex<treal>*)a,3);
            scmatrix me(3);
            rvector v(3);

            v = m.eig(me);
            cvector vc(v);

            CheckReal ((m * me(1) - me(1) * vc(1)).norm(), (treal) 0., "schmatrix::eig" , os, __LINE__, dPessimisticSp);
            CheckReal ((m * me(2) - me(2) * vc(2)).norm(), (treal) 0., "schmatrix::eig" , os, __LINE__, dPessimisticSp);
            CheckReal ((m * me(3) - me(3) * vc(3)).norm(), (treal) 0., "schmatrix::eig" , os, __LINE__, dPessimisticSp);

            CheckComplex (me(1) % me(2), tcomplex (0., 0.), "schmatrix::eig" , os, __LINE__, dPessimisticSp);
            CheckComplex (me(2) % me(3), tcomplex (0., 0.), "schmatrix::eig" , os, __LINE__, dPessimisticSp);
            CheckComplex (me(1) % me(3), tcomplex (0., 0.), "schmatrix::eig" , os, __LINE__, dPessimisticSp);
        }
        {
            std::cout.setf (std::ios::scientific | std::ios::left |
                            std::ios::showpos);
            std::cout.precision (10);
            treal a[] = {1., 0., 2., 1., -1., 2., 2., -1., 2., 0.,
                        0., 3., -1., -2., 0., -3., 3., 0.};
            schmatrix m((std::complex<treal>*)a,3);
            treal re[]={2.2,1.3,1.1,-0.9,0.2,-0.45,45.,-30.,10.,3.,1.13};
            const rvector vr(re, 11);
            schmatrix mp(3);

            mp.polynom (m, vr);

            CheckComplex (mp(1,1), tcomplex (1.231954875800000e+008, 0.), "schmatrix::polynom" , os, __LINE__, dPessimisticSp);
            CheckComplex (mp(2,1), tcomplex (1.417932391600000e+008, 7.089661958000000e+007), "schmatrix::polynom" , os, __LINE__, dPessimisticSp);
            CheckComplex (mp(3,1), tcomplex (-8.080273845999999e+007, 1.616054769200000e+008), "schmatrix::polynom" , os, __LINE__, dPessimisticSp);
            CheckComplex (mp(1,2), tcomplex (1.417932391600000e+008, -7.089661958000000e+007), "schmatrix::polynom" , os, __LINE__, dPessimisticSp);
            CheckComplex (mp(2,2), tcomplex (2.039982260400000e+008, 0.), "schmatrix::polynom" , os, __LINE__, dPessimisticSp);
            CheckComplex (mp(3,2), tcomplex (0., 2.325020965000000e+008), "schmatrix::polynom" , os, __LINE__, dPessimisticSp);
            CheckComplex (mp(1,3), tcomplex (-8.080273845999999e+007, -1.616054769200000e+008), "schmatrix::polynom" , os, __LINE__, dPessimisticSp);
            CheckComplex (mp(2,3), tcomplex (0., -2.325020965000000e+008), "schmatrix::polynom" , os, __LINE__, dPessimisticSp);
            CheckComplex (mp(3,3), tcomplex (2.649887267400000e+008, 0.), "schmatrix::polynom" , os, __LINE__, dPessimisticSp);
/*
 Column 1 

     1.231954875800000e+008                         
     1.417932391600000e+008 +7.089661958000000e+007i
    -8.080273845999999e+007 +1.616054769200000e+008i

  Column 2 

     1.417932391600000e+008 -7.089661958000000e+007i
     2.039982260400000e+008                         
                          0 +2.325020965000000e+008i

  Column 3 

    -8.080273845999999e+007 -1.616054769200000e+008i
                          0 -2.325020965000000e+008i
     2.649887267400000e+008                         
*/
        }

        {
            treal a[] = {1., 0., 2., 1., -1., 2., 2., -1., 2., 0.,
                        0., 3., -1., -2., 0., -3., 3., 0.};
            schmatrix m((std::complex<treal>*)a,3);
            treal re[]={2.2,1.3,1.1,-0.9,0.2,-0.45,45.,-30.,10.,3.,1.13};
            treal im[]={0.5,-2,0,1,3,-3.,30.,0.,-9.,0.,1.};

            const cvector vc(re, im, 11);
            scmatrix mp(3);

            mp.polynom (m, vc);

            CheckComplex (mp(1,1), tcomplex (1.231954875800000e+008, 6.128500650000000e+007), "schmatrix::polynom" , os, __LINE__, dPessimisticSp);
            CheckComplex (mp(2,1), tcomplex (1.065249031600000e+008, 1.414332915800000e+008), "schmatrix::polynom" , os, __LINE__, dPessimisticSp);
            CheckComplex (mp(3,1), tcomplex (-1.611952344600000e+008, 1.214092289200000e+008), "schmatrix::polynom" , os, __LINE__, dPessimisticSp);
            CheckComplex (mp(1,2), tcomplex (1.770615751600000e+008, -3.599475799999982e+005), "schmatrix::polynom" , os, __LINE__, dPessimisticSp);
            CheckComplex (mp(2,2), tcomplex (2.039982260400000e+008, 1.014812545000000e+008), "schmatrix::polynom" , os, __LINE__, dPessimisticSp);
            CheckComplex (mp(3,2), tcomplex (-1.156608320000000e+008, 2.325020965000000e+008), "schmatrix::polynom" , os, __LINE__, dPessimisticSp);
            CheckComplex (mp(1,3), tcomplex (-4.102424600000009e+005, -2.018017249200000e+008), "schmatrix::polynom" , os, __LINE__, dPessimisticSp);
            CheckComplex (mp(2,3), tcomplex (1.156608320000000e+008, -2.325020965000000e+008), "schmatrix::polynom" , os, __LINE__, dPessimisticSp);
            CheckComplex (mp(3,3), tcomplex (2.649887267400000e+008, 1.318216785000000e+008), "schmatrix::polynom" , os, __LINE__, dPessimisticSp);
/*
  Column 1 

     1.231954875800000e+008 +6.128500650000000e+007i
     1.065249031600000e+008 +1.414332915800000e+008i
    -1.611952344600000e+008 +1.214092289200000e+008i

  Column 2 

     1.770615751600000e+008 -3.599475799999982e+005i
     2.039982260400000e+008 +1.014812545000000e+008i
    -1.156608320000000e+008 +2.325020965000000e+008i

  Column 3 

    -4.102424600000009e+005 -2.018017249200000e+008i
     1.156608320000000e+008 -2.325020965000000e+008i
     2.649887267400000e+008 +1.318216785000000e+008i
*/

        }
        {
            treal a[] = {1., 0., 2., 1., -1., 2., 2., -1., 2., 0.,
                        0., 3., -1., -2., 0., -3., 3., 0.};
            schmatrix m((std::complex<treal>*)a,3);
            schmatrix me(3);
            me.exp(m);

            CheckComplex (me(1,1), tcomplex (2.673228708371998e+002, 0.), "schmatrix::exp" , os, __LINE__, dPessimisticSp);
            CheckComplex (me(2,1), tcomplex (3.071187567026802e+002, 1.535593783513401e+002), "schmatrix::exp" , os, __LINE__, dPessimisticSp);
            CheckComplex (me(3,1), tcomplex (-1.749365628720764e+002, 3.498731257441527e+002), "schmatrix::exp" , os, __LINE__, dPessimisticSp);
            CheckComplex (me(1,2), tcomplex (3.071187567026802e+002, -1.535593783513401e+002), "schmatrix::exp" , os, __LINE__, dPessimisticSp);
            CheckComplex (me(2,2), tcomplex (4.422594337092769e+002, 0.), "schmatrix::exp" , os, __LINE__, dPessimisticSp);
            CheckComplex (me(3,2), tcomplex (3.549798266275454e-015, 5.034325040954932e+002), "schmatrix::exp" , os, __LINE__, dPessimisticSp);
            CheckComplex (me(1,3), tcomplex (-1.749365628720763e+002, -3.498731257441526e+002), "schmatrix::exp" , os, __LINE__, dPessimisticSp);
            CheckComplex (me(2,3), tcomplex (-1.776065298147746e-014, -5.034325040954931e+002), "schmatrix::exp" , os, __LINE__, dPessimisticSp);
            CheckComplex (me(3,3), tcomplex (5.744416275398801e+002, 0.), "schmatrix::exp" , os, __LINE__, dPessimisticSp);

/*
            Column 1 

                2.673228708371998e+002 -7.105427357601002e-015i
                3.071187567026802e+002 +1.535593783513401e+002i
                -1.749365628720764e+002 +3.498731257441527e+002i

            Column 2 

                3.071187567026802e+002 -1.535593783513401e+002i
                4.422594337092769e+002 -5.489286670342458e-016i
                3.549798266275454e-015 +5.034325040954932e+002i

            Column 3 

                -1.749365628720763e+002 -3.498731257441526e+002i
                -1.776065298147746e-014 -5.034325040954931e+002i
                5.744416275398801e+002 -2.096383162906490e-014i
*/
        }
        {
            treal a[] = {1., 0., 2., 1., -1., 2., 2., -1., 2., 0.,
                        0., 3., -1., -2., 0., -3., 3., 0.};
            schmatrix m((std::complex<treal>*)a,3);
            const schmatrix mi = m.inv();
            CheckReal ((mi * m - eye_complex(3)).norm(), (treal) 0., "schmatrix::inv" , os, __LINE__, dPessimisticSp);
        }
        {
            treal a1[] = {1., -1., 2., 2., 3., -3.};
            treal a2[] = {-2., 1., 1., -2., 1., -3.};
            cvector v1((std::complex<treal>*)a1, 3);
            cvector v2((std::complex<treal>*)a2, 3);
            schmatrix mh(3);
            mh.set_real(1.); mh.set(1,3,std::complex<treal>(4.,1.));
            mh.her2k (std::complex<treal>(2.,-1.), v1, v2, 
                      std::complex<treal>(-1.,-1.));
            CheckComplex (mh(1,1), tcomplex ((treal) -11., (treal) 0.), "schmatrix::her2k" , os, __LINE__);
            CheckComplex (mh(2,3), tcomplex ((treal) 20., (treal) 23.), "schmatrix::her2k" , os, __LINE__);

            cmatrix m1(3,2);
            cmatrix m2(3,2);
            m1.set_real(2.); m1.set_imag(-1.);
            m2.set_real(-3.); m2.set_imag(1.);
            mh.her2k (false, std::complex<treal>(2.,1.), m1, m2, 
                             std::complex<treal>(3.,-2.));
            CheckComplex (mh(1,1), tcomplex ((treal) -93., (treal) 0.), "schmatrix::her2k" , os, __LINE__);
            CheckComplex (mh(2,3), tcomplex ((treal) 0., (treal) 69.), "schmatrix::her2k" , os, __LINE__);

            schmatrix mh2(2);
            mh2.her2k (true, std::complex<treal>(1.,1.), m1, m2, 
                            std::complex<treal>(2.,-3.));
            CheckComplex (mh2(1,1), tcomplex ((treal) -36., (treal) 0.), "schmatrix::her2k" , os, __LINE__);
            CheckComplex (mh2(2,1), tcomplex ((treal) -36., (treal) 0.), "schmatrix::her2k" , os, __LINE__);
        }
        {
            treal a[] = {1., -1., 2., 2., 3., -3.};
            cvector v((std::complex<treal>*)a, 3);
            schmatrix mh(3);
            mh.set_real(1.);
            mh.herk (std::complex<treal>(2.,-2.), v, 
                     std::complex<treal>(1.,1.));
            CheckComplex (mh(1,1), tcomplex ((treal) 5., (treal) 0.), "schmatrix::herk" , os, __LINE__);
            CheckComplex (mh(1,2), tcomplex ((treal) 1., (treal) -8.), "schmatrix::herk" , os, __LINE__);

            cmatrix m(3,2);
            m(1) = v;
            m(2).set(std::complex<treal>(-1.,1.));
            mh.herk (false, std::complex<treal>(2.,-1.), m, 0.);
            CheckComplex (mh(1,1), tcomplex ((treal) 8., (treal) 0.), "schmatrix::herk" , os, __LINE__);
            CheckComplex (mh(1,2), tcomplex ((treal) 4., (treal) -8.), "schmatrix::herk" , os, __LINE__);

            schmatrix mh2(2);
            mh2.herk (true, std::complex<treal>(1.,1.), m, 
                    std::complex<treal>(0.,1.));
            CheckComplex (mh2(1,1), tcomplex ((treal) 28., (treal) 0.), "schmatrix::herk" , os, __LINE__);
            CheckComplex (mh2(1,2), tcomplex ((treal) -8., (treal) 4.), "schmatrix::herk" , os, __LINE__);
            CheckComplex (mh2(2,2), tcomplex ((treal) 6., (treal) 0.), "schmatrix::herk" , os, __LINE__);
        }
        {
            treal a[] = {1., 0., 2., 1., -1., 2., 2., -1., 2., 0.,
                        0., 3., -1., -2., 0., -3., 3., 0.};
            schmatrix m((std::complex<treal>*)a,3);
            rvector v(3);
            v.set(7.7);
            m.set_main_diag(v);
            CheckComplex (m(1,1), tcomplex ((treal) 7.7, (treal) 0.), "schmatrix::set_main_diag" , os, __LINE__);
            CheckComplex (m(3,3), tcomplex ((treal) 7.7, (treal) 0.), "schmatrix::set_main_diag" , os, __LINE__);
            CheckComplex (m(1,2), tcomplex ((treal) 2., (treal) -1.), "schmatrix::set_main_diag" , os, __LINE__);
            CheckComplex (m(2,1), tcomplex ((treal) 2., (treal) 1.), "schmatrix::set_main_diag" , os, __LINE__);
        }
        {
            treal a[] = {1., 0., 2., 1., -1., 2., 2., -1., 2., 0.,
                        0., 3., -1., -2., 0., -3., 3., 0.};
            schmatrix m((std::complex<treal>*)a,3);
            cvector v(2);
            v.set(std::complex<treal>(7.,7.));
            m.set_diag(1, v);
            CheckComplex (m(1,2), tcomplex ((treal) 7., (treal) 7.), "schmatrix::set_diag" , os, __LINE__);
            CheckComplex (m(2,1), tcomplex ((treal) 7., (treal) -7.), "schmatrix::set_diag" , os, __LINE__);
            CheckComplex (m(2,2), tcomplex ((treal) 2., (treal) 0.), "schmatrix::set_diag" , os, __LINE__);
        }
        {
            treal a[] = {1., 2., 1., 2., 5., -1., 1., -1., 20.};
            const srsmatrix m(a, 3);
            srmatrix h = m.cholesky();
            CheckReal ((~h * h - m).norm(), (treal) 0., "srsmatrix::cholesky" , os, __LINE__, dPessimisticSp);
        }
        {
            treal a[] = {1., 2., 1., 2., 0., -1., 1., -1., 2.};
            const srsmatrix m(a, 3);
            srmatrix me(3);
            rvector v(3);

            v = m.eig(me);

            CheckReal ((m * me(1) - me(1) * v(1)).norm(), (treal) 0., "srsmatrix::eig" , os, __LINE__, dPessimisticSp);
            CheckReal ((m * me(2) - me(2) * v(2)).norm(), (treal) 0., "srsmatrix::eig" , os, __LINE__, dPessimisticSp);
            CheckReal ((m * me(3) - me(3) * v(3)).norm(), (treal) 0., "srsmatrix::eig" , os, __LINE__, dPessimisticSp);
        }
        {
            std::cout.setf (std::ios::scientific | std::ios::left); 
            std::cout.precision (7);
            treal a[] = {1., 2., 1., 2., 0., -1., 1., -1., 2.};
            treal av[] = {2.2, 1.3, 1.1, -0.9, 0.2,
                        -0.45, 45, -30, 10, 3, 3.2};
            const rvector v(av, 11);
            const srsmatrix m(a, 3);

            const srsmatrix mp = m.polynom (v);
            CheckReal (mp(1,1), (treal) 6.212740000000001e+004, "srsmatrix::polynom" , os, __LINE__, dPessimisticSp);
            CheckReal (mp(2,1), (treal) 2.399800000000000e+004, "srsmatrix::polynom" , os, __LINE__, dPessimisticSp);
            CheckReal (mp(3,1), (treal) 3.410055000000000e+004, "srsmatrix::polynom" , os, __LINE__, dPessimisticSp);
            CheckReal (mp(1,2), (treal) 2.399800000000000e+004, "srsmatrix::polynom" , os, __LINE__, dPessimisticSp);
            CheckReal (mp(2,2), (treal) 2.802685000000000e+004, "srsmatrix::polynom" , os, __LINE__, dPessimisticSp);
            CheckReal (mp(3,2), (treal) 1.010255000000000e+004, "srsmatrix::polynom" , os, __LINE__, dPessimisticSp);
            CheckReal (mp(1,3), (treal) 3.410055000000000e+004, "srsmatrix::polynom" , os, __LINE__, dPessimisticSp);
            CheckReal (mp(2,3), (treal) 1.010255000000000e+004, "srsmatrix::polynom" , os, __LINE__, dPessimisticSp);
            CheckReal (mp(3,3), (treal) 5.202485000000000e+004, "srsmatrix::polynom" , os, __LINE__, dPessimisticSp);

/*
            Columns 1 through 2 

                6.212740000000001e+004    2.399800000000000e+004
                2.399800000000000e+004    2.802685000000000e+004
                3.410055000000000e+004    1.010255000000000e+004

            Column 3 

                3.410055000000000e+004
                1.010255000000000e+004
                5.202485000000000e+004
*/
        }
        {
            treal a[] = {1., 2., 1., 2., 0., -1., 1., -1., 2.};
            const srsmatrix m(a, 3);
            const srsmatrix me = m.exp();

            CheckReal (me(1,1), (treal) 9.198262499129212e+000, "srsmatrix::exp" , os, __LINE__, dPessimisticSp);
            CheckReal (me(2,1), (treal) 5.558586002658865e+000, "srsmatrix::exp" , os, __LINE__, dPessimisticSp);
            CheckReal (me(3,1), (treal) 3.852443363622600e+000, "srsmatrix::exp" , os, __LINE__, dPessimisticSp);
            CheckReal (me(1,2), (treal) 5.558586002658862e+000, "srsmatrix::exp" , os, __LINE__, dPessimisticSp);
            CheckReal (me(2,2), (treal) 5.345819135506588e+000, "srsmatrix::exp" , os, __LINE__, dPessimisticSp);
            CheckReal (me(3,2), (treal) -1.706142639036258e+000, "srsmatrix::exp" , os, __LINE__, dPessimisticSp);
            CheckReal (me(1,3), (treal) 3.852443363622601e+000, "srsmatrix::exp" , os, __LINE__, dPessimisticSp);
            CheckReal (me(2,3), (treal) -1.706142639036260e+000, "srsmatrix::exp" , os, __LINE__, dPessimisticSp);
            CheckReal (me(3,3), (treal) 1.090440513816545e+001, "srsmatrix::exp" , os, __LINE__, dPessimisticSp);
/*
            Columns 1 through 2 

                9.198262499129212e+000    5.558586002658862e+000
                5.558586002658865e+000    5.345819135506588e+000
                3.852443363622600e+000   -1.706142639036258e+000

            Column 3 

                3.852443363622601e+000
               -1.706142639036260e+000
                1.090440513816545e+001
*/
        }
        {
            treal a1[] = {1., 2., 3., 4.};
            treal a2[] = {1., 2., 3., 4.};
            rvector v1(a1,4);
            rvector v2(a2,4);
            srsmatrix ms(4);
            ms.set(1.);
            ms.syr2k (2., v1, v2, 1.);
            CheckReal (ms(4,4), (treal) 65., "srsmatrix::syr2k" , os, __LINE__, dPessimisticSp);
            CheckReal (ms(1,4), (treal) 17., "srsmatrix::syr2k" , os, __LINE__, dPessimisticSp);

            rmatrix m1(4,2);
            rmatrix m2(4,2);
            m1.set(1.);
            m2.set(2.);
            ms.syr2k (false, 2., m1, m2, 0.);
            CheckReal (ms(4,4), (treal) 16., "srsmatrix::syr2k" , os, __LINE__, dPessimisticSp);
            CheckReal (ms(1,4), (treal) 16., "srsmatrix::syr2k" , os, __LINE__, dPessimisticSp);

            srsmatrix ms2(2);
            ms2.syr2k (true, 1., m1, m2, 0.);
            CheckReal (ms2(2,2), (treal) 16., "srsmatrix::syr2k" , os, __LINE__, dPessimisticSp);
            CheckReal (ms2(1,2), (treal) 16., "srsmatrix::syr2k" , os, __LINE__, dPessimisticSp);
        }
        {
            treal a[] = {1., 2., 3., 4.};
            rvector v(a,4);
            srsmatrix ms(4);
            ms.set(1.);
            ms.syrk (2., v, 1.);
            CheckReal (ms(4,4), (treal) 33., "srsmatrix::syrk" , os, __LINE__, dPessimisticSp);
            CheckReal (ms(1,4), (treal) 9., "srsmatrix::syrk" , os, __LINE__, dPessimisticSp);

            rmatrix m(4,2);
            m(1) = v;
            m(2).set(1.);
            ms.syrk (false, 2., m, 0.);
            CheckReal (ms(4,4), (treal) 34., "srsmatrix::syrk" , os, __LINE__, dPessimisticSp);
            CheckReal (ms(1,4), (treal) 10., "srsmatrix::syrk" , os, __LINE__, dPessimisticSp);

            srsmatrix ms2(2);
            ms2.syrk (true, 1., m, 0.);
            CheckReal (ms2(1,1), (treal) 30., "srsmatrix::syrk" , os, __LINE__, dPessimisticSp);
            CheckReal (ms2(1,2), (treal) 10., "srsmatrix::syrk" , os, __LINE__, dPessimisticSp);
            CheckReal (ms2(2,2), (treal) 4., "srsmatrix::syrk" , os, __LINE__, dPessimisticSp);
        }
        {
            treal a[] = {1., 2., 3., 4., 5., 6., 7., 8., 9.,
                        10., 11., 12.};
            scbmatrix ma((std::complex<treal>*)a,3,1,0);
            scbmatrix mLU(3,1,0);
            cmatrix  mb1(3,2); cvector vb1(3);
            cmatrix  mb2(3,2); cvector vb2(3);
            cmatrix  mx1(3,2); cvector vx1(3);
            cmatrix  mx2(3,2); cvector vx2(3);
            iarray   nPivots(3);
            treal   dErr = 0.;
            mb1.randomize_real(-1.,3.); mb1.randomize_imag(1.,5.);
            mb2.randomize_real(-2.,5.); mb2.randomize_imag(-3.,0.);
            vb1.randomize_real(-2.,4.); vb1.randomize_imag(-4.,1.);
            vb2.randomize_real(-3.,1.); vb2.randomize_imag(4.,5.);

            mLU.low_up(ma, nPivots);
            mx1 = ma.solve_lu (mLU, nPivots, mb1, dErr);
            CheckReal (dErr, (treal) 0., "scbmatrix::solve_lu" , os, __LINE__, dPessimisticSp);
            mx2 = ma.solve_lu (mLU, nPivots, mb2);
            CheckReal ((ma * mx1 - mb1).norm(), (treal) 0., "scbmatrix::solve_lu" , os, __LINE__, dPessimisticSp);
            CheckReal ((ma * mx2 - mb2).norm(), (treal) 0., "scbmatrix::solve_lu" , os, __LINE__, dPessimisticSp);

            vx1 = ma.solve_lu (mLU, nPivots, vb1, dErr);
            vx2 = ma.solve_lu (mLU, nPivots, vb2);
            CheckReal ((ma * vx1 - vb1).norm(), (treal) 0., "scbmatrix::solve_lu" , os, __LINE__, dPessimisticSp);
            CheckReal ((ma * vx2 - vb2).norm(), (treal) 0., "scbmatrix::solve_lu" , os, __LINE__, dPessimisticSp);
        }
        {
            treal a[] = {1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12.};
            scbmatrix m((std::complex<treal>*)a,3,1,0);
            m.resize_lu (0,1);
            m.diag(1).set(std::complex<treal>(9.,9.));
            CheckComplex (m(1,2), tcomplex ((treal) 9., (treal) 9.), "scbmatrix::resize_lu" , os, __LINE__);
            CheckComplex (m(2,1), tcomplex ((treal) 0., (treal) 0.), "scbmatrix::resize_lu" , os, __LINE__);
            CheckComplex (m(3,3), tcomplex ((treal) 9., (treal) 10.), "scbmatrix::resize_lu" , os, __LINE__);
        }
        {
            treal a[] = {1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12.};
            scbmatrix m ((std::complex<treal>*)a,3,1,0);
            CheckReal (m.real()(2,2), (treal) 5., "scbmatrix::real" , os, __LINE__, dPessimisticSp);
            CheckReal (m.imag()(2,2), (treal) 6., "scbmatrix::imag" , os, __LINE__, dPessimisticSp);
        }
        {
            treal a[] = {1., 2., 3., 4., 5., 6., 7., 8.};
            const srbmatrix m(a,4,1,0);
            scbmatrix mr(m), mi(m, false);
            CheckComplex (mr(4,3), tcomplex ((treal) 6., (treal) 0.), "scbmatrix(srbmatrix,bool)" , os, __LINE__);
            CheckComplex (mr(1,4), tcomplex ((treal) 0., (treal) 0.), "scbmatrix(srbmatrix,bool)" , os, __LINE__);
            CheckComplex (mi(4,3), tcomplex ((treal) 0., (treal) 6.), "scbmatrix(srbmatrix,bool)" , os, __LINE__);
            CheckComplex (mi(1,4), tcomplex ((treal) 0., (treal) 0.), "scbmatrix(srbmatrix,bool)" , os, __LINE__);
        }
        {
            srbmatrix mr(4,1,0), mi(4,1,0);
            mr.set(1.);
            mi.set(2.);
            const scbmatrix m(mr,mi);
            CheckComplex (m(2,1), tcomplex ((treal) 1., (treal) 2.), "scbmatrix(srbmatrix,srbmatrix)" , os, __LINE__);
            CheckComplex (m(1,2), tcomplex ((treal) 0., (treal) 0.), "scbmatrix(srbmatrix,srbmatrix)" , os, __LINE__);
        }
        {
            treal a[] = {1., 2., 3., 4., 5., 6., 7., 8.};
            srbmatrix ma(a,4,1,0);
            srbmatrix mLU(4,1,0);
            rmatrix  mb1(4,2); rvector vb1(4);
            rmatrix  mb2(4,2); rvector vb2(4);
            rmatrix  mx1(4,2); rvector vx1(4);
            rmatrix  mx2(4,2); rvector vx2(4);
            iarray   nPivots(4);
            treal   dErr = 0.;
            mb1.randomize(-1.,3.); vb1.randomize(-2.,4.);
            mb2.randomize(-2.,5.); vb2.randomize(-3.,1.);

            mLU.low_up(ma, nPivots);
            mx1 = ma.solve_lu (mLU, nPivots, mb1, dErr);
            CheckReal (dErr, (treal) 0., "srbmatrix::solve_lu" , os, __LINE__, dPessimisticSp);
            mx2 = ma.solve_lu (mLU, nPivots, mb2);
            CheckReal ((ma * mx1 - mb1).norm(), (treal) 0., "srbmatrix::solve_lu" , os, __LINE__, dPessimisticSp);
            CheckReal ((ma * mx2 - mb2).norm(), (treal) 0., "srbmatrix::solve_lu" , os, __LINE__, dPessimisticSp);

            vx1 = ma.solve_lu (mLU, nPivots, vb1, dErr);
            CheckReal (dErr, (treal) 0., "srbmatrix::solve_lu" , os, __LINE__, dPessimisticSp);
            vx2 = ma.solve_lu (mLU, nPivots, vb2);
            CheckReal ((ma * vx1 - vb1).norm(), (treal) 0., "srbmatrix::solve_lu" , os, __LINE__, dPessimisticSp);
            CheckReal ((ma * vx2 - vb2).norm(), (treal) 0., "srbmatrix::solve_lu" , os, __LINE__, dPessimisticSp);
        }
        {
            treal a[] = {1., -1., 1., 2., -2., 1., 3., -2., 1.};
            srmatrix ma(a,3);
            srmatrix mLU(3);
            rmatrix  mb1(3,2); rvector vb1(3);
            rmatrix  mb2(3,2); rvector vb2(3);
            rmatrix  mx1(3,2); rvector vx1(3);
            rmatrix  mx2(3,2); rvector vx2(3);
            iarray   nPivots(3);
            treal   dErr = 0.;
            mb1.randomize(-1.,3.); vb1.randomize(-2.,4.);
            mb2.randomize(-2.,5.); vb2.randomize(-3.,1.);

            mLU.low_up(ma, nPivots);
            mx1 = ma.solve_lu (mLU, nPivots, mb1, dErr);
            CheckReal (dErr, (treal) 0., "rmatrix::solve_lu" , os, __LINE__, dPessimisticSp);

            mx2 = ma.solve_lu (mLU, nPivots, mb2);
            CheckReal ((ma * mx1 - mb1).norm(), (treal) 0., "rmatrix::solve_lu" , os, __LINE__, dPessimisticSp);
            CheckReal ((ma * mx2 - mb2).norm(), (treal) 0., "rmatrix::solve_lu" , os, __LINE__, dPessimisticSp);

            vx1 = ma.solve_lu (mLU, nPivots, vb1, dErr);
            CheckReal (dErr, (treal) 0., "rmatrix::solve_lu" , os, __LINE__, dPessimisticSp);
            vx2 = ma.solve_lu (mLU, nPivots, vb2);
            CheckReal ((ma * vx1 - vb1).norm(), (treal) 0., "rmatrix::solve_lu" , os, __LINE__, dPessimisticSp);
            CheckReal ((ma * vx2 - vb2).norm(), (treal) 0., "rmatrix::solve_lu" , os, __LINE__, dPessimisticSp);
        }
        {
            std::complex<treal> alpha = std::complex<treal>(1.3,0.21);
            std::complex<treal> beta = std::complex<treal>(0.5,-0.1);
            cmatrix m1(2,3);
            cmatrix m2(3,2);
            schmatrix ms(2);
            cmatrix m(2,3);
            m.randomize_real(-1., 2.); m.randomize_imag(1., 3.); 
            m1.randomize_real(-1., 3.); m1.randomize_imag(1., 2.);
            m2.randomize_real(0., 2.); m2.randomize_imag(-3., -1.);
            ms.randomize_real(-3., 1.); ms.randomize_imag(-1.3, 4.);

            cmatrix mr = ms * m1 * alpha + m * beta;
            CheckReal ((mr - m.hemm (true, ms, m1, alpha, beta)).norm(), (treal) 0., "cmatrix::hemm" , os, __LINE__, dPessimisticSp);

            m.resize(3,2);
            m.randomize_real(-1.4, 1.3); m.randomize_imag(1.1, 3.); 
            cmatrix mr2 = m2 * ms * alpha + m * beta;
            CheckReal ((mr2 - m.hemm (false, ms, m2, alpha, beta)).norm(), (treal) 0., "cmatrix::hemm" , os, __LINE__, dPessimisticSp);
        }
        {
            std::complex<treal> alpha = std::complex<treal>(1.1,2.1);
            std::complex<treal> beta = std::complex<treal>(0.71,0.12);
            cmatrix m1(4,3); cmatrix m2(4,3);
            cmatrix m(3,3);
            m.randomize_real(-1., 2.); m1.randomize_real(-1., 3.); m2.randomize_real(0., 2.);
            m.randomize_imag(1., 3.); m1.randomize_imag(-2., 4.); m2.randomize_imag(-3., 2.);
            cmatrix mr = ~m1 * m2 * alpha + m * beta;
            CheckReal ((mr - m.gemm(m1, true, m2, false, alpha, beta)).norm(), (treal) 0., "cmatrix::gemm" , os, __LINE__, dPessimisticSp);
        }
        {
            std::complex<treal> alpha = std::complex<treal>(1.2,4.11);
            cmatrix m(3,2);
            cvector vc(3);
            cvector vr(2);
            m.randomize_real(-1., 2.); vc.randomize_real(-1., 3.); vr.randomize_real(0., 2.);
            m.randomize_imag(-3., 2.); vc.randomize_imag(1., 3.); vr.randomize_imag(-1., 2.);
            cmatrix mr = m + vc.rank1update_u (vr) * alpha;
            CheckReal ((mr - m.geru(alpha, vc, vr)).norm(), (treal) 0., "cmatrix::geru" , os, __LINE__, dPessimisticSp);
        }
        {
            std::complex<treal> alpha = std::complex<treal>(1.2,4.11);
            cmatrix m(3,2);
            cvector vc(3);
            cvector vr(2);
            m.randomize_real(-1., 2.); vc.randomize_real(-1., 3.); vr.randomize_real(0., 2.);
            m.randomize_imag(-3., 2.); vc.randomize_imag(1., 3.); vr.randomize_imag(-1., 2.);

            cmatrix mr = m + vc.rank1update_c (vr) * alpha;
            CheckReal ((mr - m.gerc(alpha, vc, vr)).norm(), (treal) 0., "cmatrix::gerc" , os, __LINE__, dPessimisticSp);
        }
        {
            treal a[] = {1., 2., 3., 4., 5., 6.,
                        7., 8., 9., 10., 11., 12.,
                        13., 14., 15., 16., 17., 18.};
            cmatrix m (2, 3);
            const scmatrix ms((std::complex<treal>*)a, 3);

            m.diag(-1).set(std::complex<treal>(1.,1.));
            m.diag(0).set(std::complex<treal>(2.,2.));
            m.diag(1).set(std::complex<treal>(3.,3.));
            m.diag(2).set(std::complex<treal>(4.,4.));
            CheckComplex (m(1,2), tcomplex ((treal) 3., (treal) 3.), "cmatrix::diag" , os, __LINE__);
            CheckComplex (m(2,3), tcomplex ((treal) 3., (treal) 3.), "cmatrix::diag" , os, __LINE__);
            CheckComplex (ms.diag(0)[2], tcomplex ((treal) 9., (treal) 10.), "cmatrix::diag" , os, __LINE__);
        }
        {
            treal alpha = 1.3;
            treal beta = -0.7;
            rmatrix m1(3,4);
            rmatrix m2(4,3);
            srsmatrix ms(3);
            rmatrix m(3,4);
            m.randomize(-1., 2.); m1.randomize(-1., 3.); m2.randomize(0., 2.);
            ms.randomize(-3., 1.);

            rmatrix mr1 = ms * m1 * alpha + m * beta;
            CheckReal ((mr1 - m.symm (true, ms, m1, alpha, beta)).norm(), (treal) 0., "rmatrix::symm" , os, __LINE__, dPessimisticSp);

            m.resize(4,3);
            rmatrix mr2 = m2 * ms * alpha + m * beta;
            CheckReal ((mr2 - m.symm (false, ms, m2, alpha, beta)).norm(), (treal) 0., "rmatrix::symm" , os, __LINE__, dPessimisticSp);
        }
        {
            treal alpha = 1.3;
            rmatrix m(3,4);
            rvector vc(3);
            rvector vr(4);
            m.randomize(-1., 2.); vc.randomize(-1., 3.); vr.randomize(0., 2.);

            rmatrix mr = m + vc.rank1update (vr) * alpha;
            CheckReal ((mr - m.ger(alpha, vc, vr)).norm(), (treal) 0., "rmatrix::ger" , os, __LINE__, dPessimisticSp);
        }
        {
            treal a[] = {1., 2., 3., 4., 5., 6., 7., 8., 9.};
            rmatrix  m(2,3);
            const srmatrix ms(a,3);

            m.diag(-1).set(1.);
            m.diag(0).set(2.);
            m.diag(1).set(3.);
            m.diag(2).set(4.);
            CheckReal (m(1,2), (treal) 3., "rmatrix::diag" , os, __LINE__, dPessimisticSp);
            CheckReal (m(1,3), (treal) 4., "rmatrix::diag" , os, __LINE__, dPessimisticSp);
            CheckReal (m(2,1), (treal) 1., "rmatrix::diag" , os, __LINE__, dPessimisticSp);
            CheckReal (ms.diag(0)(2), (treal) 5., "rmatrix::diag" , os, __LINE__, dPessimisticSp);
        }
        {
            std::complex<treal> alpha = std::complex<treal>(1.3,-0.7);
            std::complex<treal> beta  = std::complex<treal>(0.15,-1.09);
            scbmatrix m(3,1,0);
            cvector c(3);
            cvector v(3);
            m.randomize_real(-1., 2.); 
            m.randomize_imag(0., 1.); 
            v.randomize_real(-1., 3.); 
            v.randomize_imag(2., 4.); 
            c.randomize_real(0., 2.);
            c.randomize_imag(3., 7.);

            cvector vr1 = m * v * alpha + c * beta;
            CheckReal ((vr1 - c.gbmv(false, m, alpha, v, beta)).norm(), (treal) 0., "cvector::gbmv" , os, __LINE__, dPessimisticSp);
            cvector vr2 = c * m * alpha + v * beta;
            CheckReal ((vr2 - v.gbmv(true, m, alpha, c, beta)).norm(), (treal) 0., "cvector::gbmv" , os, __LINE__, dPessimisticSp);
        }
        {
            std::complex<treal> alpha = std::complex<treal>(1.3,-0.7);
            std::complex<treal> beta  = std::complex<treal>(0.15,-1.09);
            cmatrix m(3,2);
            cvector c(3);
            cvector v(2);
            m.randomize_real(-1., 2.); 
            m.randomize_imag(0., 1.); 
            v.randomize_real(-1., 3.); 
            v.randomize_imag(2., 4.); 
            c.randomize_real(0., 2.);
            c.randomize_imag(3., 7.);

            cvector vr1 = m * v * alpha + c * beta;
            CheckReal ((vr1 - c.gemv(false, m, alpha, v, beta)).norm(), (treal) 0., "cvector::gemv" , os, __LINE__, dPessimisticSp);
            cvector vr2 = c * m * alpha + v * beta;
            CheckReal ((vr2 - v.gemv(true, m, alpha, c, beta)).norm(), (treal) 0., "cvector::gemv" , os, __LINE__, dPessimisticSp);
        }
        {
            cvector vc(3);
            vc.set(std::complex<treal>(1.,1.));
            vc.imag()(1) = 7.77;
            CheckComplex (vc[1], tcomplex ((treal) 1., (treal) 7.77), "cvector::imag" , os, __LINE__, dPessimisticSp);
            CheckComplex (vc[2], tcomplex ((treal) 1., (treal) 1.), "cvector::imag" , os, __LINE__, dPessimisticSp);
        }
        {
            cvector vc(3);
            vc.set(std::complex<treal>(1.,1.));
            vc.real()(1) = 7.77;
            CheckComplex (vc[1], tcomplex ((treal) 7.77, (treal) 1.), "cvector::real" , os, __LINE__, dPessimisticSp);
            CheckComplex (vc[2], tcomplex ((treal) 1., (treal) 1.), "cvector::real" , os, __LINE__, dPessimisticSp);
        }
        {
            cvector v(3);
            v.set_real(1.);
            CheckComplex (v[2], tcomplex ((treal) 1., (treal) 0.), "cvector::set_real" , os, __LINE__, dPessimisticSp);
        }
        {
            rvector v(3);
            cvector vc(3);
            v(1) = 1.; v(2) = 2.; v(3) = 3.;
            vc.assign_imag(v);
            CheckComplex (vc[2], tcomplex ((treal) 0., (treal) 2.), "cvector::assign_imag" , os, __LINE__, dPessimisticSp);
        }
        {
            treal alpha = 1.3;
            treal beta = -0.7;
            srbmatrix m(3, 1, 0);
            rvector c(3);
            rvector v(3);
            m.randomize(-1., 2.); v.randomize(-1., 3.); c.randomize(0., 2.);

            rvector vr1 = m * v * alpha + c * beta;
            CheckReal ((vr1 - c.gbmv(false, m, alpha, v, beta)).norm(), (treal) 0., "rvector::gbmv" , os, __LINE__, dPessimisticSp);
            rvector vr2 = c * m * alpha + v * beta;
            CheckReal ((vr2 - v.gbmv(true, m, alpha, c, beta)).norm(), (treal) 0., "rvector::gbmv" , os, __LINE__, dPessimisticSp);
        }
        {
            treal alpha = 1.3;
            treal beta = -0.7;
            rmatrix m(4,3);
            rvector c(4);
            rvector v(3);
            m.randomize(-1., 2.); v.randomize(-1., 3.); c.randomize(0., 2.);

            rvector vr1 = m * v * alpha + c * beta;
            CheckReal ((vr1 - c.gemv(false, m, alpha, v, beta)).norm(), (treal) 0., "rvector::gemv" , os, __LINE__, dPessimisticSp);
            rvector vr2 = c * m * alpha + v * beta;
            CheckReal ((vr2 - v.gemv(true, m, alpha, c, beta)).norm(), (treal) 0., "rvector::gemv" , os, __LINE__, dPessimisticSp);
        }
        {
            srsmatrix m(3);
            srmatrix me(3);
            rvector v(3);
            m.randomize(1., 3.);

            v.eig (m, me);
            CheckReal ((m * me(1) - me(1) * v(1)).norm(), (treal) 0., "srsmatrix::eig" , os, __LINE__, dPessimisticSp);
            CheckReal ((m * me(2) - me(2) * v(2)).norm(), (treal) 0., "srsmatrix::eig" , os, __LINE__, dPessimisticSp);
            CheckReal ((m * me(3) - me(3) * v(3)).norm(), (treal) 0., "srsmatrix::eig" , os, __LINE__, dPessimisticSp);

            CheckReal (me(1) * me(2), (treal) 0., "srsmatrix::eig" , os, __LINE__, dPessimisticSp);
            CheckReal (me(1) * me(3), (treal) 0., "srsmatrix::eig" , os, __LINE__, dPessimisticSp);
            CheckReal (me(3) * me(2), (treal) 0., "srsmatrix::eig" , os, __LINE__, dPessimisticSp);

            schmatrix mc(3);
            scmatrix mce(3);
            mc.randomize_real(1., 3.);
            mc.randomize_imag(1., 3.);

            v.eig (mc, mce);

            CheckReal ((mc * mce(1) - mce(1) * v(1)).norm(), (treal) 0., "schmatrix::eig" , os, __LINE__, dPessimisticSp);
            CheckReal ((mc * mce(2) - mce(2) * v(2)).norm(), (treal) 0., "schmatrix::eig" , os, __LINE__, dPessimisticSp);
            CheckReal ((mc * mce(3) - mce(3) * v(3)).norm(), (treal) 0., "schmatrix::eig" , os, __LINE__, dPessimisticSp);

            CheckComplex (mce(1) % mce(2), tcomplex ((treal) 0., (treal) 0.), "srsmatrix::eig" , os, __LINE__, dPessimisticSp);
            CheckComplex (mce(1) % mce(3), tcomplex ((treal) 0., (treal) 0.), "srsmatrix::eig" , os, __LINE__, dPessimisticSp);
            CheckComplex (mce(3) % mce(2), tcomplex ((treal) 0., (treal) 0.), "srsmatrix::eig" , os, __LINE__, dPessimisticSp);
        }
        {
            treal m[] = {1., -1., 1., 2., -2., 1., 3., -2., 1.};
            treal b1[] = {1., 2., 3.};
            treal b2[] = {0., -1., -2.};
            srmatrix ma(m, 3);
            srmatrix mLU(3);
            rvector  vb1(b1, 3);
            rvector  vb2(b2, 3);
            rvector  vx1(3);
            rvector  vx2(3);
            iarray   nPivots(3);
            treal   dErr = 0.;

            mLU.low_up(ma, nPivots);
            vx1.solve_lu (ma, mLU, nPivots, vb1, dErr);
            CheckReal (dErr, (treal) 0., "rmatrix::solve_lu" , os, __LINE__, dPessimisticSp);
            vx2.solve_lu (ma, mLU, nPivots, vb2);
            CheckReal ((ma * vx1 - vb1).norm(), (treal) 0., "rmatrix::solve_lu" , os, __LINE__, dPessimisticSp);
            CheckReal ((ma * vx2 - vb2).norm(), (treal) 0., "rmatrix::solve_lu" , os, __LINE__, dPessimisticSp);
        }
        {
            scmatrix ma(3);
            scmatrix mLU(3);
            cmatrix  mb1(3,2);
            cmatrix  mb2(3,2);
            cmatrix  mx1(3,2);
            cmatrix  mx2(3,2);
            iarray   nPivots(3);
            treal   dErr = 0.;
            ma.randomize_real(0.,10.); ma.randomize_imag(0.,10.);
            mb1.randomize_real(0.,10.); mb1.randomize_imag(0.,10.);
            mb2.randomize_real(0.,10.); mb2.randomize_imag(0.,10.);

            mLU.low_up(ma, nPivots);
            mx1.solve_lu (ma, mLU, nPivots, mb1, dErr);
            CheckReal (dErr, (treal) 0., "cmatrix::solve_lu" , os, __LINE__, dPessimisticSp);
            mx2.solve_lu (ma, mLU, nPivots, mb2);
            CheckReal ((ma * mx1 - mb1).norm(), (treal) 0., "cmatrix::solve_lu" , os, __LINE__, dPessimisticSp);
            CheckReal ((ma * mx2 - mb2).norm(), (treal) 0., "cmatrix::solve_lu" , os, __LINE__, dPessimisticSp);
        }
        {
            rvector v(5);
            v.set(3.);
            CheckReal (v[3], (treal) 3., "rvector::set" , os, __LINE__);
        }
        {
            const treal a[] = {1., 2., 3., 4., 5., 6., 7.,};
            rvector v (5);
            rvector v2 (4);

            v.assign(a);
            CheckReal (v[3], (treal) 3., "rvector::assign" , os, __LINE__);
            v2.assign(a, 2);
            CheckReal (v2[3], (treal) 5., "rvector::assign" , os, __LINE__);
        }
        {
            rvector v (5, 1.5);
            CheckReal (v[3], (treal) 1.5, "rvector (int, treal)" , os, __LINE__);
        }
        {
            rmatrix m (100, 200);
            srmatrix ms (m, 30, 40, 5); // 5x5 submatrix
            CheckInt (ms.ld(), 100, "srmatrix::ld" , os, __LINE__);
        }
        {
            iarray a(5);
            iarray::iterator pos = a.begin() + 2;
            a.insert(pos, 88);
            CheckInt (a[3], 88, "iarray::begin, iarray::insert" , os, __LINE__);
            pos = a.begin() + 1;
            a.erase(pos);
            CheckInt (a[2], 88, "iarray::begin, iarray::erase" , os, __LINE__);
            CheckInt (a[3], 0, "iarray::begin, iarray::erase" , os, __LINE__);
        }
        {
            iarray a(5);
            a.push_back(88);
            CheckInt (a[6], 88, "iarray::push_back" , os, __LINE__);
            a.pop_back();
            CheckInt (a[5], 0, "iarray::pop_back" , os, __LINE__);
            CheckInt (a.size(), 5, "iarray::pop_back" , os, __LINE__);
        }
        {
            iarray a(5);
            a[1] = 1; a[2] = 2; a[3] = 3; a[4] = 4; a[5] = 5;
            CheckInt (a.at(0), 1, "iarray::at" , os, __LINE__);
            CheckInt (a.at(4), 5, "iarray::at" , os, __LINE__);
        }
        {
            iarray a(5);
            a[1] = 1; a[2] = 2; a[3] = 3; a[4] = 4; a[5] = 5;

            int val = 5;
            for (iarray::reverse_iterator it = a.rbegin(); it != a.rend(); ++it)
            {
                CheckInt (*it, val, "iarray::reverse_iterator" , os, __LINE__);
                --val;
            }

            CheckInt (a.front(), 1, "iarray::front" , os, __LINE__);
            CheckInt (a.back(), 5, "iarray::back" , os, __LINE__);
        }

        {
            const int a[] = {1, 2, 3, 4};
            iarray v (a, 3);
            v.resize(2);
            CheckInt (v[2], 2, "iarray.resize" , os, __LINE__);
            v.resize(4);
            CheckInt (v[4], 0, "iarray.resize" , os, __LINE__);
        }
        {
            iarray a(5);
            a.set(3);
            iarray b(a);
            CheckInt (b[4], 3, "iarray copy ctr" , os, __LINE__);
        }
        {
            iarray a(5), b(5);
            a.set(3);
            b = a;
            CheckInt (b[4], 3, "iarray assignment" , os, __LINE__);
        }
        {
            iarray a(5);
            a.set(3);
            CheckInt (a[4], 3, "iarray.set" , os, __LINE__);
        }
        {
            iarray a;
            a.resize(10);
            CheckInt (a.size(), 10, "iarray.resize" , os, __LINE__);
        }

        {
            const int a[] = {1, 2, 3, 4};
            const iarray v (a+1, a+3);
            CheckInt (v[2], 3, "iarray (*,*)" , os, __LINE__);
        }

        {
            iarray a(10);
            a[2] = 1;
            CheckInt (a.get()[1], 1, "iarray.get", os, __LINE__);
        }

        {
            srsmatrix ssm1(4);
            srsmatrix ssm2(4);

            rmatrix m1(4,4);
            rmatrix m2(4,4);
            rmatrix m3(4,4);

            ssm1.set(1.);
            m1.set(1.);

            m2 = ssm1 * m1;
            CheckReal (m2.norminf(),  (treal) 16. ,  "srsmatrix * rmatrix",  os, __LINE__);

            m2 = m1 + ssm1;
            CheckReal (m2(3,4),  (treal) 2. ,  "srsmatrix + rmatrix",  os, __LINE__);

            m2 = m1 * ssm1;
            CheckReal (m2.norminf(),  (treal) 16. ,  "rmatrix * srsmatrix",  os, __LINE__);
        }




        rvector vs1(5);
        vs1[1] = 1.; vs1[2] = 2.; vs1[3] = 3.; vs1[4] = 4.; vs1[5] = 5.;

        rvector::iterator it = vs1.begin() + 1;
        rvector::iterator ite = vs1.erase(it);

        CheckReal (vs1[1],  (treal) 1. ,  "rvector.insert",  os, __LINE__);
        CheckReal (vs1[2],  (treal) 3. ,  "rvector.insert",  os, __LINE__);
        CheckReal (vs1[3],  (treal) 4. ,  "rvector.insert",  os, __LINE__);

        ite = vs1.insert(ite, 10.);

        CheckReal (vs1[1],  (treal) 1. ,  "rvector.insert",  os, __LINE__);
        CheckReal (vs1[2],  (treal) 10. ,  "rvector.insert",  os, __LINE__);
        CheckReal (vs1[3],  (treal) 3. ,  "rvector.insert",  os, __LINE__);

        vs1.push_back(9.);
        CheckReal (vs1[5],  (treal) 5. ,  "rvector.push_back",  os, __LINE__);
        CheckReal (vs1[6],  (treal) 9. ,  "rvector.push_back",  os, __LINE__);
        CheckReal (*std::max_element(vs1.begin(), vs1.end()),  (treal) 10. ,  "rvector.max_element",  os, __LINE__);

        std::sort(vs1.begin(), vs1.end());
        CheckReal (vs1[6],  (treal) 10. ,  "std::sort",  os, __LINE__);

        std::reverse(vs1.begin(), vs1.end());
        CheckReal (vs1[1],  (treal) 10. ,  "std::reverse",  os, __LINE__);
        CheckReal (vs1[2],  (treal) 9. ,  "std::reverse",  os, __LINE__);
        CheckReal (vs1[6],  (treal) 1. ,  "std::reverse",  os, __LINE__);

        {
            // N > 64*M bug fixed
            rmatrix A(600,4);
            srmatrix U(600), V(4);
            const rvector singVal = A.svd(U,V);

            rmatrix singValM (A);
            singValM.set(0.);
            singValM.diag(0) = singVal;

            CheckReal    ((A * ~V - U * singValM).norm(),  (treal) 0.,  "rmatrix svd", os, __LINE__, dPessimisticSp);
            CheckReal    ((~A * U - ~(singValM * V)).norm(),  (treal) 0.,  "rmatrix svd", os, __LINE__, dPessimisticSp);
        }

        {
            // Gantmaher, p. 33
            rmatrix mA(3,4);
            mA(1,1) =  1.; mA(1,2) = -1.; mA(1,3) =  2.; mA(1,4) =  0.;
            mA(2,1) = -1.; mA(2,2) =  2.; mA(2,3) = -3.; mA(2,4) =  1.;
            mA(3,1) =  0.; mA(3,2) =  1.; mA(3,3) = -1.; mA(3,4) =  1.;

            // lower rank case
            rmatrix mX = mA.pinv(dPessimisticSp);
            CheckReal ((mA * mX * mA - mA).norm2(), (treal) 0.,  "pinv, lower rank, m < n", os, __LINE__, dPessimisticSp);

            // m > n
            mA.transpose();
            rmatrix mX2 = mA.pinv(dPessimisticSp);
            CheckReal ((mA * mX2 * mA - mA).norm2(), (treal) 0.,  "pinv, lower rank, m > n", os, __LINE__, dPessimisticSp);

            // full rank case
            mA.transpose();
            mA(1,3) = 4.;
            mX.pinv (mA);
            CheckReal ((mA * mX * mA - mA).norm2(), (treal) 0.,  "pinv, full rank, m < n", os, __LINE__, dPessimisticSp);

            // m > n
            mA.transpose();
            mX2.pinv(mA);
            CheckReal ((mA * mX2 * mA - mA).norm2(), (treal) 0.,  "pinv, full rank, m > n", os, __LINE__, dPessimisticSp);
        }

        {
            cmatrix mA(3,4), mX(4,3);
            mA.randomize_real(-2., 11.);
            mA.randomize_imag(-9., 7.);

            mX.pinv (mA);
            CheckReal ((mA * mX * mA - mA).norm2(), (treal) 0.,  "complex pinv, m < n", os, __LINE__, dPessimisticSp);

            // m > n
            mA.conj();
            cmatrix mX2 = mA.pinv();
            CheckReal ((mA * mX2 * mA - mA).norm2(), (treal) 0.,  "complex pinv, m > n", os, __LINE__, dPessimisticSp);
        }

        {
            srbmatrix mA (40, 1, 2);
            mA.diag(0).randomize(-1.,1.);
            mA.diag(-1).randomize(-1.,1.);
            mA.diag(1).randomize(-1.,1.);

            rmatrix mX = mA.pinv(dPessimisticSp);
            CheckReal ((mA * mX * mA - mA).norm2(), (treal) 0.,  "srbmatrix pinv", os, __LINE__, dVeryPessimisticSp);

            mA.transpose();
            mX.pinv (mA);
            CheckReal ((mA * mX * mA - mA).norm2(), (treal) 0.,  "srbmatrix pinv", os, __LINE__, dVeryPessimisticSp);
        }

        {
            scbmatrix mA (40, 1, 2);
            mA.diag(0).randomize_real(-1.,1.);
            mA.diag(-1).randomize_real(-1.,1.);
            mA.diag(1).randomize_real(-1.,1.);
            mA.diag(0).randomize_imag(-1.,2.);
            mA.diag(-1).randomize_imag(-1.,2.);
            mA.diag(1).randomize_imag(-1.,2.);

            scmatrix mX = mA.pinv(dPessimisticSp);
            CheckReal ((mA * mX * mA - mA).norm2(), (treal) 0.,  "scbmatrix pinv", os, __LINE__, dVeryPessimisticSp);

            mA.conj();
            mX.pinv (mA);
            CheckReal ((mA * mX * mA - mA).norm2(), (treal) 0.,  "scbmatrix pinv", os, __LINE__, dPessimisticSp);
        }

        // 5.4.1
        {
            treal m[] = {1., -1., 1., 2., -2., 1.,
                          3., -2., 1., 0., -2., 1.};
            rmatrix mA(m,4,3);
            rmatrix mSigma(4,3);
            rvector v(3);
            srmatrix mU(4), mVH(3);
            v.svd(mA, mU, mVH);
            mSigma.diag(0) = v;
//            std::cout << mU << std::endl;
//            std::cout << mVH << std::endl;
//            std::cout << mSigma << std::endl;
//            std::cout << (mA * ~mVH - mU * mSigma).norm() << std::endl;
//            std::cout << (~mA * mU - ~(mSigma * mVH)).norm() << std::endl;
            CheckReal ((mA * ~mVH - mU * mSigma).norm(), 0.,  "rmatrix svd",  os, __LINE__, dPessimisticSp);
            CheckReal ((~mA * mU - ~(mSigma * mVH)).norm(), 0.,  "rmatrix svd",  os, __LINE__, dPessimisticSp);
        }

        {
            treal m[] = {1., -1., 1., 2., -2., 1.,
                          3., -2., 1., 0., -2., 1.};
            cmatrix mA((std::complex<treal>*) m, 2, 3);
            cmatrix mSigma(2,3);
            rvector v(2);
            scmatrix mU(2), mVH(3);
            
            v = mA.svd(mU, mVH);
            mSigma.diag(0) = cvector(v);

//            std::cout << mU << std::endl;
//            std::cout << mVH << std::endl;
//            std::cout << mSigma << std::endl;
//            std::cout << (mA * ~mVH - mU * mSigma).norm() << std::endl;
//            std::cout << (~mA * mU - ~(mSigma * mVH)).norm() << std::endl;
            CheckReal ((mA * ~mVH - mU * mSigma).norm(), 0.,  "rmatrix svd",  os, __LINE__, dPessimisticSp);
            CheckReal ((~mA * mU - ~(mSigma * mVH)).norm(), 0.,  "rmatrix svd",  os, __LINE__, dPessimisticSp);
        }


        {
            cvm::srmatrix s(9);
            cvm::srbmatrix m(3,0,1);
            for (int i = 1; i <= 9; ++i) {
                s[i].set((treal)i);
            }
            m.assign(s(9));     // should be 1,2,..9
            CheckReal (m(1,1), (treal) 2., "rmatrix.assign(vector)",  os, __LINE__);
            CheckReal (m(1,2), (treal) 3., "rmatrix.assign(vector)",  os, __LINE__);
            CheckReal (m(3,3), (treal) 6., "rmatrix.assign(vector)",  os, __LINE__);
            m.assign(s[9]);     // should be 9,9,..9
            CheckReal (m(1,1), (treal) 9., "rmatrix.assign(vector)",  os, __LINE__);
            CheckReal (m(1,2), (treal) 9., "rmatrix.assign(vector)",  os, __LINE__);
            CheckReal (m(3,3), (treal) 9., "rmatrix.assign(vector)",  os, __LINE__);
        }

        {
            cvm::scmatrix s(9);
            cvm::scbmatrix m(3,1,1);
            for (int i = 1; i <= 9; ++i) {
                s[i].set(tcomplex((treal)i,(treal)-i));
            }
            m.assign(s(9));     // should be 1,2,..9
            CheckComplex (m(1,1), tcomplex((treal) 2., (treal) -2.), "scbmatrix.assign(vector)",  os, __LINE__);
            CheckComplex (m(1,2), tcomplex((treal) 4., (treal) -4.), "scbmatrix.assign(vector)",  os, __LINE__);
            CheckComplex (m(3,3), tcomplex((treal) 8., (treal) -8.), "scbmatrix.assign(vector)",  os, __LINE__);
            CheckComplex (m(1,3), tcomplex((treal) 0., (treal) 0.),  "scbmatrix.assign(vector)",  os, __LINE__);
            m.assign(s[9]);     // should be 9,9,..9
            CheckComplex (m(1,1), tcomplex((treal) 9., (treal) -9.), "scbmatrix.assign(vector)",  os, __LINE__);
            CheckComplex (m(1,2), tcomplex((treal) 9., (treal) -9.), "scbmatrix.assign(vector)",  os, __LINE__);
            CheckComplex (m(3,3), tcomplex((treal) 9., (treal) -9.), "scbmatrix.assign(vector)",  os, __LINE__);
            CheckComplex (m(1,3), tcomplex((treal) 0., (treal) 0.),  "scbmatrix.assign(vector)",  os, __LINE__);
        }


        {
            cvm::srmatrix s(9);
            cvm::rmatrix mbig(30,30);
            cvm::rmatrix m(mbig,4,7,3,3);
            for (int i = 1; i <= 9; ++i) {
                s[i].set((treal)i);
            }
            m.assign(s(9));     // should be 1,2,..9
            CheckReal (m(1,1), (treal) 1., "rmatrix.assign(vector)",  os, __LINE__);
            CheckReal (m(1,2), (treal) 4., "rmatrix.assign(vector)",  os, __LINE__);
            CheckReal (m(3,3), (treal) 9., "rmatrix.assign(vector)",  os, __LINE__);
            m.assign(s[9]);     // should be 9,9,..9
            CheckReal (m(1,1), (treal) 9., "rmatrix.assign(vector)",  os, __LINE__);
            CheckReal (m(1,2), (treal) 9., "rmatrix.assign(vector)",  os, __LINE__);
            CheckReal (m(3,3), (treal) 9., "rmatrix.assign(vector)",  os, __LINE__);
        }

        {
            cvm::srmatrix s(9);
            cvm::srmatrix mbig(20);
            cvm::srmatrix m(mbig,4,7,3);
            for (int i = 1; i <= 9; ++i) {
                s[i].set((treal)i);
            }
            m.assign(s(9));     // should be 1,2,..9
            CheckReal (m(1,1), (treal) 1., "srmatrix.assign(vector)",  os, __LINE__);
            CheckReal (m(1,2), (treal) 4., "srmatrix.assign(vector)",  os, __LINE__);
            CheckReal (m(3,3), (treal) 9., "srmatrix.assign(vector)",  os, __LINE__);
            m.assign(s[9]);     // should be 9,9,..9
            CheckReal (m(1,1), (treal) 9., "srmatrix.assign(vector)",  os, __LINE__);
            CheckReal (m(1,2), (treal) 9., "srmatrix.assign(vector)",  os, __LINE__);
            CheckReal (m(3,3), (treal) 9., "srmatrix.assign(vector)",  os, __LINE__);
        }

        {
            cvm::scmatrix s(9);
            cvm::cmatrix mbig(30,30);
            cvm::cmatrix m(mbig,4,7,3,3);
            for (int i = 1; i <= 9; ++i) {
                s[i].set(tcomplex((treal)i,(treal)-i));
            }
            m.assign(s(9));     // should be 1,2,..9
            CheckComplex (m(1,1), tcomplex((treal) 1., (treal) -1.), "cmatrix.assign(vector)",  os, __LINE__);
            CheckComplex (m(1,2), tcomplex((treal) 4., (treal) -4.), "cmatrix.assign(vector)",  os, __LINE__);
            CheckComplex (m(3,3), tcomplex((treal) 9., (treal) -9.), "cmatrix.assign(vector)",  os, __LINE__);
            m.assign(s[9]);     // should be 9,9,..9
            CheckComplex (m(1,1), tcomplex((treal) 9., (treal) -9.), "cmatrix.assign(vector)",  os, __LINE__);
            CheckComplex (m(1,2), tcomplex((treal) 9., (treal) -9.), "cmatrix.assign(vector)",  os, __LINE__);
            CheckComplex (m(3,3), tcomplex((treal) 9., (treal) -9.), "cmatrix.assign(vector)",  os, __LINE__);
        }

        {
            cvm::scmatrix s(9);
            cvm::scmatrix mbig(20);
            cvm::scmatrix m(mbig,4,7,3);
            for (int i = 1; i <= 9; ++i) {
                s[i].set(tcomplex((treal)i,(treal)-i));
            }
            m.assign(s(9));     // should be 1,2,..9
            CheckComplex (m(1,1), tcomplex((treal) 1., (treal) -1.), "scmatrix.assign(vector)",  os, __LINE__);
            CheckComplex (m(1,2), tcomplex((treal) 4., (treal) -4.), "scmatrix.assign(vector)",  os, __LINE__);
            CheckComplex (m(3,3), tcomplex((treal) 9., (treal) -9.), "scmatrix.assign(vector)",  os, __LINE__);
            m.assign(s[9]);     // should be 9,9,..9
            CheckComplex (m(1,1), tcomplex((treal) 9., (treal) -9.), "scmatrix.assign(vector)",  os, __LINE__);
            CheckComplex (m(1,2), tcomplex((treal) 9., (treal) -9.), "scmatrix.assign(vector)",  os, __LINE__);
            CheckComplex (m(3,3), tcomplex((treal) 9., (treal) -9.), "scmatrix.assign(vector)",  os, __LINE__);
        }

        {
            cvm::srmatrix s(9);
            cvm::srsmatrix m(3);
            s(9,2) = (treal) 1.;
            s(9,4) = (treal) 1.;
            s(9,9) = (treal) 5.;
            m.assign(s[9]);
            CheckReal (m(1,1), (treal) 0., "srsmatrix.assign(vector)",  os, __LINE__);
            CheckReal (m(1,2), (treal) 1., "srsmatrix.assign(vector)",  os, __LINE__);
            CheckReal (m(2,1), (treal) 1., "srsmatrix.assign(vector)",  os, __LINE__);
            CheckReal (m(3,3), (treal) 5., "srsmatrix.assign(vector)",  os, __LINE__);
        }


        // 5.4.2
        {
            // bug fix check
            cmatrix a(3,4);
            CheckInt (a.ld(), 3, "a.ld()", os, __LINE__);
            a.resize(0,0);
            CheckInt (a.ld(), 0, "a.ld()", os, __LINE__);
        }
        {
            std::vector<cvector> vcv;
            vcv.reserve(5);
            vcv.push_back(cvector(10));
            vcv.push_back(cvector());
            vcv[0][1] = tcomplex((treal) 1., (treal) -1.);
            CheckComplex (vcv[0](1), tcomplex((treal) 1., (treal) -1.), "std::vector<cvector>[][]",  os, __LINE__);
        }
        {
            std::vector<cmatrix> vcm;
            vcm.reserve(5);
            vcm.push_back(cmatrix(10,20));
            vcm[0][1][2] = tcomplex((treal) 1., (treal) -1.);
            CheckComplex (vcm[0](1,2), tcomplex((treal) 1., (treal) -1.), "std::vector<cmatrix>[][][]",  os, __LINE__);
        }
        {
            std::vector<rmatrix> vcm;
            vcm.reserve(5);
            vcm.push_back(srmatrix(10));
            vcm.push_back(srmatrix());
            vcm[0][1][2] = (treal) 7.77;
            CheckReal (vcm[0](1,2), (treal) 7.77, "std::vector<srmatrix>[][][]",  os, __LINE__);
        }

        // 5.5 QR stuff
        {
            treal a[] = {1., 2., 3., 4., 5., 6.};
            const cvm::rmatrix mh(a, 2, 3);
            const cvm::rmatrix mv(a, 3, 2);
            cvm::srmatrix s2(2), s3(3);
            cvm::rmatrix  h(2,3), v(3,2);

            mh.qr(h,s3);
//            std::cout << (eye_real(2) - ~rmatrix(h,1,1,2,2) * rmatrix(h,1,1,2,2)).norm()
//                      << " " << (mh - h * s3).norm() << std::endl;
            CheckReal ((eye_real(2) - ~rmatrix(h,1,1,2,2) * rmatrix(h,1,1,2,2)).norm(), (treal) 0., "cvm::rmatrix QR",  os, __LINE__, dPessimisticSp);
            CheckReal ((mh - h * s3).norm(), (treal) 0., "cvm::rmatrix QR",  os, __LINE__, dPessimisticSp);
            mh.qr(s2,h);
//            std::cout << (eye_real(2) - ~s2 * s2).norm()
//                      << " " << (mh - s2 * h).norm() << std::endl;
            CheckReal ((eye_real(2) - ~s2 * s2).norm(), (treal) 0., "cvm::rmatrix QR",  os, __LINE__, dPessimisticSp);
            CheckReal ((mh - s2 * h).norm(), (treal) 0., "cvm::rmatrix QR",  os, __LINE__, dPessimisticSp);
            mv.qr(v,s2);
//            std::cout << (eye_real(2) - ~v * v).norm()
//                      << " " << (mv - v * s2).norm() << std::endl;
            CheckReal ((eye_real(2) - ~v * v).norm(), (treal) 0., "cvm::rmatrix QR",  os, __LINE__, dPessimisticSp);
            CheckReal ((mv - v * s2).norm(), (treal) 0., "cvm::rmatrix QR",  os, __LINE__, dPessimisticSp);
            mv.qr(s3,v);
//            std::cout << (eye_real(3) - ~s3 * s3).norm()
//                      << " " << (mv - s3 * v).norm() << std::endl;
            CheckReal ((eye_real(3) - ~s3 * s3).norm(), (treal) 0., "cvm::rmatrix QR",  os, __LINE__, dPessimisticSp);
            CheckReal ((mv - s3 * v).norm(), (treal) 0., "cvm::rmatrix QR",  os, __LINE__, dPessimisticSp);
        }
        {
            treal a[] = {1., 2., 3., 4., 5., 6., 7., 8., 9.};
            const cvm::srmatrix m(a, 3);
            cvm::srmatrix q(3), r(3);

            m.qr(q,r);
//            std::cout << (eye_real(3) - ~q * q).norm()
//                      << " " << (m - q * r).norm() << std::endl;
            CheckReal ((eye_real(3) - ~q * q).norm(), (treal) 0., "cvm::srmatrix QR",  os, __LINE__, dPessimisticSp);
            CheckReal ((m - q * r).norm(), (treal) 0., "cvm::srmatrix QR",  os, __LINE__, dPessimisticSp);
        }
        {
            treal ar[] = {1., 2., 3., 4., 5., 6.};
            treal ai[] = {1., -1., 2., -2., 3., -3.};
            const cvm::cmatrix mh(ar, ai, 2, 3);
            const cvm::cmatrix mv(ar, ai, 3, 2);
            cvm::scmatrix s2(2), s3(3);
            cvm::cmatrix  h(2,3), v(3,2);

            mh.qr(h,s3);
//            std::cout << (eye_complex(2) - ~cmatrix(h,1,1,2,2) * cmatrix(h,1,1,2,2)).norm()
//                      << " " << (mh - h * s3).norm() << std::endl;
            CheckReal ((eye_complex(2) - ~cmatrix(h,1,1,2,2) * cmatrix(h,1,1,2,2)).norm(), (treal) 0., "cvm::cmatrix QR",  os, __LINE__, dPessimisticSp);
            CheckReal ((mh - h * s3).norm(), (treal) 0., "cvm::cmatrix QR",  os, __LINE__, dPessimisticSp);
            mh.qr(s2,h);
//            std::cout << (eye_complex(2) - ~s2 * s2).norm()
//                      << " " << (mh - s2 * h).norm() << std::endl;
            CheckReal ((eye_complex(2) - ~s2 * s2).norm(), (treal) 0., "cvm::cmatrix QR",  os, __LINE__, dPessimisticSp);
            CheckReal ((mh - s2 * h).norm(), (treal) 0., "cvm::cmatrix QR",  os, __LINE__, dPessimisticSp);
            mv.qr(v,s2);
//            std::cout << (eye_complex(2) - ~v * v).norm()
//                      << " " << (mv - v * s2).norm() << std::endl;
            CheckReal ((eye_complex(2) - ~v * v).norm(), (treal) 0., "cvm::cmatrix QR",  os, __LINE__, dPessimisticSp);
            CheckReal ((mv - v * s2).norm(), (treal) 0., "cvm::cmatrix QR",  os, __LINE__, dPessimisticSp);
            mv.qr(s3,v);
//            std::cout << (eye_complex(3) - ~s3 * s3).norm()
//                      << " " << (mv - s3 * v).norm() << std::endl;
            CheckReal ((eye_complex(3) - ~s3 * s3).norm(), (treal) 0., "cvm::cmatrix QR",  os, __LINE__, dPessimisticSp);
            CheckReal ((mv - s3 * v).norm(), (treal) 0., "cvm::cmatrix QR",  os, __LINE__, dPessimisticSp);
        }
        {
            treal ar[] = {1., 2., 3., 4., 5., 6., 7., 8., 9.};
            treal ai[] = {1., -1., 2., -2., 3., -3., 4., -4., 5.};
            const cvm::scmatrix m(ar, ai, 3);
            cvm::scmatrix q(3), r(3);

            m.qr(q,r);
//            std::cout << (eye_complex(3) - ~q * q).norm()
//                      << " " << (m - q * r).norm() << std::endl;
            CheckReal ((eye_complex(3) - ~q * q).norm(), (treal) 0., "cvm::scmatrix QR",  os, __LINE__, dPessimisticSp);
            CheckReal ((m - q * r).norm(), (treal) 0., "cvm::scmatrix QR",  os, __LINE__, dPessimisticSp);
        }


        {
            treal a[] = {1., 4., 7., 2., 5., 8., 3., 6., 0.};
            cvm::srmatrix A(a, 3);

            cvm::srmatrix Q(3);
            cvm::srmatrix R(3);

            A.qr(Q,R);

            CheckReal ((A - Q * R).norm(), (treal) 0., "cvm::srmatrix QR",  os, __LINE__, dPessimisticSp);
            CheckReal ((eye_real(3) - ~Q * Q).norm(), (treal) 0., "cvm::srmatrix QR - Q",  os, __LINE__, dPessimisticSp);
        }
        {
            const int m = 10;
            srbmatrix A (m, 2, 3);
            A.randomize(-10., 10.);

            cvm::srmatrix Q(m);
            cvm::srmatrix R(m);

            A.qr(Q,R);

            CheckReal ((Q * R - A).norm(), (treal) 0., "cvm::srbmatrix QR",  os, __LINE__, dPessimisticSp);
            CheckReal ((eye_real(m) - ~Q * Q).norm(), (treal) 0., "cvm::srbmatrix QR - Q",  os, __LINE__, dPessimisticSp);
        }
        {
            const int m = 10;
            srsmatrix A (m);
            A.randomize(-10., 10.);

            cvm::srmatrix Q(m);
            cvm::srmatrix R(m);

            A.qr(Q,R);

            CheckReal ((Q * R - A).norm(), (treal) 0., "cvm::srsmatrix QR",  os, __LINE__, dPessimisticSp);
            CheckReal ((eye_real(m) - ~Q * Q).norm(), (treal) 0., "cvm::srsmatrix QR - Q",  os, __LINE__, dPessimisticSp);
        }

        {
            const int m = 10;
            scmatrix A (m);
            A.randomize_real(-10., 10.);
            A.randomize_imag(-10., 10.);

            cvm::scmatrix Q(m);
            cvm::scmatrix R(m);

            A.qr(Q,R);

            CheckReal ((Q * R - A).norm(), (treal) 0., "cvm::scmatrix QR",  os, __LINE__, dPessimisticSp);
            CheckReal ((eye_complex(m) - ~Q * Q).norm(), (treal) 0., "cvm::scmatrix QR - Q",  os, __LINE__, dPessimisticSp);
        }
        {
            const int m = 10;
            scbmatrix A (m,2,3);
            A.randomize_real(-10., 10.);
            A.randomize_imag(-10., 10.);

            cvm::scmatrix Q(m);
            cvm::scmatrix R(m);

            A.qr(Q,R);

            CheckReal ((Q * R - A).norm(), (treal) 0., "cvm::scbmatrix QR",  os, __LINE__, dPessimisticSp);
            CheckReal ((eye_complex(m) - ~Q * Q).norm(), (treal) 0., "cvm::scbmatrix QR - Q",  os, __LINE__, dPessimisticSp);
        }
        {
            const int m = 10;
            schmatrix A (m);
            A.randomize_real(-10., 10.);
            A.randomize_imag(-10., 10.);

            cvm::scmatrix Q(m);
            cvm::scmatrix R(m);

            A.qr(Q,R);

            CheckReal ((Q * R - A).norm(), (treal) 0., "cvm::schmatrix QR",  os, __LINE__, dPessimisticSp);
            CheckReal ((eye_complex(m) - ~Q * Q).norm(), (treal) 0., "cvm::schmatrix QR - Q",  os, __LINE__, dPessimisticSp);
        }


        // Case 1: economy mode, A is (m x n) and Q is (m x n) and R is (n x n)
        {
            int m = 10;
            int n = 5;
            cvm::rmatrix A(m,n);
            A.randomize(-10.0, 10.0);

            cvm::rmatrix Q(m,n);
            cvm::srmatrix R(n);
            A.qr(Q,R);

            CheckReal ((Q * R - A).norm(), (treal) 0., "cvm::rmatrix QR economy",  os, __LINE__, dPessimisticSp);
            CheckReal ((eye_real(n) - ~Q * Q).norm(), (treal) 0., "cvm::rmatrix QR - Q economy",  os, __LINE__, dPessimisticSp);
        }
        // Case 2: full mode, A is (m x n) and Q is (m x m) and R is (m x n)
        {
            int m = 10;
            int n = 5;
            cvm::rmatrix A(m,n);
            A.randomize(-10.0, 10.0);

            cvm::srmatrix Q(m);
            cvm::rmatrix R(m,n);
            A.qr(Q,R);

            CheckReal ((Q * R - A).norm(), (treal) 0., "cvm::rmatrix QR full",  os, __LINE__, dPessimisticSp);
            CheckReal ((eye_real(m) - ~Q * Q).norm(), (treal) 0., "cvm::rmatrix QR - Q full",  os, __LINE__, dPessimisticSp);
        }

        // Case 1: economy mode, A is (m x n) and Q is (m x n) and R is (n x n)
        {
            int m = 10;
            int n = 5;
            cvm::cmatrix A(m,n);
            A.randomize_real(-10.0, 10.0);
            A.randomize_imag(-10.0, 10.0);

            cvm::cmatrix Q(m,n);
            cvm::scmatrix R(n);
            A.qr(Q,R);

            CheckReal ((Q * R - A).norm(), (treal) 0., "cvm::cmatrix QR economy",  os, __LINE__, dPessimisticSp);
            CheckReal ((eye_complex(n) - ~Q * Q).norm(), (treal) 0., "cvm::cmatrix QR - Q economy",  os, __LINE__, dPessimisticSp);
        }
        // Case 2: full mode, A is (m x n) and Q is (m x m) and R is (m x n)
        {
            int m = 10;
            int n = 5;
            cvm::cmatrix A(m,n);
            A.randomize_real(-10.0, 10.0);
            A.randomize_imag(-10.0, 10.0);

            cvm::scmatrix Q(m);
            cvm::cmatrix R(m,n);
            A.qr(Q,R);

            CheckReal ((Q * R - A).norm(), (treal) 0., "cvm::cmatrix QR full",  os, __LINE__, dPessimisticSp);
            CheckReal ((eye_complex(m) - ~Q * Q).norm(), (treal) 0., "cvm::cmatrix QR - Q full",  os, __LINE__, dPessimisticSp);
        }

        // 5.5 left eigenvalues
        {
            scmatrix m(3);
            scmatrix e(3);
            cvector cv(3), cv1(3);
            m(1,1)=tcomplex(0.1,0.01); m(1,2)=tcomplex(0.5,0.05); m(1,3)=tcomplex(0.9,0.09);
            m(2,1)=tcomplex(0.2,0.02); m(2,2)=tcomplex(0.6,0.06); m(2,3)=tcomplex(1.0,0.1);
            m(3,1)=tcomplex(0.3,0.03); m(3,2)=tcomplex(0.7,0.07); m(3,3)=tcomplex(1.0,-1.0);

            cv.eig (m, e);
            cv1 = m.eig (e);
            CheckReal    ((cv - cv1).norm(),   0,  "scmatrix eig", os, __LINE__, dPessimisticSp);
            CheckReal    ((m * e(1) - e(1) * cv1(1)).norm(),   0,  "scmatrix eig", os, __LINE__, dPessimisticSp);
            CheckReal    ((m * e(2) - e(2) * cv1(2)).norm(),   0,  "scmatrix eig", os, __LINE__, dPessimisticSp);
            CheckReal    ((m * e(3) - e(3) * cv1(3)).norm(),   0,  "scmatrix eig", os, __LINE__, dPessimisticSp);

            cv.eig (m, e, false);
            cv1 = m.eig (e , false);
            CheckBool    (cv == cv1,   true,  "scmatrix eig, left", os, __LINE__);

            CheckReal    ((~e(1) * m - ~e(1) * cv1(1)).norm(),   0,  "scmatrix eig, left", os, __LINE__, dPessimisticSp);
            CheckReal    ((~e(2) * m - ~e(2) * cv1(2)).norm(),   0,  "scmatrix eig, left", os, __LINE__, dPessimisticSp);
            CheckReal    ((~e(3) * m - ~e(3) * cv1(3)).norm(),   0,  "scmatrix eig, left", os, __LINE__, dPessimisticSp);
        }
        {
            srmatrix m(3);
            scmatrix e(3);
            cvector cv(3), cv1(3);
            m(1,1)=0.1; m(1,2)=0.5; m(1,3)=0.9;
            m(2,1)=0.2; m(2,2)=0.6; m(2,3)=1.0;
            m(3,1)=0.3; m(3,2)=0.7; m(3,3)=1.0;

            cv.eig (m, e);
            cv1 = m.eig (e);
            CheckReal    ((cv - cv1).norm(),   0,  "scmatrix eig", os, __LINE__, dPessimisticSp);
            CheckReal    ((scmatrix(m) * e(1) - e(1) * cv1(1)).norm(),   0,  "scmatrix eig", os, __LINE__, dPessimisticSp);
            CheckReal    ((scmatrix(m) * e(2) - e(2) * cv1(2)).norm(),   0,  "scmatrix eig", os, __LINE__, dPessimisticSp);
            CheckReal    ((scmatrix(m) * e(3) - e(3) * cv1(3)).norm(),   0,  "scmatrix eig", os, __LINE__, dPessimisticSp);

            cv.eig (m, e, false);
            cv1 = m.eig (e , false);
            CheckBool    (cv == cv1,   true,  "scmatrix eig, left", os, __LINE__);

            CheckReal    ((~e(1) * scmatrix(m) - ~e(1) * cv1(1)).norm(),   0,  "scmatrix eig, left", os, __LINE__, dPessimisticSp);
            CheckReal    ((~e(2) * scmatrix(m) - ~e(2) * cv1(2)).norm(),   0,  "scmatrix eig, left", os, __LINE__, dPessimisticSp);
            CheckReal    ((~e(3) * scmatrix(m) - ~e(3) * cv1(3)).norm(),   0,  "scmatrix eig, left", os, __LINE__, dPessimisticSp);
        }

        // 5.5.1 coverage
        {
            cvector v(5);
            v.set_real((treal) 3.45);
            v.set_imag((treal) -4.17);
            CheckComplex (v(1), tcomplex((treal) 3.45, (treal) -4.17), "cvector set_real set_image",  os, __LINE__);
            CheckComplex (v[5], tcomplex((treal) 3.45, (treal) -4.17), "cvector set_real set_image",  os, __LINE__);

            cmatrix m(4,5);
            m.set_real((treal) 3.45);
            m.set_imag((treal) -4.17);
            CheckComplex (m(1,3), tcomplex((treal) 3.45, (treal) -4.17), "cmatrix set_real set_image",  os, __LINE__);
            CheckComplex (m[4][5], tcomplex((treal) 3.45, (treal) -4.17), "cmatrix set_real set_image",  os, __LINE__);

            scmatrix sm(5);
            sm.set_real((treal) 3.45);
            sm.set_imag((treal) -4.17);
            CheckComplex (sm(1,3), tcomplex((treal) 3.45, (treal) -4.17), "scmatrix set_real set_image",  os, __LINE__);
            CheckComplex (sm[5][5], tcomplex((treal) 3.45, (treal) -4.17), "scmatrix set_real set_image",  os, __LINE__);

            scbmatrix bm(5,1,2);
            bm.set_real((treal) 3.45);
            bm.set_imag((treal) -4.17);
            CheckComplex (bm(1,3), tcomplex((treal) 3.45, (treal) -4.17), "scbmatrix set_real set_image",  os, __LINE__);
            CheckComplex (bm[5][5], tcomplex((treal) 3.45, (treal) -4.17), "scbmatrix set_real set_image",  os, __LINE__);

            schmatrix hm(5);
            hm.set_real((treal) 3.45);
            CheckComplex (hm(1,3), tcomplex((treal) 3.45, (treal) 0.), "schmatrix set_real",  os, __LINE__);

            bool ex = false;
            try {
                hm.set_imag((treal) -4.17);
            }
            catch (cvm::cvmexception&)
            {
                ex = true;
            }
            CheckBool (ex, true, "schmatrix set_image not allowed",  os, __LINE__);
        }

        {
            schmatrix hm(5);
            hm.randomize_real((treal)3., (treal)7.);
            hm.randomize_imag((treal)-5., (treal)4.);
            cvector x(5), b(5);
            b.randomize_real((treal)-4., (treal)9.);
            b.randomize_imag((treal)-2., (treal)1.);

            x = hm.solve (b);
            CheckReal    ((hm * x - b).norm(),   0,  "schmatrix solve", os, __LINE__, dPessimisticSp);
            treal err;
            x = hm.solve (b, err);
            CheckReal    ((hm * x - b).norm(),   0,  "schmatrix solve", os, __LINE__, dPessimisticSp);
            CheckReal    (err,   0,  "schmatrix solve", os, __LINE__, dVeryPessimisticSp);

            cmatrix mb(5,6), mx(5,6);
            mb.randomize_real((treal)-4., (treal)9.);
            mb.randomize_imag((treal)-2., (treal)1.);
            
            mx = hm.solve (mb);
            CheckReal    ((hm * mx - mb).norm(),   0,  "schmatrix solve", os, __LINE__, dPessimisticSp);

            mx = hm.solve (mb, err);
            CheckReal    ((hm * mx - mb).norm(),   0,  "schmatrix solve", os, __LINE__, dPessimisticSp);
            CheckReal    (err,   0,  "schmatrix solve", os, __LINE__, dVeryPessimisticSp);

            schmatrix im(hm.inv());
            CheckReal    ((im * hm - eye_complex(5)).norm(),   0,  "srsmatrix inv", os, __LINE__, dPessimisticSp);
            im.inv(hm);
            CheckReal    ((im * hm - eye_complex(5)).norm(),   0,  "srsmatrix inv", os, __LINE__, dPessimisticSp);

            rvector ev(5), ev1(5), ev2(5);
            scmatrix evect(5);
            ev.eig (hm, evect);
            ev1 = hm.eig(evect);
            ev2 = hm.eig();

            CheckReal    ((ev - ev1).norm(),   0,  "srsmatrix eig", os, __LINE__, dPessimisticSp);
            CheckReal    ((ev - ev2).norm(),   0,  "srsmatrix eig", os, __LINE__, dPessimisticSp);
            CheckReal    ((hm * cvector(evect(1)) - evect(1) * ev1(1)).norm(),   0,  "srsmatrix eig", os, __LINE__, dPessimisticSp);
            CheckReal    ((hm * cvector(evect(2)) - evect(2) * ev1(2)).norm(),   0,  "srsmatrix eig", os, __LINE__, dPessimisticSp);
            CheckReal    ((hm * cvector(evect(3)) - evect(3) * ev1(3)).norm(),   0,  "srsmatrix eig", os, __LINE__, dPessimisticSp);
            CheckReal    ((hm * cvector(evect(4)) - evect(4) * ev1(4)).norm(),   0,  "srsmatrix eig", os, __LINE__, dPessimisticSp);
            CheckReal    ((hm * cvector(evect(5)) - evect(5) * ev1(5)).norm(),   0,  "srsmatrix eig", os, __LINE__, dPessimisticSp);
        }

        {
            srsmatrix m(5);
            m.randomize((treal)3., (treal)7.);
            rvector x(5), b(5);
            b.randomize((treal)-4., (treal)9.);

            x = m.solve (b);
            CheckReal    ((m * x - b).norm(),   0,  "srsmatrix solve", os, __LINE__, dPessimisticSp);
            treal err;
            x = m.solve (b, err);
            CheckReal    ((m * x - b).norm(),   0,  "srsmatrix solve", os, __LINE__, dPessimisticSp);
            CheckReal    (err,   0,  "srsmatrix solve", os, __LINE__, dVeryPessimisticSp);

            rmatrix mb(5,6), mx(5,6);
            mb.randomize((treal)-4., (treal)9.);
            
            mx = m.solve (mb);
            CheckReal    ((m * mx - mb).norm(),   0,  "srsmatrix solve", os, __LINE__, dPessimisticSp);

            mx = m.solve (mb, err);
            CheckReal    ((m * mx - mb).norm(),   0,  "srsmatrix solve", os, __LINE__, dPessimisticSp);
            CheckReal    (err,   0,  "srsmatrix solve", os, __LINE__, dVeryPessimisticSp);

            srsmatrix im(m.inv());
            CheckReal    ((im * m - eye_real(5)).norm(),   0,  "srsmatrix inv", os, __LINE__, dPessimisticSp);
            im.inv(m);
            CheckReal    ((im * m - eye_real(5)).norm(),   0,  "srsmatrix inv", os, __LINE__, dPessimisticSp);

            rvector ev(5), ev1(5), ev2(5);
            srmatrix evect(5);
            ev.eig (m, evect);
            ev1 = m.eig(evect);
            ev2 = m.eig();

            CheckReal    ((ev - ev1).norm(),   0,  "srsmatrix eig", os, __LINE__, dPessimisticSp);
            CheckReal    ((ev - ev2).norm(),   0,  "srsmatrix eig", os, __LINE__, dPessimisticSp);
            CheckReal    ((m * evect(1) - evect(1) * ev1(1)).norm(),   0,  "srsmatrix eig", os, __LINE__, dPessimisticSp);
            CheckReal    ((m * evect(2) - evect(2) * ev1(2)).norm(),   0,  "srsmatrix eig", os, __LINE__, dPessimisticSp);
            CheckReal    ((m * evect(3) - evect(3) * ev1(3)).norm(),   0,  "srsmatrix eig", os, __LINE__, dPessimisticSp);
            CheckReal    ((m * evect(4) - evect(4) * ev1(4)).norm(),   0,  "srsmatrix eig", os, __LINE__, dPessimisticSp);
            CheckReal    ((m * evect(5) - evect(5) * ev1(5)).norm(),   0,  "srsmatrix eig", os, __LINE__, dPessimisticSp);

        }


        // positive-definite:
        {
/*
mpositive = zeros(3,3);
m(1,1) = 1;       m(1,2) = 2+i;     m(1,3) = 1-2*i;
m(2,1) = 2-i;     m(2,2) = 15;      m(2,3) = -1-3*i;
m(3,1) = 1+2*i;   m(3,2) = -1+3*i;  m(3,3) = 20;
det (m)
*/
            treal r[] = {1., 2., 1., 2., 15., -1., 1., -1., 20.};
            treal i[] = {0., -1., 2., 1., 0., 3., -2., -3., 0.};
            const schmatrix m(r, i, 3);
            scmatrix c(3);
            c.cholesky(m);
            CheckReal    ((~c * c - m).norm(),   0,  "schmatrix cholesky", os, __LINE__, dPessimisticSp);
            CheckComplex (m.det(), tcomplex((treal) 145., (treal) 0.), "schmatrix det (positive_defiite)",  os, __LINE__, dPessimisticSp);
        }
        // not positive-definite:
        {
            treal r[] = {1., 2., 1., 2., 5., -1., 1., -1., 20.};
            treal i[] = {0., -1., 2., 1., 0., 3., -2., -3., 0.};
            const schmatrix m(r, i, 3);
            CheckComplex (m.det(), tcomplex((treal) -5., (treal) 0.), "schmatrix det (not positive_defiite)",  os, __LINE__, dPessimisticSp);
        }



        {
            LockIt l;
            os << "***END***" << std::endl;
            std::cout << "ALL TESTS SUCEEDED" << std::endl;
        }
    }
    catch (std::exception& e) {
        LockIt l;
        std::cout << "Exception: " << e.what () << std::endl;
    }

//    cvmExit ();

    return NULL;
}



int main(int argc, char* argv[])
{
    os << "rvector   " << sizeof (rvector) << std::endl;
    os << "cvector   " << sizeof (cvector) << std::endl;
    os << "rmatrix   " << sizeof (rmatrix) << std::endl;
    os << "cmatrix   " << sizeof (cmatrix) << std::endl;
    os << "srmatrix  " << sizeof (srmatrix) << std::endl;
    os << "scmatrix  " << sizeof (scmatrix) << std::endl;
    os << "srbmatrix " << sizeof (srbmatrix) << std::endl;
    os << "scbmatrix " << sizeof (scbmatrix) << std::endl;
    os << "srsmatrix " << sizeof (srsmatrix) << std::endl;
    os << "schmatrix " << sizeof (schmatrix) << std::endl;

    int i;
    int nThreads = 1;
    int nRuns = 1;

    for (i = 1; i < argc; ++i)
    {
        if (argv[i] != NULL && strlen(argv[i]) > 2)
        {
            char cSep = argv[i][0];
            char cKey = argv[i][1];
            if (strchr ("-/", cSep))
            {
                if (cKey == 't' || cKey == 'T')
                {
                    nThreads = atol(argv[i] + 2);
                    if (nThreads <= 0) nThreads = 1;
                }
                if (cKey == 'r' || cKey == 'R')
                {
                    nRuns = atol(argv[i] + 2);
                    if (nRuns <= 0) nRuns = 1;
                }
            }
        }
    }

#if !defined (WIN32) && !defined (_WIN32)
    if (nThreads > 2)
    {
        std::cout << "Warning! Running 2 threads only!" << std::endl;
    }
#else
    unsigned int unThreadId = 1;
#endif

    clock_t start, finish;
    start = clock();

    for (i = 0; i < nRuns; ++i)
    {
#if defined (WIN32) || defined (_WIN32)

        if (nThreads > 1)
        {
            int i;
            std::vector<HANDLE> tha(nThreads);

            for (i = 0; i < nThreads; i++)
            {
                tha[i] = (HANDLE) _beginthreadex (NULL, 0, TestBody, NULL, 0, &unThreadId);
                unThreadId++;
            }

            WaitForMultipleObjects(nThreads, &tha[0], TRUE, INFINITE);

            for (i = nThreads - 1; i >= 0; i--)
            {
                ::CloseHandle(tha[i]);
            }
        }
        else
        {
            TestBody (NULL);
        }

#else

        pthread_t tid = 0;
        void* tr = NULL;

        if (nThreads > 1)
        {
            LockIt l;
            if (pthread_create (&tid, NULL, TestBody, NULL) != 0)
            {
                std::cout << "Failed to create thread!" << std::endl;
            }
            else
            {
                std::cout << tid << std::endl;
            }
        }

        TestBody (NULL);

        if (nThreads > 1)
        {
            if (pthread_join (tid, &tr) != 0)
            {
                std::cout << "Failed to join thread!" << std::endl;
            }
        }
#endif
    }

    finish = clock();
    double duration = (double) (finish - start) / CLOCKS_PER_SEC;    
    std::cout << "TOTAL TIME " << duration << " sec." << std::endl;
    return 1;
}
