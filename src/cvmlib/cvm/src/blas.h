/* CVM Class Library */
/* Copyright (C), Sergei Nikolaev, 1992-2007, http://cvmlib.com */
/* BLAS and LAPACK functions declaration */
/* Also contains my fortran utilities */

#ifndef _BLAS_H
#define _BLAS_H

#ifdef __cplusplus
extern "C" {
#endif

#if defined (_MSC_VER)
#    if defined (CVM_ACML)
#        define CVM_FTN_CALL
#        define CVM_STD_CALL
#    else
#        define CVM_FTN_CALL __stdcall
#        define CVM_STD_CALL __stdcall
#    endif
#else
#    if defined (__BORLANDC__)

#        define CVM_FTN_CALL __stdcall
#        define CVM_STD_CALL

// blas & lapack stuff

#        define ISAMAX  isamax
#        define IDAMAX  idamax
#        define ISAMIN  isamin
#        define IDAMIN  idamin

#        define ICAMAX  icamax
#        define IZAMAX  izamax
#        define ICAMIN  icamin
#        define IZAMIN  izamin

#        define SNRM2   snrm2
#        define DNRM2   dnrm2
#        define SCNRM2  scnrm2
#        define DZNRM2  dznrm2

#        define SSWAP   sswap
#        define DSWAP   dswap
#        define CSWAP   cswap
#        define ZSWAP   zswap

#        define SDOT    sdot
#        define DDOT    ddot

// complex dot wrappers:

#    if defined (CVM_COMPLEX_NUMBER_RETURNED)
#        define VCDOTU  cdotu
#        define VZDOTU  zdotu
#        define VCDOTC  cdotc
#        define VZDOTC  zdotc
#    else
#        define VCDOTU  vcdotu
#        define VZDOTU  vzdotu
#        define VCDOTC  vcdotc
#        define VZDOTC  vzdotc
#    endif

#        define SAXPY   saxpy
#        define DAXPY   daxpy
#        define CAXPY   caxpy
#        define ZAXPY   zaxpy

#        define SCOPY   scopy
#        define DCOPY   dcopy
#        define CCOPY   ccopy
#        define ZCOPY   zcopy

#        define SSCAL   sscal
#        define DSCAL   dscal
#        define CSCAL   cscal
#        define ZSCAL   zscal

#        define CSSCAL  csscal
#        define ZDSCAL  zdscal

#        define SGER    sger
#        define DGER    dger

#        define CGERU   cgeru
#        define ZGERU   zgeru

#        define CGERC   cgerc
#        define ZGERC   zgerc

#        define SGEMV   sgemv
#        define DGEMV   dgemv
#        define CGEMV   cgemv
#        define ZGEMV   zgemv

#        define SGBMV   sgbmv
#        define DGBMV   dgbmv
#        define CGBMV   cgbmv
#        define ZGBMV   zgbmv

#        define SGEMM   sgemm
#        define DGEMM   dgemm
#        define CGEMM   cgemm
#        define ZGEMM   zgemm

#        define SGETRF  sgetrf
#        define DGETRF  dgetrf
#        define CGETRF  cgetrf
#        define ZGETRF  zgetrf

#        define SGBTRF  sgbtrf
#        define DGBTRF  dgbtrf
#        define CGBTRF  cgbtrf
#        define ZGBTRF  zgbtrf

#        define SGETRS  sgetrs
#        define DGETRS  dgetrs
#        define CGETRS  cgetrs
#        define ZGETRS  zgetrs

#        define SGBTRS  sgbtrs
#        define DGBTRS  dgbtrs
#        define CGBTRS  cgbtrs
#        define ZGBTRS  zgbtrs

#        define SGERFS  sgerfs
#        define DGERFS  dgerfs
#        define CGERFS  cgerfs
#        define ZGERFS  zgerfs

#        define SGBRFS  sgbrfs
#        define DGBRFS  dgbrfs
#        define CGBRFS  cgbrfs
#        define ZGBRFS  zgbrfs

#        define SGETRI  sgetri
#        define DGETRI  dgetri
#        define CGETRI  cgetri
#        define ZGETRI  zgetri

#        define SGEBRD  sgebrd
#        define DGEBRD  dgebrd
#        define CGEBRD  cgebrd
#        define ZGEBRD  zgebrd

#        define SGBBRD  sgbbrd
#        define DGBBRD  dgbbrd
#        define CGBBRD  cgbbrd
#        define ZGBBRD  zgbbrd

#        define SORGBR  sorgbr
#        define DORGBR  dorgbr

#        define CUNGBR  cungbr
#        define ZUNGBR  zungbr

#        define SBDSQR  sbdsqr
#        define DBDSQR  dbdsqr
#        define CBDSQR  cbdsqr
#        define ZBDSQR  zbdsqr

#        define SGEBAL  sgebal
#        define DGEBAL  dgebal
#        define CGEBAL  cgebal
#        define ZGEBAL  zgebal

#        define SGEHRD  sgehrd
#        define DGEHRD  dgehrd
#        define CGEHRD  cgehrd
#        define ZGEHRD  zgehrd

#        define SORGHR  sorghr
#        define DORGHR  dorghr

#        define CUNGHR  cunghr
#        define ZUNGHR  zunghr

#        define SHSEQR  shseqr
#        define DHSEQR  dhseqr
#        define CHSEQR  chseqr
#        define ZHSEQR  zhseqr

#        define STREVC  strevc
#        define DTREVC  dtrevc
#        define CTREVC  ctrevc
#        define ZTREVC  ztrevc

#        define SGEBAK  sgebak
#        define DGEBAK  dgebak
#        define CGEBAK  cgebak
#        define ZGEBAK  zgebak

#        define SGECON  sgecon
#        define DGECON  dgecon
#        define CGECON  cgecon
#        define ZGECON  zgecon

#        define SSPMV   sspmv
#        define DSPMV   dspmv

#        define SSYMM   ssymm
#        define DSYMM   dsymm
#        define CSYMM   csymm
#        define ZSYMM   zsymm
#        define CHEMM   chemm
#        define ZHEMM   zhemm

#        define SPOTRF  spotrf
#        define DPOTRF  dpotrf
#        define CPOTRF  cpotrf
#        define ZPOTRF  zpotrf

#        define SSYTRF  ssytrf
#        define DSYTRF  dsytrf
#        define CSYTRF  csytrf
#        define ZSYTRF  zsytrf
#        define CHETRF  chetrf
#        define ZHETRF  zhetrf

#        define SPOTRS  spotrs
#        define DPOTRS  dpotrs
#        define CPOTRS  cpotrs
#        define ZPOTRS  zpotrs

#        define SPORFS  sporfs
#        define DPORFS  dporfs
#        define CPORFS  cporfs
#        define ZPORFS  zporfs

#        define SSYTRS  ssytrs
#        define DSYTRS  dsytrs
#        define CSYTRS  csytrs
#        define ZSYTRS  zsytrs
#        define CHETRS  chetrs
#        define ZHETRS  zhetrs

#        define SSYRFS  ssyrfs
#        define DSYRFS  dsyrfs
#        define CSYRFS  csyrfs
#        define ZSYRFS  zsyrfs
#        define CHERFS  cherfs
#        define ZHERFS  zherfs

#        define SPOTRI  spotri
#        define DPOTRI  dpotri
#        define CPOTRI  cpotri
#        define ZPOTRI  zpotri

#        define SSYTRI  ssytri
#        define DSYTRI  dsytri
#        define CSYTRI  csytri
#        define ZSYTRI  zsytri
#        define CHETRI  chetri
#        define ZHETRI  zhetri

#        define SSYEVD  ssyevd
#        define DSYEVD  dsyevd
#        define CHEEVD  cheevd
#        define ZHEEVD  zheevd

#        define SPOEQU  spoequ
#        define DPOEQU  dpoequ
#        define CPOEQU  cpoequ
#        define ZPOEQU  zpoequ

#        define SSYMV   ssymv
#        define DSYMV   dsymv
#        define CHEMV   chemv
#        define ZHEMV   zhemv

#        define SSYRK   ssyrk
#        define DSYRK   dsyrk
#        define CSYRK   csyrk
#        define ZSYRK   zsyrk
#        define CHERK   cherk
#        define ZHERK   zherk

#        define SSYR2K  ssyr2k
#        define DSYR2K  dsyr2k
#        define CSYR2K  csyr2k
#        define ZSYR2K  zsyr2k
#        define CHER2K  cher2k
#        define ZHER2K  zher2k

#        define SGEQRF  sgeqrf
#        define DGEQRF  dgeqrf
#        define CGEQRF  cgeqrf
#        define ZGEQRF  zgeqrf

#        define SORGQR  sorgqr
#        define DORGQR  dorgqr
#        define CUNGQR  cungqr
#        define ZUNGQR  zungqr

#    else      // !__BORLANDC__

#        define  CVM_FTN_CALL
#        define  CVM_STD_CALL

// my fortran stuff

#        define DPOLY   dpoly_
#        define SPOLY   spoly_
#        define CPOLY   cpoly_
#        define ZPOLY   zpoly_

#        define NPOLY   npoly_

#        define SMEXP   smexp_
#        define DMEXP   dmexp_
#        define CMEXP   cmexp_
#        define ZMEXP   zmexp_

#        define SMEXPC  smexpc_
#        define DMEXPC  dmexpc_
#        define CMEXPC  cmexpc_
#        define ZMEXPC  zmexpc_

// blas & lapack stuff

#        define ISAMAX  isamax_
#        define IDAMAX  idamax_
#        define ISAMIN  isamin_
#        define IDAMIN  idamin_

#        define ICAMAX  icamax_
#        define IZAMAX  izamax_
#        define ICAMIN  icamin_
#        define IZAMIN  izamin_

#        define SNRM2   snrm2_
#        define DNRM2   dnrm2_
#        define SCNRM2  scnrm2_
#        define DZNRM2  dznrm2_

#        define SSWAP   sswap_
#        define DSWAP   dswap_
#        define CSWAP   cswap_
#        define ZSWAP   zswap_

#        define SDOT    sdot_
#        define DDOT    ddot_

// complex dot wrappers

#    if defined (CVM_COMPLEX_NUMBER_RETURNED)
#        define VCDOTU  cdotu_
#        define VZDOTU  zdotu_
#        define VCDOTC  cdotc_
#        define VZDOTC  zdotc_
#    else
#        define VCDOTU  vcdotu_
#        define VZDOTU  vzdotu_
#        define VCDOTC  vcdotc_
#        define VZDOTC  vzdotc_
#    endif

#        define SAXPY   saxpy_
#        define DAXPY   daxpy_
#        define CAXPY   caxpy_
#        define ZAXPY   zaxpy_

#        define SCOPY   scopy_
#        define DCOPY   dcopy_
#        define CCOPY   ccopy_
#        define ZCOPY   zcopy_

#        define SSCAL   sscal_
#        define DSCAL   dscal_
#        define CSCAL   cscal_
#        define ZSCAL   zscal_

#        define CSSCAL  csscal_
#        define ZDSCAL  zdscal_

#        define SGER    sger_
#        define DGER    dger_

#        define CGERU   cgeru_
#        define ZGERU   zgeru_

#        define CGERC   cgerc_
#        define ZGERC   zgerc_

#        define SGEMV   sgemv_
#        define DGEMV   dgemv_
#        define CGEMV   cgemv_
#        define ZGEMV   zgemv_

#        define SGBMV   sgbmv_
#        define DGBMV   dgbmv_
#        define CGBMV   cgbmv_
#        define ZGBMV   zgbmv_

#        define SGEMM   sgemm_
#        define DGEMM   dgemm_
#        define CGEMM   cgemm_
#        define ZGEMM   zgemm_

#        define SGETRF  sgetrf_
#        define DGETRF  dgetrf_
#        define CGETRF  cgetrf_
#        define ZGETRF  zgetrf_

#        define SGBTRF  sgbtrf_
#        define DGBTRF  dgbtrf_
#        define CGBTRF  cgbtrf_
#        define ZGBTRF  zgbtrf_

#        define SGETRS  sgetrs_
#        define DGETRS  dgetrs_
#        define CGETRS  cgetrs_
#        define ZGETRS  zgetrs_

#        define SGBTRS  sgbtrs_
#        define DGBTRS  dgbtrs_
#        define CGBTRS  cgbtrs_
#        define ZGBTRS  zgbtrs_

#        define SGERFS  sgerfs_
#        define DGERFS  dgerfs_
#        define CGERFS  cgerfs_
#        define ZGERFS  zgerfs_

#        define SGBRFS  sgbrfs_
#        define DGBRFS  dgbrfs_
#        define CGBRFS  cgbrfs_
#        define ZGBRFS  zgbrfs_

#        define SGETRI  sgetri_
#        define DGETRI  dgetri_
#        define CGETRI  cgetri_
#        define ZGETRI  zgetri_

#        define SGEBRD  sgebrd_
#        define DGEBRD  dgebrd_
#        define CGEBRD  cgebrd_
#        define ZGEBRD  zgebrd_

#        define SGBBRD  sgbbrd_
#        define DGBBRD  dgbbrd_
#        define CGBBRD  cgbbrd_
#        define ZGBBRD  zgbbrd_

#        define SORGBR  sorgbr_
#        define DORGBR  dorgbr_

#        define CUNGBR  cungbr_
#        define ZUNGBR  zungbr_

#        define SBDSQR  sbdsqr_
#        define DBDSQR  dbdsqr_
#        define CBDSQR  cbdsqr_
#        define ZBDSQR  zbdsqr_

#        define SGEBAL  sgebal_
#        define DGEBAL  dgebal_
#        define CGEBAL  cgebal_
#        define ZGEBAL  zgebal_

#        define SGEHRD  sgehrd_
#        define DGEHRD  dgehrd_
#        define CGEHRD  cgehrd_
#        define ZGEHRD  zgehrd_

#        define SORGHR  sorghr_
#        define DORGHR  dorghr_

#        define CUNGHR  cunghr_
#        define ZUNGHR  zunghr_

#        define SHSEQR  shseqr_
#        define DHSEQR  dhseqr_
#        define CHSEQR  chseqr_
#        define ZHSEQR  zhseqr_

#        define STREVC  strevc_
#        define DTREVC  dtrevc_
#        define CTREVC  ctrevc_
#        define ZTREVC  ztrevc_

#        define SGEBAK  sgebak_
#        define DGEBAK  dgebak_
#        define CGEBAK  cgebak_
#        define ZGEBAK  zgebak_

#        define SGECON  sgecon_
#        define DGECON  dgecon_
#        define CGECON  cgecon_
#        define ZGECON  zgecon_

#        define SSPMV   sspmv_
#        define DSPMV   dspmv_

#        define SSYMM   ssymm_
#        define DSYMM   dsymm_
#        define CSYMM   csymm_
#        define ZSYMM   zsymm_
#        define CHEMM   chemm_
#        define ZHEMM   zhemm_

#        define SPOTRF  spotrf_
#        define DPOTRF  dpotrf_
#        define CPOTRF  cpotrf_
#        define ZPOTRF  zpotrf_

#        define SSYTRF  ssytrf_
#        define DSYTRF  dsytrf_
#        define CSYTRF  csytrf_
#        define ZSYTRF  zsytrf_
#        define CHETRF  chetrf_
#        define ZHETRF  zhetrf_

#        define SPOTRS  spotrs_
#        define DPOTRS  dpotrs_
#        define CPOTRS  cpotrs_
#        define ZPOTRS  zpotrs_

#        define SPORFS  sporfs_
#        define DPORFS  dporfs_
#        define CPORFS  cporfs_
#        define ZPORFS  zporfs_

#        define SSYTRS  ssytrs_
#        define DSYTRS  dsytrs_
#        define CSYTRS  csytrs_
#        define ZSYTRS  zsytrs_
#        define CHETRS  chetrs_
#        define ZHETRS  zhetrs_

#        define SSYRFS  ssyrfs_
#        define DSYRFS  dsyrfs_
#        define CSYRFS  csyrfs_
#        define ZSYRFS  zsyrfs_
#        define CHERFS  cherfs_
#        define ZHERFS  zherfs_

#        define SPOTRI  spotri_
#        define DPOTRI  dpotri_
#        define CPOTRI  cpotri_
#        define ZPOTRI  zpotri_

#        define SSYTRI  ssytri_
#        define DSYTRI  dsytri_
#        define CSYTRI  csytri_
#        define ZSYTRI  zsytri_
#        define CHETRI  chetri_
#        define ZHETRI  zhetri_

#        define SSYEVD  ssyevd_
#        define DSYEVD  dsyevd_
#        define CHEEVD  cheevd_
#        define ZHEEVD  zheevd_

#        define SPOEQU  spoequ_
#        define DPOEQU  dpoequ_
#        define CPOEQU  cpoequ_
#        define ZPOEQU  zpoequ_

#        define SSYMV   ssymv_
#        define DSYMV   dsymv_
#        define CHEMV   chemv_
#        define ZHEMV   zhemv_

#        define SSYRK   ssyrk_
#        define DSYRK   dsyrk_
#        define CSYRK   csyrk_
#        define ZSYRK   zsyrk_
#        define CHERK   cherk_
#        define ZHERK   zherk_

#        define SSYR2K  ssyr2k_
#        define DSYR2K  dsyr2k_
#        define CSYR2K  csyr2k_
#        define ZSYR2K  zsyr2k_
#        define CHER2K  cher2k_
#        define ZHER2K  zher2k_

#        define SGEQRF  sgeqrf_
#        define DGEQRF  dgeqrf_
#        define CGEQRF  cgeqrf_
#        define ZGEQRF  zgeqrf_

#        define SORGQR  sorgqr_
#        define DORGQR  dorgqr_
#        define CUNGQR  cungqr_
#        define ZUNGQR  zungqr_

#    endif      // !__BORLANDC__
#endif      // !_MSC_VER


void  CVM_FTN_CALL SPOLY       (const int* m,
                                const float* a,
                                const int* lda,
                                const int* n,
                                const float* v,
                                      float* p,
                                const int* ldp,
                                      float* r);

void  CVM_FTN_CALL DPOLY       (const int* m,
                                const double* a,
                                const int* lda,
                                const int* n,
                                const double* v,
                                      double* p,
                                const int* ldp,
                                      double* r);

void  CVM_FTN_CALL CPOLY       (const int* m,
                                const std::complex<float>* a,
                                const int* lda,
                                const int* n,
                                const std::complex<float>* v,
                                      std::complex<float>* p,
                                const int* ldp,
                                      std::complex<float>* r);

void  CVM_FTN_CALL ZPOLY       (const int* m,
                                const std::complex<double>* a,
                                const int* lda,
                                const int* n,
                                const std::complex<double>* v,
                                      std::complex<double>* p,
                                const int* ldp,
                                      std::complex<double>* r);

int   CVM_FTN_CALL NPOLY       (const int* m,
                                const int* n);

int   CVM_STD_CALL ISAMAX      (const int* n,
                                const float* x, 
                                const int* incx);
int   CVM_STD_CALL IDAMAX      (const int* n,
                                const double* x, 
                                const int* incx);

int   CVM_STD_CALL ISAMIN      (const int* n,
                                const float* x,
                                const int* incx);
int   CVM_STD_CALL IDAMIN      (const int* n,
                                const double* x,
                                const int* incx);

int   CVM_STD_CALL ICAMAX      (const int* n,
                                const std::complex<float>* x,
                                const int* incx);
int   CVM_STD_CALL IZAMAX      (const int* n,
                                const std::complex<double>* x,
                                const int* incx);

int   CVM_STD_CALL ICAMIN      (const int* n,
                                const std::complex<float>* x,
                                const int* incx);
int   CVM_STD_CALL IZAMIN      (const int* n,
                                const std::complex<double>* x,
                                const int* incx);

float  CVM_STD_CALL SNRM2      (const int* n,
                                const float* x,
                                const int* incx);
double CVM_STD_CALL DNRM2      (const int* n,
                                const double* x,
                                const int* incx);

float  CVM_STD_CALL SCNRM2     (const int* n,
                                const std::complex<float>* x,
                                const int* incx);
double CVM_STD_CALL DZNRM2     (const int* n,
                                const std::complex<double>* x,
                                const int* incx);

void  CVM_STD_CALL SSWAP       (const int* n,
                                      float* x,
                                const int* incx,
                                      float* y,
                                const int* incy);
void  CVM_STD_CALL DSWAP       (const int* n,
                                      double* x,
                                const int* incx,
                                      double* y,
                                const int* incy);

void  CVM_STD_CALL CSWAP       (const int* n, 
                                      std::complex<float>* x, 
                                const int* incx, 
                                      std::complex<float>* y, 
                                const int* incy);
void  CVM_STD_CALL ZSWAP       (const int* n, 
                                      std::complex<double>* x, 
                                const int* incx, 
                                      std::complex<double>* y, 
                                const int* incy);

float  CVM_STD_CALL SDOT       (const int* n,
                                const float* x, 
                                const int* incx, 
                                const float* y, 
                                const int* incy); 
double CVM_STD_CALL DDOT       (const int* n,
                                const double* x, 
                                const int* incx, 
                                const double* y, 
                                const int* incy); 

// complex dot wrappers
void  CVM_STD_CALL VCDOTU      (std::complex<float>* dot,
                                const int* n,
                                const std::complex<float>* x, 
                                const int* incx, 
                                const std::complex<float>* y,
                                const int* incy);
void  CVM_STD_CALL VZDOTU      (std::complex<double>* dot,
                                const int* n,
                                const std::complex<double>* x, 
                                const int* incx, 
                                const std::complex<double>* y,
                                const int* incy);
void  CVM_STD_CALL VCDOTC      (std::complex<float>* dot,
                                const int* n,
                                const std::complex<float>* x, 
                                const int* incx, 
                                const std::complex<float>* y, 
                                const int* incy);
void  CVM_STD_CALL VZDOTC      (std::complex<double>* dot,
                                const int* n,
                                const std::complex<double>* x, 
                                const int* incx, 
                                const std::complex<double>* y, 
                                const int* incy);

void  CVM_STD_CALL SAXPY       (const int* n,
                                const float* a, 
                                const float* x, 
                                const int* incx, 
                                      float* y, 
                                const int* incy); 
void  CVM_STD_CALL DAXPY       (const int* n,
                                const double* a, 
                                const double* x, 
                                const int* incx, 
                                      double* y, 
                                const int* incy); 

void  CVM_STD_CALL CAXPY       (const int* n, 
                                const std::complex<float>* a,
                                const std::complex<float>* x, 
                                const int* incx, 
                                      std::complex<float>* y,
                                const int* incy); 
void  CVM_STD_CALL ZAXPY       (const int* n, 
                                const std::complex<double>* a,
                                const std::complex<double>* x, 
                                const int* incx, 
                                      std::complex<double>* y,
                                const int* incy); 

void  CVM_STD_CALL DCOPY       (const int* n,
                                const double* x,
                                const int* incx, 
                                      double* y, 
                                const int* incy); 
void  CVM_STD_CALL SCOPY       (const int* n,
                                const float* x,
                                const int* incx, 
                                      float* y, 
                                const int* incy); 

void  CVM_STD_CALL CCOPY       (const int* n,
                                const std::complex<float>* x, 
                                const int* incx, 
                                      std::complex<float>* y, 
                                const int* incy); 
void  CVM_STD_CALL ZCOPY       (const int* n,
                                const std::complex<double>* x, 
                                const int* incx, 
                                      std::complex<double>* y, 
                                const int* incy); 

void  CVM_STD_CALL SSCAL       (const int* n,
                                const float* a, 
                                      float* x, 
                                const int* incx);
void  CVM_STD_CALL DSCAL       (const int* n,
                                const double* a, 
                                      double* x, 
                                const int* incx);

void  CVM_STD_CALL CSCAL       (const int* n,
                                const std::complex<float>* a,
                                      std::complex<float>* x, 
                                const int* incx); 
void  CVM_STD_CALL ZSCAL       (const int* n,
                                const std::complex<double>* a,
                                      std::complex<double>* x, 
                                const int* incx); 

void  CVM_STD_CALL CSSCAL      (const int* n,
                                const float* a, 
                                      std::complex<float>* x, 
                                const int* incx);
void  CVM_STD_CALL ZDSCAL      (const int* n,
                                const double* a, 
                                      std::complex<double>* x, 
                                const int* incx);

void  CVM_STD_CALL SGER        (const int* m,
                                const int* n,
                                const float* alpha,
                                const float* x,
                                const int* incx,
                                const float* y,
                                const int* incy,
                                      float* a,
                                const int* lda);
void  CVM_STD_CALL DGER        (const int* m,
                                const int* n,
                                const double* alpha,
                                const double* x,
                                const int* incx,
                                const double* y,
                                const int* incy,
                                      double* a,
                                const int* lda);

void  CVM_STD_CALL CGERU       (const int* m,
                                const int* n,
                                const std::complex<float>* alpha,
                                const std::complex<float>* x,
                                const int* incx,
                                const std::complex<float>* y,
                                const int* incy,
                                      std::complex<float>* a,
                                const int* lda);
void  CVM_STD_CALL ZGERU       (const int* m,
                                const int* n,
                                const std::complex<double>* alpha,
                                const std::complex<double>* x,
                                const int* incx,
                                const std::complex<double>* y,
                                const int* incy,
                                      std::complex<double>* a,
                                const int* lda);

void  CVM_STD_CALL CGERC       (const int* m,
                                const int* n,
                                const std::complex<float>* alpha,
                                const std::complex<float>* x,
                                const int* incx,
                                const std::complex<float>* y,
                                const int* incy,
                                      std::complex<float>* a,
                                const int* lda);
void  CVM_STD_CALL ZGERC       (const int* m,
                                const int* n,
                                const std::complex<double>* alpha,
                                const std::complex<double>* x,
                                const int* incx,
                                const std::complex<double>* y,
                                const int* incy,
                                      std::complex<double>* a,
                                const int* lda);

void  CVM_STD_CALL SGEMV       (const char* trans,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int transsz,
#endif
                                const int* m,
                                const int* n,
                                const float* alpha,
                                const float* a,
                                const int* lda,
                                const float* x,
                                const int* incx,
                                const float* beta,
                                      float* y,
                                const int* incy);
void  CVM_STD_CALL DGEMV       (const char* trans,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int transsz,
#endif
                                const int* m,
                                const int* n,
                                const double* alpha,
                                const double* a,
                                const int* lda,
                                const double* x,
                                const int* incx,
                                const double* beta,
                                      double* y,
                                const int* incy);

void  CVM_STD_CALL CGEMV       (const char* trans,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int transsz,
#endif
                                const int* m,
                                const int* n,
                                const std::complex<float>* alpha,
                                const std::complex<float>* a,
                                const int* lda,
                                const std::complex<float>* x,
                                const int* incx,
                                const std::complex<float>* beta,
                                      std::complex<float>* y,
                                const int* incy);
void  CVM_STD_CALL ZGEMV       (const char* trans,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int transsz,
#endif
                                const int* m,
                                const int* n,
                                const std::complex<double>* alpha,
                                const std::complex<double>* a,
                                const int* lda,
                                const std::complex<double>* x,
                                const int* incx,
                                const std::complex<double>* beta,
                                      std::complex<double>* y,
                                const int* incy);

void  CVM_STD_CALL SGBMV       (const char* trans,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int transsz,
#endif
                                const int* m,
                                const int* n,
                                const int* kl,
                                const int* ku,
                                const float* alpha,
                                const float* a,
                                const int* lda,
                                const float* x,
                                const int* incx,
                                const float* beta,
                                      float* y,
                                const int* incy);
void  CVM_STD_CALL DGBMV       (const char* trans,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int transsz,
#endif
                                const int* m,
                                const int* n,
                                const int* kl,
                                const int* ku,
                                const double* alpha,
                                const double* a,
                                const int* lda,
                                const double* x,
                                const int* incx,
                                const double* beta,
                                      double* y,
                                const int* incy);
void  CVM_STD_CALL CGBMV       (const char* trans,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int transsz,
#endif
                                const int* m,
                                const int* n,
                                const int* kl,
                                const int* ku,
                                const std::complex<float>* alpha,
                                const std::complex<float>* a,
                                const int* lda,
                                const std::complex<float>* x,
                                const int* incx,
                                const std::complex<float>* beta,
                                      std::complex<float>* y,
                                const int* incy);
void  CVM_STD_CALL ZGBMV       (const char* trans,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int transsz,
#endif
                                const int* m,
                                const int* n,
                                const int* kl,
                                const int* ku,
                                const std::complex<double>* alpha,
                                const std::complex<double>* a,
                                const int* lda,
                                const std::complex<double>* x,
                                const int* incx,
                                const std::complex<double>* beta,
                                      std::complex<double>* y,
                                const int* incy);

void  CVM_STD_CALL DGEMM       (const char* transa,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int transasz,
#endif
                                const char* transb,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int transbsz,
#endif
                                const int* m,
                                const int* n,
                                const int* k,
                                const double* alpha,
                                const double* a,
                                const int* lda,
                                const double* b,
                                const int* ldb,
                                const double* beta,
                                      double* c,
                                const int* ldc);

void  CVM_STD_CALL SGEMM       (const char* transa,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int transasz,
#endif
                                const char* transb,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int transbsz,
#endif
                                const int* m,
                                const int* n,
                                const int* k,
                                const float* alpha,
                                const float* a,
                                const int* lda,
                                const float* b,
                                const int* ldb,
                                const float* beta,
                                      float* c,
                                const int* ldc);

void  CVM_STD_CALL CGEMM       (const char* transa,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int transasz,
#endif
                                const char* transb,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int transbsz,
#endif
                                const int* m,
                                const int* n,
                                const int* k,
                                const std::complex<float>* alpha,
                                const std::complex<float>* a,
                                const int* lda,
                                const std::complex<float>* b,
                                const int* ldb,
                                const std::complex<float>* beta,
                                      std::complex<float>* c,
                                const int* ldc);
void  CVM_STD_CALL ZGEMM       (const char* transa,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int transasz,
#endif
                                const char* transb,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int transbsz,
#endif
                                const int* m,
                                const int* n,
                                const int* k,
                                const std::complex<double>* alpha,
                                const std::complex<double>* a,
                                const int* lda,
                                const std::complex<double>* b,
                                const int* ldb,
                                const std::complex<double>* beta,
                                      std::complex<double>* c,
                                const int* ldc);

void  CVM_FTN_CALL SMEXP       (const int* m,
                                const float* a,
                                const int* lda,
                                      float* ea,
                                const int* lde,
                                      float* r,
                                      int* ir,
                                      int* nr,
                                      int* ni,
                                      int* nq,
                                      int* j,
                                const int* issymm,                  // LOGICAL is 4-byte int
                                      float* work,
                                const int* lwork);
                                      
void  CVM_FTN_CALL DMEXP       (const int* m,
                                const double* a,
                                const int* lda,
                                      double* ea,
                                const int* lde,
                                      double* r,
                                      int* ir,
                                      int* nr,
                                      int* ni,
                                      int* nq,
                                      int* j,
                                const int* issymm,                  // LOGICAL is 4-byte int
                                      double* work,
                                const int* lwork);

void  CVM_FTN_CALL CMEXP       (const int* m,
                                const std::complex<float>* a,
                                const int* lda,
                                      std::complex<float>* ea,
                                const int* lde,
                                      std::complex<float>* r,
                                      int* ir,
                                      int* nr,
                                      int* ni,
                                      int* nq,
                                      int* j,
                                const int* issymm,                  // LOGICAL is 4-byte int
                                      std::complex<float>* work,
                                const int* lwork);

void  CVM_FTN_CALL ZMEXP       (const int* m,
                                const std::complex<double>* a,
                                const int* lda,
                                      std::complex<double>* ea,
                                const int* lde,
                                      std::complex<double>* r,
                                      int* ir,
                                      int* nr,
                                      int* ni,
                                      int* nq,
                                      int* j,
                                const int* issymm,                  // LOGICAL is 4-byte int
                                      std::complex<double>* work,
                                const int* lwork);

void  CVM_FTN_CALL SMEXPC      (const int* m,
                                const float* a,
                                const int* lda,
                                const float* tol,
                                      int* nr,
                                      int* ni,
                                      int* nq,
                                      int* j);

void  CVM_FTN_CALL DMEXPC      (const int* m,
                                const double* a,
                                const int* lda,
                                const double* tol,
                                      int* nr,
                                      int* ni,
                                      int* nq,
                                      int* j);

void  CVM_FTN_CALL CMEXPC      (const int* m,
                                const std::complex<float>* a,
                                const int* lda,
                                const float* tol,
                                      int* nr,
                                      int* ni,
                                      int* nq,
                                      int* j);

void  CVM_FTN_CALL ZMEXPC      (const int* m,
                                const std::complex<double>* a,
                                const int* lda,
                                const double* tol,
                                      int* nr,
                                      int* ni,
                                      int* nq,
                                      int* j);

void  CVM_STD_CALL SGETRF      (const int* m,
                                const int* n, 
                                      float* a,
                                const int* lda,
                                      int* ipiv, 
                                      int* info);
void  CVM_STD_CALL DGETRF      (const int* m,
                                const int* n,
                                      double* a,
                                const int* lda,
                                      int* ipiv,
                                      int* info);
void  CVM_STD_CALL CGETRF      (const int* m,
                                const int* n, 
                                      std::complex<float>* a,
                                const int* lda,
                                      int* ipiv, 
                                      int* info);
void  CVM_STD_CALL ZGETRF      (const int* m,
                                const int* n, 
                                      std::complex<double>* a,
                                const int* lda,
                                      int* ipiv, 
                                      int* info);

void  CVM_STD_CALL SGBTRF      (const int* m,
                                const int* n,
                                const int* kl,
                                const int* ku,
                                      float* a,
                                const int* lda,
                                      int* ipiv,
                                      int* info);
void  CVM_STD_CALL DGBTRF      (const int* m,
                                const int* n,
                                const int* kl,
                                const int* ku,
                                      double* a,
                                const int* lda,
                                      int* ipiv,
                                      int* info);
void  CVM_STD_CALL CGBTRF      (const int* m,
                                const int* n,
                                const int* kl,
                                const int* ku,
                                      std::complex<float>* a,
                                const int* lda,
                                      int* ipiv,
                                      int* info);
void  CVM_STD_CALL ZGBTRF      (const int* m,
                                const int* n,
                                const int* kl,
                                const int* ku,
                                      std::complex<double>* a,
                                const int* lda,
                                      int* ipiv,
                                      int* info);

void  CVM_STD_CALL SGETRS      (const char* trans,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int transsz,
#endif
                                const int* n,
                                const int* nrhs, 
                                const float* a,
                                const int* lda,
                                      int* ipiv, 
                                      float* b,
                                const int* ldb,
                                      int* info);
void  CVM_STD_CALL DGETRS      (const char* trans,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int transsz,
#endif
                                const int* n,
                                const int* nrhs, 
                                const double* a,
                                const int* lda,
                                      int* ipiv, 
                                      double* b,
                                const int* ldb,
                                      int* info);

void  CVM_STD_CALL SGBTRS      (const char* trans,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int transsz,
#endif
                                const int* n,
                                const int* kl,
                                const int* ku,
                                const int* nrhs,
                                const float* a,
                                const int* lda,
                                      int* ipiv,
                                      float* b,
                                const int* ldb,
                                      int* info);
void  CVM_STD_CALL DGBTRS      (const char* trans,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int transsz,
#endif
                                const int* n,
                                const int* kl,
                                const int* ku,
                                const int* nrhs,
                                const double* a,
                                const int* lda,
                                      int* ipiv,
                                      double* b,
                                const int* ldb,
                                      int* info);

void  CVM_STD_CALL SGERFS      (const char* trans,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int transsz,
#endif
                                const int* n,
                                const int* nrhs, 
                                const float* a,
                                const int* lda,
                                const float* af,
                                const int* ldaf,
                                      int* ipiv, 
                                const float* b,
                                const int* ldb,
                                      float* x,
                                const int* ldx,
                                      float* ferr,
                                      float* berr,
                                      float* work,
                                      int* iwork,
                                      int* info);
void  CVM_STD_CALL DGERFS      (const char* trans,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int transsz,
#endif
                                const int* n,
                                const int* nrhs, 
                                const double* a,
                                const int* lda,
                                const double* af,
                                const int* ldaf,
                                      int* ipiv, 
                                const double* b,
                                const int* ldb,
                                      double* x,
                                const int* ldx,
                                      double* ferr,
                                      double* berr,
                                      double* work,
                                      int* iwork,
                                      int* info);

void  CVM_STD_CALL SGBRFS      (const char* trans,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int transsz,
#endif
                                const int* n,
                                const int* kl,
                                const int* ku,
                                const int* nrhs,
                                const float* a,
                                const int* lda,
                                const float* af,
                                const int* ldaf,
                                      int* ipiv,
                                const float* b,
                                const int* ldb,
                                      float* x,
                                const int* ldx,
                                      float* ferr,
                                      float* berr,
                                      float* work,
                                      int* iwork,
                                      int* info);
void  CVM_STD_CALL DGBRFS      (const char* trans,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int transsz,
#endif
                                const int* n,
                                const int* kl,
                                const int* ku,
                                const int* nrhs,
                                const double* a,
                                const int* lda,
                                const double* af,
                                const int* ldaf,
                                      int* ipiv,
                                const double* b,
                                const int* ldb,
                                      double* x,
                                const int* ldx,
                                      double* ferr,
                                      double* berr,
                                      double* work,
                                      int* iwork,
                                      int* info);

void  CVM_STD_CALL CGETRS      (const char* trans,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int transsz,
#endif
                                const int* n,
                                const int* nrhs, 
                                const std::complex<float>* a,
                                const int* lda,
                                      int* ipiv, 
                                      std::complex<float>* b,
                                const int* ldb,
                                      int* info);
void  CVM_STD_CALL ZGETRS      (const char* trans,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int transsz,
#endif
                                const int* n,
                                const int* nrhs, 
                                const std::complex<double>* a,
                                const int* lda,
                                      int* ipiv, 
                                      std::complex<double>* b,
                                const int* ldb,
                                      int* info);

void  CVM_STD_CALL CGBTRS      (const char* trans,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int transsz,
#endif
                                const int* n,
                                const int* kl,
                                const int* ku,
                                const int* nrhs,
                                const std::complex<float>* a,
                                const int* lda,
                                      int* ipiv,
                                      std::complex<float>* b,
                                const int* ldb,
                                      int* info);
void  CVM_STD_CALL ZGBTRS      (const char* trans,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int transsz,
#endif
                                const int* n,
                                const int* kl,
                                const int* ku,
                                const int* nrhs,
                                const std::complex<double>* a,
                                const int* lda,
                                      int* ipiv,
                                      std::complex<double>* b,
                                const int* ldb,
                                      int* info);

void  CVM_STD_CALL CGERFS      (const char* trans,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int transsz,
#endif
                                const int* n,
                                const int* nrhs, 
                                const std::complex<float>* a,
                                const int* lda,
                                const std::complex<float>* af,
                                const int* ldaf,
                                      int* ipiv, 
                                const std::complex<float>* b,
                                const int* ldb,
                                      std::complex<float>* x,
                                const int* ldx,
                                      float* ferr,
                                      float* berr,
                                      std::complex<float>* work,
                                      float* rwork,
                                      int* info);
void  CVM_STD_CALL ZGERFS      (const char* trans,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int transsz,
#endif
                                const int* n,
                                const int* nrhs, 
                                const std::complex<double>* a,
                                const int* lda,
                                const std::complex<double>* af,
                                const int* ldaf,
                                      int* ipiv, 
                                const std::complex<double>* b,
                                const int* ldb,
                                      std::complex<double>* x,
                                const int* ldx,
                                      double* ferr,
                                      double* berr,
                                      std::complex<double>* work,
                                      double* rwork,
                                      int* info);

void  CVM_STD_CALL CGBRFS      (const char* trans,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int transsz,
#endif
                                const int* n,
                                const int* kl,
                                const int* ku,
                                const int* nrhs,
                                const std::complex<float>* a,
                                const int* lda,
                                const std::complex<float>* af,
                                const int* ldaf,
                                      int* ipiv,
                                const std::complex<float>* b,
                                const int* ldb,
                                      std::complex<float>* x,
                                const int* ldx,
                                      float* ferr,
                                      float* berr,
                                      std::complex<float>* work,
                                      float* rwork,
                                      int* info);
void  CVM_STD_CALL ZGBRFS      (const char* trans,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int transsz,
#endif
                                const int* n,
                                const int* kl,
                                const int* ku,
                                const int* nrhs,
                                const std::complex<double>* a,
                                const int* lda,
                                const std::complex<double>* af,
                                const int* ldaf,
                                      int* ipiv,
                                const std::complex<double>* b,
                                const int* ldb,
                                      std::complex<double>* x,
                                const int* ldx,
                                      double* ferr,
                                      double* berr,
                                      std::complex<double>* work,
                                      double* rwork,
                                      int* info);

void  CVM_STD_CALL SGETRI      (const int* n, 
                                      float* a,
                                const int* lda,
                                      int* ipiv, 
                                      float* work,
                                const int* lwork,
                                      int* info);
void  CVM_STD_CALL DGETRI      (const int* n, 
                                      double* a,
                                const int* lda,
                                      int* ipiv, 
                                      double* work,
                                const int* lwork,
                                      int* info);

void  CVM_STD_CALL CGETRI      (const int* n, 
                                      std::complex<float>* a,
                                const int* lda,
                                      int* ipiv, 
                                      std::complex<float>* work,
                                const int* lwork,
                                      int* info);
void  CVM_STD_CALL ZGETRI      (const int* n, 
                                      std::complex<double>* a,
                                const int* lda,
                                      int* ipiv, 
                                      std::complex<double>* work,
                                const int* lwork,
                                      int* info);

void  CVM_STD_CALL SGEBRD      (const int* m,
                                const int* n,
                                      float* a,
                                const int* lda,
                                      float* d,
                                      float* e,
                                      float* tauq,
                                      float* taup,
                                      float* work,
                                const int* lwork,
                                      int* info);
void  CVM_STD_CALL DGEBRD      (const int* m,
                                const int* n,
                                      double* a,
                                const int* lda,
                                      double* d,
                                      double* e,
                                      double* tauq,
                                      double* taup,
                                      double* work,
                                const int* lwork,
                                      int* info);

void  CVM_STD_CALL CGEBRD      (const int* m,
                                const int* n,
                                      std::complex<float>* a,
                                const int* lda,
                                      float* d,
                                      float* e,
                                      std::complex<float>* tauq,
                                      std::complex<float>* taup,
                                      std::complex<float>* work,
                                const int* lwork,
                                      int* info);
void  CVM_STD_CALL ZGEBRD      (const int* m,
                                const int* n,
                                      std::complex<double>* a,
                                const int* lda,
                                      double* d,
                                      double* e,
                                      std::complex<double>* tauq,
                                      std::complex<double>* taup,
                                      std::complex<double>* work,
                                const int* lwork,
                                      int* info);

void  CVM_STD_CALL SGBBRD      (const char* vect,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int vectsz,
#endif
                                const int* m,
                                const int* n,
                                const int* ncc,
                                const int* kl,
                                const int* ku,
                                      float* ab,
                                const int* ldab,
                                      float* d,
                                      float* e,
                                      float* q,
                                const int* ldq,
                                      float* pt,
                                const int* ldpt,
                                      float* c,
                                const int* ldc,
                                      float* work,
                                      int* info);

void  CVM_STD_CALL DGBBRD      (const char* vect,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int vectsz,
#endif
                                const int* m,
                                const int* n,
                                const int* ncc,
                                const int* kl,
                                const int* ku,
                                      double* ab,
                                const int* ldab,
                                      double* d,
                                      double* e,
                                      double* q,
                                const int* ldq,
                                      double* pt,
                                const int* ldpt,
                                      double* c,
                                const int* ldc,
                                      double* work,
                                      int* info);

void  CVM_STD_CALL CGBBRD      (const char* vect,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int vectsz,
#endif
                                const int* m,
                                const int* n,
                                const int* ncc,
                                const int* kl,
                                const int* ku,
                                      std::complex<float>* ab,
                                const int* ldab,
                                      float* d,
                                      float* e,
                                      std::complex<float>* q,
                                const int* ldq,
                                      std::complex<float>* pt,
                                const int* ldpt,
                                      std::complex<float>* c,
                                const int* ldc,
                                      std::complex<float>* work,
                                      float* rwork,
                                      int* info);

void  CVM_STD_CALL ZGBBRD      (const char* vect,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int vectsz,
#endif
                                const int* m,
                                const int* n,
                                const int* ncc,
                                const int* kl,
                                const int* ku,
                                      std::complex<double>* ab,
                                const int* ldab,
                                      double* d,
                                      double* e,
                                      std::complex<double>* q,
                                const int* ldq,
                                      std::complex<double>* pt,
                                const int* ldpt,
                                      std::complex<double>* c,
                                const int* ldc,
                                      std::complex<double>* work,
                                      double* rwork,
                                      int* info);

void  CVM_STD_CALL SORGBR      (const char* uplo,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int uplosz,
#endif
                                const int* m,
                                const int* n,
                                const int* k,
                                      float* a,
                                const int* lda,
                                      float* tau,
                                      float* work,
                                const int* lwork,
                                      int* info);
void  CVM_STD_CALL DORGBR      (const char* uplo,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int uplosz,
#endif
                                const int* m,
                                const int* n,
                                const int* k,
                                      double* a,
                                const int* lda,
                                      double* tau,
                                      double* work,
                                const int* lwork,
                                      int* info);

void  CVM_STD_CALL CUNGBR      (const char* uplo,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int uplosz,
#endif
                                const int* m,
                                const int* n,
                                const int* k,
                                      std::complex<float>* a,
                                const int* lda,
                                      std::complex<float>* tau,
                                      std::complex<float>* work,
                                const int* lwork,
                                      int* info);
void  CVM_STD_CALL ZUNGBR      (const char* uplo,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int uplosz,
#endif
                                const int* m,
                                const int* n,
                                const int* k,
                                      std::complex<double>* a,
                                const int* lda,
                                      std::complex<double>* tau,
                                      std::complex<double>* work,
                                const int* lwork,
                                      int* info);

void  CVM_STD_CALL SBDSQR      (const char* uplo,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int uplosz,
#endif
                                const int* n,
                                const int* ncvt,
                                const int* nru,
                                const int* ncc,
                                      float* d,
                                      float* e,
                                      float* vt,
                                const int* ldvt,
                                      float* u,
                                const int* ldu,
                                      float* c,
                                const int* ldc,
                                      float* work,
                                      int* info);
void  CVM_STD_CALL DBDSQR      (const char* uplo,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int uplosz,
#endif
                                const int* n,
                                const int* ncvt,
                                const int* nru,
                                const int* ncc,
                                      double* d,
                                      double* e,
                                      double* vt,
                                const int* ldvt,
                                      double* u,
                                const int* ldu,
                                      double* c,
                                const int* ldc,
                                      double* work,
                                      int* info);

void  CVM_STD_CALL CBDSQR      (const char* uplo,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int uplosz,
#endif
                                const int* n,
                                const int* ncvt,
                                const int* nru,
                                const int* ncc,
                                      float* d,
                                      float* e,
                                      std::complex<float>* vt,
                                const int* ldvt,
                                      std::complex<float>* u,
                                const int* ldu,
                                      std::complex<float>* c,
                                const int* ldc,
                                      float* work,
                                      int* info);
void  CVM_STD_CALL ZBDSQR      (const char* uplo,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int uplosz,
#endif
                                const int* n,
                                const int* ncvt,
                                const int* nru,
                                const int* ncc,
                                      double* d,
                                      double* e,
                                      std::complex<double>* vt,
                                const int* ldvt,
                                      std::complex<double>* u,
                                const int* ldu,
                                      std::complex<double>* c,
                                const int* ldc,
                                      double* work,
                                      int* info);

void  CVM_STD_CALL SGEBAL      (const char* job,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int jobsz,
#endif
                                const int* n,
                                      float* a,
                                const int* lda,

                                      int* ilo,
                                      int* ihi,
                                      float* scale,
                                      int* info);
void  CVM_STD_CALL DGEBAL      (const char* job,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int jobsz,
#endif
                                const int* n,
                                      double* a,
                                const int* lda,

                                      int* ilo,
                                      int* ihi,
                                      double* scale,
                                      int* info);

void  CVM_STD_CALL CGEBAL      (const char* job,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int jobsz,
#endif
                                const int* n,
                                      std::complex<float>* a,
                                const int* lda,
                                      int* ilo,
                                      int* ihi,
                                      float* scale,
                                      int* info);
void  CVM_STD_CALL ZGEBAL      (const char* job,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int jobsz,
#endif
                                const int* n,
                                      std::complex<double>* a,
                                const int* lda,
                                      int* ilo,
                                      int* ihi,
                                      double* scale,
                                      int* info);

void  CVM_STD_CALL SGEHRD      (const int* n,
                                const int* ilo,
                                const int* ihi,
                                      float* a,
                                const int* lda,
                                      float* tau,
                                      float* work,
                                const int* lwork,
                                      int* info);
void  CVM_STD_CALL DGEHRD      (const int* n,
                                const int* ilo,
                                const int* ihi,
                                      double* a,
                                const int* lda,
                                      double* tau,
                                      double* work,
                                const int* lwork,
                                      int* info);

void  CVM_STD_CALL CGEHRD      (const int* n,
                                const int* ilo,
                                const int* ihi,
                                      std::complex<float>* a,
                                const int* lda,
                                      std::complex<float>* tau,
                                      std::complex<float>* work,
                                const int* lwork,
                                      int* info);
void  CVM_STD_CALL ZGEHRD      (const int* n,
                                const int* ilo,
                                const int* ihi,
                                      std::complex<double>* a,
                                const int* lda,
                                      std::complex<double>* tau,
                                      std::complex<double>* work,
                                const int* lwork,
                                      int* info);

void  CVM_STD_CALL SORGHR      (const int* n,
                                const int* ilo,
                                const int* ihi,
                                      float* a,
                                const int* lda,
                                      float* tau,
                                      float* work,
                                const int* lwork,
                                      int* info);
void  CVM_STD_CALL DORGHR      (const int* n,
                                const int* ilo,
                                const int* ihi,
                                      double* a,
                                const int* lda,
                                      double* tau,
                                      double* work,
                                const int* lwork,
                                      int* info);

void  CVM_STD_CALL CUNGHR      (const int* n,
                                const int* ilo,
                                const int* ihi,
                                      std::complex<float>* a,
                                const int* lda,
                                      std::complex<float>* tau,
                                      std::complex<float>* work,
                                const int* lwork,
                                      int* info);
void  CVM_STD_CALL ZUNGHR      (const int* n,
                                const int* ilo,
                                const int* ihi,
                                      std::complex<double>* a,
                                const int* lda,
                                      std::complex<double>* tau,
                                      std::complex<double>* work,
                                const int* lwork,
                                      int* info);

void  CVM_STD_CALL SHSEQR      (const char* job,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int jobsz,
#endif
                                const char* compz,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int compzsz,
#endif
                                const int* n,
                                const int* ilo,
                                const int* ihi,
                                      float* h,
                                const int* ldh,
                                      float* wr,
                                      float* wi,
                                      float* z,
                                const int* ldz,
                                      float* work,
                                const int* lwork,
                                      int* info);
void  CVM_STD_CALL DHSEQR      (const char* job,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int jobsz,
#endif
                                const char* compz,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int compzsz,
#endif
                                const int* n,
                                const int* ilo,
                                const int* ihi,
                                      double* h,
                                const int* ldh,
                                      double* wr,
                                      double* wi,
                                      double* z,
                                const int* ldz,
                                      double* work,
                                const int* lwork,
                                      int* info);

void  CVM_STD_CALL CHSEQR      (const char* job,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int jobsz,
#endif
                                const char* compz,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int compzsz,
#endif
                                const int* n,
                                const int* ilo,
                                const int* ihi,
                                      std::complex<float>* h,
                                const int* ldh,
                                      std::complex<float>* w,
                                      std::complex<float>* z,
                                const int* ldz,
                                      std::complex<float>* work,
                                const int* lwork,
                                      int* info);
void  CVM_STD_CALL ZHSEQR      (const char* job,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int jobsz,
#endif
                                const char* compz,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int compzsz,
#endif
                                const int* n,
                                const int* ilo,
                                const int* ihi,
                                      std::complex<double>* h,
                                const int* ldh,
                                      std::complex<double>* w,
                                      std::complex<double>* z,
                                const int* ldz,
                                      std::complex<double>* work,
                                const int* lwork,
                                      int* info);

void  CVM_STD_CALL STREVC      (const char* job,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int jobsz,
#endif
                                const char* howmny,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int howmnysz,
#endif
                                const int* select,                                  // CAUTION! We assume that LOGICAL datatype is 4-byte long
                                const int* n,
                                      float* t,
                                const int* ldt,
                                      float* vl,
                                const int* ldvl,
                                      float* vr,
                                const int* ldvr,
                                const int* mm,
                                      int* m,
                                      float* work,
                                      int* info);
void  CVM_STD_CALL DTREVC      (const char* job,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int jobsz,
#endif
                                const char* howmny,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int howmnysz,
#endif
                                const int* select,                                  // CAUTION! We assume that LOGICAL datatype is 4-byte long
                                const int* n,
                                      double* t,
                                const int* ldt,
                                      double* vl,
                                const int* ldvl,
                                      double* vr,
                                const int* ldvr,
                                const int* mm,
                                      int* m,
                                      double* work,
                                      int* info);

void  CVM_STD_CALL CTREVC      (const char* job,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int jobsz,
#endif
                                const char* howmny,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int howmnysz,
#endif
                                const int* select,                                  // CAUTION! We assume that LOGICAL datatype is 4-byte long
                                const int* n,
                                      std::complex<float>* t,
                                const int* ldt,
                                      std::complex<float>* vl,
                                const int* ldvl,
                                      std::complex<float>* vr,
                                const int* ldvr,
                                const int* mm,
                                      int* m,
                                      std::complex<float>* work,
                                      float* rwork,
                                      int* info);
void  CVM_STD_CALL ZTREVC      (const char* job,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int jobsz,
#endif
                                const char* howmny,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int howmnysz,
#endif
                                const int* select,                                  // CAUTION! We assume that LOGICAL datatype is 4-byte long
                                const int* n,
                                      std::complex<double>* t,
                                const int* ldt,
                                      std::complex<double>* vl,
                                const int* ldvl,
                                      std::complex<double>* vr,
                                const int* ldvr,
                                const int* mm,
                                      int* m,
                                      std::complex<double>* work,
                                      double* rwork,
                                      int* info);

void  CVM_STD_CALL SGEBAK      (const char* job,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int jobsz,
#endif
                                const char* side,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int sidesz,
#endif
                                const int* n,
                                const int* ilo,
                                const int* ihi,
                                      float* scale,
                                const int* m,
                                      float* v,
                                const int* ldv,
                                      int* info);
void  CVM_STD_CALL DGEBAK      (const char* job,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int jobsz,
#endif
                                const char* side,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int sidesz,
#endif
                                const int* n,
                                const int* ilo,
                                const int* ihi,
                                      double* scale,
                                const int* m,
                                      double* v,
                                const int* ldv,
                                      int* info);

void  CVM_STD_CALL CGEBAK      (const char* job,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int jobsz,
#endif
                                const char* side,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int sidesz,
#endif
                                const int* n,
                                const int* ilo,
                                const int* ihi,
                                      float* scale,
                                const int* m,
                                      std::complex<float>* v,
                                const int* ldv,
                                      int* info);
void  CVM_STD_CALL ZGEBAK      (const char* job,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int jobsz,
#endif
                                const char* side,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int sidesz,
#endif
                                const int* n,
                                const int* ilo,
                                const int* ihi,
                                      double* scale,
                                const int* m,
                                      std::complex<double>* v,
                                const int* ldv,
                                      int* info);

void  CVM_STD_CALL SGECON      (const char* norm,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int normsz,
#endif
                                const int* n,
                                      float* a,                                     // const
                                const int* lda,
                                const float* anorm,
                                      float* rcond,
                                      float* work,
                                      int* iwork,
                                      int* info);
void  CVM_STD_CALL DGECON      (const char* norm,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int normsz,
#endif
                                const int* n,
                                      double* a,                                    // const
                                const int* lda,
                                const double* anorm,
                                      double* rcond,
                                      double* work,
                                      int* iwork,
                                      int* info);

void  CVM_STD_CALL CGECON      (const char* norm,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int normsz,
#endif
                                const int* n,
                                      std::complex<float>* a,                       // const
                                const int* lda,
                                const float* anorm,
                                      float* rcond,
                                      std::complex<float>* work,
                                      float* rwork,
                                      int* info);
void  CVM_STD_CALL ZGECON      (const char* norm,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int normsz,
#endif
                                const int* n,
                                      std::complex<double>* a,                      // const
                                const int* lda,
                                const double* anorm,
                                      double* rcond,
                                      std::complex<double>* work,
                                      double* rwork,
                                      int* info);

void  CVM_STD_CALL SSPMV       (const char* uplo,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int uplosz,
#endif
                                const int* n,
                                const float* alpha,
                                const float* ap,
                                const float* x,
                                const int* incx,
                                const float* beta,
                                      float* y,
                                const int* incy);

void  CVM_STD_CALL DSPMV       (const char* uplo,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int uplosz,
#endif
                                const int* n,
                                const double* alpha,
                                const double* ap,
                                const double* x,
                                const int* incx,
                                const double* beta,
                                      double* y,
                                const int* incy);

void  CVM_STD_CALL SSYMM       (const char* side,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int sidesz,
#endif
                                const char* uplo,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int uplosz,
#endif
                                const int* m,
                                const int* n,
                                const float* alpha,
                                const float* a,
                                const int* lda,
                                const float* b,
                                const int* ldb,
                                const float* beta,
                                      float* c,
                                const int* ldc);

void  CVM_STD_CALL DSYMM       (const char* side,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int sidesz,
#endif
                                const char* uplo,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int uplosz,
#endif
                                const int* m,
                                const int* n,
                                const double* alpha,
                                const double* a,
                                const int* lda,
                                const double* b,
                                const int* ldb,
                                const double* beta,
                                      double* c,
                                const int* ldc);

void  CVM_STD_CALL CSYMM       (const char* side,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int sidesz,
#endif
                                const char* uplo,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int uplosz,
#endif
                                const int* m,
                                const int* n,
                                const std::complex<float>* alpha,
                                const std::complex<float>* a,
                                const int* lda,
                                const std::complex<float>* b,
                                const int* ldb,
                                const std::complex<float>* beta,
                                      std::complex<float>* c,
                                const int* ldc);

void  CVM_STD_CALL ZSYMM       (const char* side,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int sidesz,
#endif
                                const char* uplo,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int uplosz,
#endif
                                const int* m,
                                const int* n,
                                const std::complex<double>* alpha,
                                const std::complex<double>* a,
                                const int* lda,
                                const std::complex<double>* b,
                                const int* ldb,
                                const std::complex<double>* beta,
                                      std::complex<double>* c,
                                const int* ldc);

void  CVM_STD_CALL CHEMM       (const char* side,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int sidesz,
#endif
                                const char* uplo,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int uplosz,
#endif
                                const int* m,
                                const int* n,
                                const std::complex<float>* alpha,
                                const std::complex<float>* a,
                                const int* lda,
                                const std::complex<float>* b,
                                const int* ldb,
                                const std::complex<float>* beta,
                                      std::complex<float>* c,
                                const int* ldc);

void  CVM_STD_CALL ZHEMM       (const char* side,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int sidesz,
#endif
                                const char* uplo,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int uplosz,
#endif
                                const int* m,
                                const int* n,
                                const std::complex<double>* alpha,
                                const std::complex<double>* a,
                                const int* lda,
                                const std::complex<double>* b,
                                const int* ldb,
                                const std::complex<double>* beta,
                                      std::complex<double>* c,
                                const int* ldc);

void  CVM_STD_CALL SPOTRF      (const char* uplo,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int uplosz,
#endif
                                const int* n,
                                      float* a,
                                const int* lda,
                                      int* info);

void  CVM_STD_CALL DPOTRF      (const char* uplo,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int uplosz,
#endif
                                const int* n,
                                      double* a,
                                const int* lda,
                                      int* info);

void  CVM_STD_CALL CPOTRF      (const char* uplo,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int uplosz,
#endif
                                const int* n,
                                      std::complex<float>* a,
                                const int* lda,
                                      int* info);

void  CVM_STD_CALL ZPOTRF      (const char* uplo,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int uplosz,
#endif
                                const int* n,
                                      std::complex<double>* a,
                                const int* lda,
                                      int* info);

void  CVM_STD_CALL SSYTRF      (const char* uplo,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int uplosz,
#endif
                                const int* n,
                                      float* a,
                                const int* lda,
                                      int* ipiv,
                                      float* work,
                                const int* lwork,
                                      int* info);

void  CVM_STD_CALL DSYTRF      (const char* uplo,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int uplosz,
#endif
                                const int* n,
                                      double* a,
                                const int* lda,
                                      int* ipiv,
                                      double* work,
                                const int* lwork,
                                      int* info);

void  CVM_STD_CALL CSYTRF      (const char* uplo,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int uplosz,
#endif
                                const int* n,
                                      std::complex<float>* a,
                                const int* lda,
                                      int* ipiv,
                                      std::complex<float>* work,
                                const int* lwork,
                                      int* info);

void  CVM_STD_CALL ZSYTRF      (const char* uplo,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int uplosz,
#endif
                                const int* n,
                                      std::complex<double>* a,
                                const int* lda,
                                      int* ipiv,
                                      std::complex<double>* work,
                                const int* lwork,
                                      int* info);

void  CVM_STD_CALL CHETRF      (const char* uplo,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int uplosz,
#endif
                                const int* n,
                                      std::complex<float>* a,
                                const int* lda,
                                      int* ipiv,
                                      std::complex<float>* work,
                                const int* lwork,
                                      int* info);

void  CVM_STD_CALL ZHETRF      (const char* uplo,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int uplosz,
#endif
                                const int* n,
                                      std::complex<double>* a,
                                const int* lda,
                                      int* ipiv,
                                      std::complex<double>* work,
                                const int* lwork,
                                      int* info);

void  CVM_STD_CALL SPOTRS      (const char* uplo,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int uplosz,
#endif
                                const int* n,
                                const int* nrhs,
                                const float* a,
                                const int* lda,
                                      float* b,
                                const int* ldb,
                                      int* info);

void  CVM_STD_CALL DPOTRS      (const char* uplo,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int uplosz,
#endif
                                const int* n,
                                const int* nrhs,
                                const double* a,
                                const int* lda,
                                      double* b,
                                const int* ldb,
                                      int* info);

void  CVM_STD_CALL CPOTRS      (const char* uplo,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int uplosz,
#endif
                                const int* n,
                                const int* nrhs,
                                const std::complex<float>* a,
                                const int* lda,
                                      std::complex<float>* b,
                                const int* ldb,
                                      int* info);

void  CVM_STD_CALL ZPOTRS      (const char* uplo,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int uplosz,
#endif
                                const int* n,
                                const int* nrhs,
                                const std::complex<double>* a,
                                const int* lda,
                                      std::complex<double>* b,
                                const int* ldb,
                                      int* info);

void  CVM_STD_CALL SPORFS      (const char* uplo,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int uplosz,
#endif
                                const int* n,
                                const int* nrhs,
                                const float* a,
                                const int* lda,
                                const float* af,
                                const int* ldaf,
                                const float* b,
                                const int* ldb,
                                      float* x,
                                const int* ldx,
                                      float* ferr,
                                      float* berr,
                                      float* work,
                                      int* iwork,
                                      int* info);

void  CVM_STD_CALL DPORFS      (const char* uplo,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int uplosz,
#endif
                                const int* n,
                                const int* nrhs,
                                const double* a,
                                const int* lda,
                                const double* af,
                                const int* ldaf,
                                const double* b,
                                const int* ldb,
                                      double* x,
                                const int* ldx,
                                      double* ferr,
                                      double* berr,
                                      double* work,
                                      int* iwork,
                                      int* info);

void  CVM_STD_CALL CPORFS      (const char* uplo,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int uplosz,
#endif
                                const int* n,
                                const int* nrhs,
                                const std::complex<float>* a,
                                const int* lda,
                                const std::complex<float>* af,
                                const int* ldaf,
                                const std::complex<float>* b,
                                const int* ldb,
                                      std::complex<float>* x,
                                const int* ldx,
                                      float* ferr,
                                      float* berr,
                                      std::complex<float>* work,
                                      float* rwork,
                                      int* info);

void  CVM_STD_CALL ZPORFS      (const char* uplo,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int uplosz,
#endif
                                const int* n,
                                const int* nrhs,
                                const std::complex<double>* a,
                                const int* lda,
                                const std::complex<double>* af,
                                const int* ldaf,
                                const std::complex<double>* b,
                                const int* ldb,
                                      std::complex<double>* x,
                                const int* ldx,
                                      double* ferr,
                                      double* berr,
                                      std::complex<double>* work,
                                      double* rwork,
                                      int* info);

void  CVM_STD_CALL SSYTRS      (const char* uplo,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int uplosz,
#endif
                                const int* n,
                                const int* nrhs,
                                const float* a,
                                const int* lda,
                                const int* ipiv,
                                      float* b,
                                const int* ldb,
                                      int* info);

void  CVM_STD_CALL DSYTRS      (const char* uplo,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int uplosz,
#endif
                                const int* n,
                                const int* nrhs,
                                const double* a,
                                const int* lda,
                                const int* ipiv,
                                      double* b,
                                const int* ldb,
                                      int* info);

void  CVM_STD_CALL CSYTRS      (const char* uplo,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int uplosz,
#endif
                                const int* n,
                                const int* nrhs,
                                const std::complex<float>* a,
                                const int* lda,
                                const int* ipiv,
                                      std::complex<float>* b,
                                const int* ldb,
                                      int* info);

void  CVM_STD_CALL ZSYTRS      (const char* uplo,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int uplosz,
#endif
                                const int* n,
                                const int* nrhs,
                                const std::complex<double>* a,
                                const int* lda,
                                const int* ipiv,
                                      std::complex<double>* b,
                                const int* ldb,
                                      int* info);

void  CVM_STD_CALL CHETRS      (const char* uplo,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int uplosz,
#endif
                                const int* n,
                                const int* nrhs,
                                const std::complex<float>* a,
                                const int* lda,
                                const int* ipiv,
                                      std::complex<float>* b,
                                const int* ldb,
                                      int* info);

void  CVM_STD_CALL ZHETRS      (const char* uplo,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int uplosz,
#endif
                                const int* n,
                                const int* nrhs,
                                const std::complex<double>* a,
                                const int* lda,
                                const int* ipiv,
                                      std::complex<double>* b,
                                const int* ldb,
                                      int* info);

void  CVM_STD_CALL SSYRFS      (const char* uplo,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int uplosz,
#endif
                                const int* n,
                                const int* nrhs,
                                const float* a,
                                const int* lda,
                                const float* af,
                                const int* ldaf,
                                const int* ipiv,
                                const float* b,
                                const int* ldb,
                                      float* x,
                                const int* ldx,
                                      float* ferr,
                                      float* berr,
                                      float* work,
                                      int* iwork,
                                      int* info);

void  CVM_STD_CALL DSYRFS      (const char* uplo,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int uplosz,
#endif
                                const int* n,
                                const int* nrhs,
                                const double* a,
                                const int* lda,
                                const double* af,
                                const int* ldaf,
                                const int* ipiv,
                                const double* b,
                                const int* ldb,
                                      double* x,
                                const int* ldx,
                                      double* ferr,
                                      double* berr,
                                      double* work,
                                      int* iwork,
                                      int* info);

void  CVM_STD_CALL CSYRFS      (const char* uplo,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int uplosz,
#endif
                                const int* n,
                                const int* nrhs,
                                const std::complex<float>* a,
                                const int* lda,
                                const std::complex<float>* af,
                                const int* ldaf,
                                const int* ipiv,
                                const std::complex<float>* b,
                                const int* ldb,
                                      std::complex<float>* x,
                                const int* ldx,
                                      float* ferr,
                                      float* berr,
                                      std::complex<float>* work,
                                      float* rwork,
                                      int* info);

void  CVM_STD_CALL ZSYRFS      (const char* uplo,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int uplosz,
#endif
                                const int* n,
                                const int* nrhs,
                                const std::complex<double>* a,
                                const int* lda,
                                const std::complex<double>* af,
                                const int* ldaf,
                                const int* ipiv,
                                const std::complex<double>* b,
                                const int* ldb,
                                      std::complex<double>* x,
                                const int* ldx,
                                      double* ferr,
                                      double* berr,
                                      std::complex<double>* work,
                                      double* rwork,
                                      int* info);

void  CVM_STD_CALL CHERFS      (const char* uplo,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int uplosz,
#endif
                                const int* n,
                                const int* nrhs,
                                const std::complex<float>* a,
                                const int* lda,
                                const std::complex<float>* af,
                                const int* ldaf,
                                const int* ipiv,
                                const std::complex<float>* b,
                                const int* ldb,
                                      std::complex<float>* x,
                                const int* ldx,
                                      float* ferr,
                                      float* berr,
                                      std::complex<float>* work,
                                      float* rwork,
                                      int* info);

void  CVM_STD_CALL ZHERFS      (const char* uplo,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int uplosz,
#endif
                                const int* n,
                                const int* nrhs,
                                const std::complex<double>* a,
                                const int* lda,
                                const std::complex<double>* af,
                                const int* ldaf,
                                const int* ipiv,
                                const std::complex<double>* b,
                                const int* ldb,
                                      std::complex<double>* x,
                                const int* ldx,
                                      double* ferr,
                                      double* berr,
                                      std::complex<double>* work,
                                      double* rwork,
                                      int* info);

void  CVM_STD_CALL SPOTRI      (const char* uplo,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int uplosz,
#endif
                                const int* n,
                                      float* a,
                                const int* lda,
                                      int* info);

void  CVM_STD_CALL DPOTRI      (const char* uplo,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int uplosz,
#endif
                                const int* n,
                                      double* a,
                                const int* lda,
                                      int* info);

void  CVM_STD_CALL CPOTRI      (const char* uplo,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int uplosz,
#endif
                                const int* n,
                                      std::complex<float>* a,
                                const int* lda,
                                      int* info);

void  CVM_STD_CALL ZPOTRI      (const char* uplo,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int uplosz,
#endif
                                const int* n,
                                      std::complex<double>* a,
                                const int* lda,
                                      int* info);

void  CVM_STD_CALL SSYTRI      (const char* uplo,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int uplosz,
#endif
                                const int* n,
                                      float* a,
                                const int* lda,
                                const int* ipiv,
                                      float* work,
                                      int* info);

void  CVM_STD_CALL DSYTRI      (const char* uplo,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int uplosz,
#endif
                                const int* n,
                                      double* a,
                                const int* lda,
                                const int* ipiv,
                                      double* work,
                                      int* info);

void  CVM_STD_CALL CSYTRI      (const char* uplo,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int uplosz,
#endif
                                const int* n,
                                      std::complex<float>* a,
                                const int* lda,
                                const int* ipiv,
                                      std::complex<float>* work,
                                      int* info);

void  CVM_STD_CALL ZSYTRI      (const char* uplo,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int uplosz,
#endif
                                const int* n,
                                      std::complex<double>* a,
                                const int* lda,
                                const int* ipiv,
                                      std::complex<double>* work,
                                      int* info);

void  CVM_STD_CALL CHETRI      (const char* uplo,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int uplosz,
#endif
                                const int* n,
                                      std::complex<float>* a,
                                const int* lda,
                                const int* ipiv,
                                      std::complex<float>* work,
                                      int* info);

void  CVM_STD_CALL ZHETRI      (const char* uplo,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int uplosz,
#endif
                                const int* n,
                                      std::complex<double>* a,
                                const int* lda,
                                const int* ipiv,
                                      std::complex<double>* work,
                                      int* info);

void  CVM_STD_CALL SSYEVD      (const char* job,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int jobsz,
#endif
                                const char* uplo,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int uplosz,
#endif
                                const int* n,
                                      float* a,
                                const int* lda,
                                      float* w,
                                      float* work,
                                      int* lwork,
                                      int* iwork,
                                      int* liwork,
                                      int* info);

void  CVM_STD_CALL DSYEVD      (const char* job,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int jobsz,
#endif
                                const char* uplo,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int uplosz,
#endif
                                const int* n,
                                      double* a,
                                const int* lda,
                                      double* w,
                                      double* work,
                                      int* lwork,
                                      int* iwork,
                                      int* liwork,
                                      int* info);

void  CVM_STD_CALL CHEEVD      (const char* job,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int jobsz,
#endif
                                const char* uplo,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int uplosz,
#endif
                                const int* n,
                                      std::complex<float>* a,
                                const int* lda,
                                      float* w,
                                      std::complex<float>* work,
                                      int* lwork,
                                      float* rwork,
                                      int* lrwork,
                                      int* iwork,
                                      int* liwork,
                                      int* info);

void  CVM_STD_CALL ZHEEVD      (const char* job,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int jobsz,
#endif
                                const char* uplo,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int uplosz,
#endif
                                const int* n,
                                      std::complex<double>* a,
                                const int* lda,
                                      double* w,
                                      std::complex<double>* work,
                                      int* lwork,
                                      double* rwork,
                                      int* lrwork,
                                      int* iwork,
                                      int* liwork,
                                      int* info);

void  CVM_STD_CALL SPOEQU      (const int* n,
                                const float* a,
                                const int* lda,
                                      float* s,
                                      float* scond,
                                      float* amax,
                                      int* info);

void  CVM_STD_CALL DPOEQU      (const int* n,
                                const double* a,
                                const int* lda,
                                      double* s,
                                      double* scond,
                                      double* amax,
                                      int* info);

void  CVM_STD_CALL CPOEQU      (const int* n,
                                const std::complex<float>* a,
                                const int* lda,
                                      float* s,
                                      float* scond,
                                      float* amax,
                                      int* info);

void  CVM_STD_CALL ZPOEQU      (const int* n,
                                const std::complex<double>* a,
                                const int* lda,
                                      double* s,
                                      double* scond,
                                      double* amax,
                                      int* info);

void  CVM_STD_CALL SSYMV       (const char* uplo,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int uplosz,
#endif
                                const int* n,
                                const float* alpha,
                                const float* a,
                                const int* lda,
                                const float* x,
                                const int* incx,
                                const float* beta,
                                      float* y,
                                const int* incy);

void  CVM_STD_CALL DSYMV       (const char* uplo,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int uplosz,
#endif
                                const int* n,
                                const double* alpha,
                                const double* a,
                                const int* lda,
                                const double* x,
                                const int* incx,
                                const double* beta,
                                      double* y,
                                const int* incy);

void  CVM_STD_CALL CHEMV       (const char* uplo,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int uplosz,
#endif
                                const int* n,
                                const std::complex<float>* alpha,
                                const std::complex<float>* a,
                                const int* lda,
                                const std::complex<float>* x,
                                const int* incx,
                                const std::complex<float>* beta,
                                      std::complex<float>* y,
                                const int* incy);

void  CVM_STD_CALL ZHEMV       (const char* uplo,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int uplosz,
#endif
                                const int* n,
                                const std::complex<double>* alpha,
                                const std::complex<double>* a,
                                const int* lda,
                                const std::complex<double>* x,
                                const int* incx,
                                const std::complex<double>* beta,
                                      std::complex<double>* y,
                                const int* incy);


void  CVM_STD_CALL SSYRK       (const char* uplo,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int uplosz,
#endif
                                const char* trans,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int transsz,
#endif
                                const int* n,
                                const int* k,
                                const float* alpha,
                                const float* a,
                                const int* lda,
                                const float* beta,
                                      float* c,
                                const int* ldc);

void  CVM_STD_CALL DSYRK       (const char* uplo,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int uplosz,
#endif
                                const char* trans,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int transsz,
#endif
                                const int* n,
                                const int* k,
                                const double* alpha,
                                const double* a,
                                const int* lda,
                                const double* beta,
                                      double* c,
                                const int* ldc);

void  CVM_STD_CALL CSYRK       (const char* uplo,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int uplosz,
#endif
                                const char* trans,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int transsz,
#endif
                                const int* n,
                                const int* k,
                                const std::complex<float>* alpha,
                                const std::complex<float>* a,
                                const int* lda,
                                const std::complex<float>* beta,
                                      std::complex<float>* c,
                                const int* ldc);

void  CVM_STD_CALL ZSYRK       (const char* uplo,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int uplosz,
#endif
                                const char* trans,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int transsz,
#endif
                                const int* n,
                                const int* k,
                                const std::complex<double>* alpha,
                                const std::complex<double>* a,
                                const int* lda,
                                const std::complex<double>* beta,
                                      std::complex<double>* c,
                                const int* ldc);

void  CVM_STD_CALL CHERK       (const char* uplo,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int uplosz,
#endif
                                const char* trans,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int transsz,
#endif
                                const int* n,
                                const int* k,
                                const std::complex<float>* alpha,
                                const std::complex<float>* a,
                                const int* lda,
                                const std::complex<float>* beta,
                                      std::complex<float>* c,
                                const int* ldc);

void  CVM_STD_CALL ZHERK       (const char* uplo,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int uplosz,
#endif
                                const char* trans,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int transsz,
#endif
                                const int* n,
                                const int* k,
                                const std::complex<double>* alpha,
                                const std::complex<double>* a,
                                const int* lda,
                                const std::complex<double>* beta,
                                      std::complex<double>* c,
                                const int* ldc);

void  CVM_STD_CALL SSYR2K      (const char* uplo,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int uplosz,
#endif
                                const char* trans,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int transsz,
#endif
                                const int* n,
                                const int* k,
                                const float* alpha,
                                const float* a,
                                const int* lda,
                                const float* b,
                                const int* ldb,
                                const float* beta,
                                      float* c,
                                const int* ldc);

void  CVM_STD_CALL DSYR2K      (const char* uplo,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int uplosz,
#endif
                                const char* trans,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int transsz,
#endif
                                const int* n,
                                const int* k,
                                const double* alpha,
                                const double* a,
                                const int* lda,
                                const double* b,
                                const int* ldb,
                                const double* beta,
                                      double* c,
                                const int* ldc);

void  CVM_STD_CALL CSYR2K      (const char* uplo,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int uplosz,
#endif
                                const char* trans,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int transsz,
#endif
                                const int* n,
                                const int* k,
                                const std::complex<float>* alpha,
                                const std::complex<float>* a,
                                const int* lda,
                                const std::complex<float>* b,
                                const int* ldb,
                                const std::complex<float>* beta,
                                      std::complex<float>* c,
                                const int* ldc);

void  CVM_STD_CALL ZSYR2K      (const char* uplo,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int uplosz,
#endif
                                const char* trans,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int transsz,
#endif
                                const int* n,
                                const int* k,
                                const std::complex<double>* alpha,
                                const std::complex<double>* a,
                                const int* lda,
                                const std::complex<double>* b,
                                const int* ldb,
                                const std::complex<double>* beta,
                                      std::complex<double>* c,
                                const int* ldc);

void  CVM_STD_CALL CHER2K      (const char* uplo,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int uplosz,
#endif
                                const char* trans,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int transsz,
#endif
                                const int* n,
                                const int* k,
                                const std::complex<float>* alpha,
                                const std::complex<float>* a,
                                const int* lda,
                                const std::complex<float>* b,
                                const int* ldb,
                                const std::complex<float>* beta,
                                      std::complex<float>* c,
                                const int* ldc);

void  CVM_STD_CALL ZHER2K      (const char* uplo,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int uplosz,
#endif
                                const char* trans,
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                                const unsigned int transsz,
#endif
                                const int* n,
                                const int* k,
                                const std::complex<double>* alpha,
                                const std::complex<double>* a,
                                const int* lda,
                                const std::complex<double>* b,
                                const int* ldb,
                                const std::complex<double>* beta,
                                      std::complex<double>* c,
                                const int* ldc);

// QR subroutines
void CVM_FTN_CALL SGEQRF       (const int* m,
                                const int* n,
                                      float* a,
                                const int* lda,
                                      float* tau,
                                      float* work,
                                const int* lwork,
                                      int* info);

void CVM_FTN_CALL DGEQRF       (const int* m,
                                const int* n,
                                      double* a,
                                const int* lda,
                                      double* tau,
                                      double* work,
                                const int* lwork,
                                      int* info);

void CVM_FTN_CALL CGEQRF       (const int* m,
                                const int* n,
                                      std::complex<float>* a,
                                const int* lda,
                                      std::complex<float>* tau,
                                      std::complex<float>* work,
                                const int* lwork,
                                      int* info);

void CVM_FTN_CALL ZGEQRF       (const int* m,
                                const int* n,
                                      std::complex<double>* a,
                                const int* lda,
                                      std::complex<double>* tau,
                                      std::complex<double>* work,
                                const int* lwork,
                                      int* info);

void CVM_FTN_CALL SORGQR       (const int* m,
                                const int* n,
                                const int* k,
                                      float* a,
                                const int* lda,
                                const float* tau,
                                      float* work,
                                const int* lwork,
                                      int* info);

void CVM_FTN_CALL DORGQR       (const int* m,
                                const int* n,
                                const int* k,
                                      double* a,
                                const int* lda,
                                const double* tau,
                                      double* work,
                                const int* lwork,
                                      int* info);

void CVM_FTN_CALL CUNGQR       (const int* m,
                                const int* n,
                                const int* k,
                                      std::complex<float>* a,
                                const int* lda,
                                const std::complex<float>* tau,
                                      std::complex<float>* work,
                                const int* lwork,
                                      int* info);

void CVM_FTN_CALL ZUNGQR       (const int* m,
                                const int* n,
                                const int* k,
                                      std::complex<double>* a,
                                const int* lda,
                                const std::complex<double>* tau,
                                      std::complex<double>* work,
                                const int* lwork,
                                      int* info);

#ifdef __cplusplus
}                       // extern "C" 
#endif

#endif                  // !_BLAS_H
