C     Matrix infinity norm routines
C
C     Copyright (C), Sergei Nikolaev, 1992-2006, http://cvmlib.com
C
C     Input/Output parameters:
C
C     M   - rows (int)(input)
C     N   - columns (int)(input)
C     A   - matrix (real)(input)
C     LDA - leading dimesion of A (int)(input)

      REAL*4 FUNCTION SINFNM (M, N, A, LDA) 
CDEC$ IF DEFINED (FTN_EXPORTS)
CDEC$     ATTRIBUTES DLLEXPORT::SINFNM
CDEC$ ENDIF
      INTEGER*4 M, N, LDA
      REAL*4    A(LDA*N)
      INTEGER*4 I
      REAL*4    S
      INTEGER*4 ISAMAX

      SINFNM = 0.
      IF (M .EQ. LDA) THEN
          SINFNM = ABS (A (ISAMAX (M * N, A, 1)))
      ELSE
          DO 10 I = 0, N-1
              S = ABS (A (ISAMAX (M, A(I*LDA+1), 1)))
              IF (S .GT. SINFNM) SINFNM = S
10        CONTINUE
      END IF
      RETURN
      END !FUNCTION SINFNM


      REAL*8 FUNCTION DINFNM (M, N, A, LDA) 
CDEC$ IF DEFINED (FTN_EXPORTS)
CDEC$     ATTRIBUTES DLLEXPORT::DINFNM
CDEC$ ENDIF
      INTEGER*4 M, N, LDA
      REAL*8    A(LDA*N)
      INTEGER*4 I
      REAL*8    S
      INTEGER*4 IDAMAX

      DINFNM = 0.D0
      IF (M .EQ. LDA) THEN
          DINFNM = DABS (A (IDAMAX (M * N, A, 1)))
      ELSE
          DO 20 I = 0, N-1
              S = DABS (A (IDAMAX (M, A(I*LDA+1), 1)))
              IF (S .GT. DINFNM) DINFNM = S
20        CONTINUE
      END IF
      RETURN
      END !FUNCTION DINFNM


      REAL*4 FUNCTION CINFNM (M, N, A, LDA) 
CDEC$ IF DEFINED (FTN_EXPORTS)
CDEC$     ATTRIBUTES DLLEXPORT::CINFNM
CDEC$ ENDIF
      INTEGER*4 M, N, LDA
      COMPLEX*8 A(LDA*N)
      INTEGER*4 I
      REAL*4    S
      INTEGER*4 ICAMAX

      CINFNM = (0., 0.)
      IF (M .EQ. LDA) THEN
          CINFNM = CABS (A (ICAMAX (M * N, A, 1)))
      ELSE
          DO 20 I = 0, N-1
              S = CABS (A (ICAMAX (M, A(I*LDA+1), 1)))
              IF (S .GT. CINFNM) CINFNM = S
20        CONTINUE
      END IF
      RETURN
      END !FUNCTION CINFNM


      REAL*8 FUNCTION ZINFNM (M, N, A, LDA) 
CDEC$ IF DEFINED (FTN_EXPORTS)
CDEC$     ATTRIBUTES DLLEXPORT::ZINFNM
CDEC$ ENDIF
      INTEGER*4  M, N, LDA
      COMPLEX*16 A(LDA*N)
      INTEGER*4  I
      REAL*8     S
      INTEGER*4  IZAMAX

      ZINFNM = (0.D0, 0.D0)
      IF (M .EQ. LDA) THEN
          ZINFNM = ZABS (A (IZAMAX (M * N, A, 1)))
      ELSE
          DO 20 I = 0, N-1
              S = ZABS (A (IZAMAX (M, A(I*LDA+1), 1)))
              IF (S .GT. ZINFNM) ZINFNM = S
20        CONTINUE
      END IF
      RETURN
      END !FUNCTION ZINFNM
