C     Matrix axpy routines
C
C     Copyright (C), Sergei Nikolaev, 1992-2006, http://cvmlib.com
C
C     The ?axpy routines perform a vector-vector operation defined as
C     y := a*x + y
C
C
C     Input/Output parameters:
C
C     M   - rows (int)(input)
C     N   - columns (int)(input)
C     A   - multiplier (real)(input)
C     X   - source matrix (real)(input)
C     LDX - leading dimesion of X (int)(input)
C     Y   - destination matrix (real)(output)
C     LDY - leading dimesion of Y (int)(input)

      SUBROUTINE SAXPYM (M, N, A, X, LDX, Y, LDY)
CDEC$ IF DEFINED (FTN_EXPORTS)
CDEC$     ATTRIBUTES DLLEXPORT::SAXPYM
CDEC$ ENDIF
      INTEGER*4 M, N, LDX, LDY
      REAL*4    A, X(LDX*N), Y(LDY*N)
      INTEGER*4 I

      IF (M .EQ. LDX .AND. M .EQ. LDY) THEN
          CALL SAXPY (M * N, A, X, 1, Y, 1)
      ELSE
          DO 10 I = 0, N-1
              CALL SAXPY (M, A, X(I*LDX+1), 1, Y(I*LDY+1), 1)
10        CONTINUE
      END IF
      RETURN
      END !SUBROUTINE SAXPYM


      SUBROUTINE DAXPYM (M, N, A, X, LDX, Y, LDY)
CDEC$ IF DEFINED (FTN_EXPORTS)
CDEC$     ATTRIBUTES DLLEXPORT::DAXPYM
CDEC$ ENDIF
      INTEGER*4 M, N, LDX, LDY
      REAL*8    A, X(LDX*N), Y(LDY*N)
      INTEGER*4 I

      IF (M .EQ. LDX .AND. M .EQ. LDY) THEN
          CALL DAXPY (M * N, A, X, 1, Y, 1)
      ELSE
          DO 20 I = 0, N-1
              CALL DAXPY (M, A, X(I*LDX+1), 1, Y(I*LDY+1), 1)
20        CONTINUE
      END IF
      RETURN
      END !SUBROUTINE DAXPYM


      SUBROUTINE CAXPYM (M, N, A, X, LDX, Y, LDY)
CDEC$ IF DEFINED (FTN_EXPORTS)
CDEC$     ATTRIBUTES DLLEXPORT::CAXPYM
CDEC$ ENDIF
      INTEGER*4 M, N, LDX, LDY
      COMPLEX*8 A, X(LDX*N), Y(LDY*N)
      INTEGER*4 I

      IF (M .EQ. LDX .AND. M .EQ. LDY) THEN
          CALL CAXPY (M * N, A, X, 1, Y, 1)
      ELSE
          DO 30 I = 0, N-1
              CALL CAXPY (M, A, X(I*LDX+1), 1, Y(I*LDY+1), 1)
30        CONTINUE
      END IF
      RETURN
      END !SUBROUTINE CAXPYM


      SUBROUTINE ZAXPYM (M, N, A, X, LDX, Y, LDY)
CDEC$ IF DEFINED (FTN_EXPORTS)
CDEC$     ATTRIBUTES DLLEXPORT::ZAXPYM
CDEC$ ENDIF
      INTEGER*4  M, N, LDX, LDY
      COMPLEX*16 A, X(LDX*N), Y(LDY*N)
      INTEGER*4  I

      IF (M .EQ. LDX .AND. M .EQ. LDY) THEN
          CALL ZAXPY (M * N, A, X, 1, Y, 1)
      ELSE
          DO 40 I = 0, N-1
              CALL ZAXPY (M, A, X(I*LDX+1), 1, Y(I*LDY+1), 1)
40        CONTINUE
      END IF
      RETURN
      END !SUBROUTINE ZAXPYM

