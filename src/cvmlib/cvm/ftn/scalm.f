C     Matrix scaling routines
C
C     Copyright (C), Sergei Nikolaev, 1992-2006, http://cvmlib.com
C
C     Input/Output parameters:
C
C     M   - rows (int)(input)
C     N   - columns (int)(input)
C     S   - scale factor (real)(input)
C     A   - matrix to be scaled (real)(input, output)
C     LDA - leading dimesion of A (int)(input)

      SUBROUTINE SSCALM (M, N, S, A, LDA) 
CDEC$ IF DEFINED (FTN_EXPORTS)
CDEC$     ATTRIBUTES DLLEXPORT::SSCALM
CDEC$ ENDIF
      INTEGER*4 M, N, LDA
      REAL*4    S
      REAL*4    A(LDA*N)
      INTEGER*4 I

      IF (M .EQ. LDA) THEN
          CALL SSCAL (M * N, S, A, 1)
      ELSE
          DO 10 I = 0, N-1
              CALL SSCAL (M, S, A(I*LDA+1), 1)
10        CONTINUE
      END IF
      RETURN
      END !SUBROUTINE SSCALM

      SUBROUTINE DSCALM (M, N, S, A, LDA) 
CDEC$ IF DEFINED (FTN_EXPORTS)
CDEC$     ATTRIBUTES DLLEXPORT::DSCALM
CDEC$ ENDIF
      INTEGER*4 M, N, LDA
      REAL*8    S
      REAL*8    A(LDA*N)
      INTEGER*4 I

      IF (M .EQ. LDA) THEN
          CALL DSCAL (M * N, S, A, 1)
      ELSE
          DO 10 I = 0, N-1
              CALL DSCAL (M, S, A(I*LDA+1), 1)
10        CONTINUE
      END IF
      RETURN
      END !SUBROUTINE DSCALM

      SUBROUTINE CSCALM (M, N, S, A, LDA) 
CDEC$ IF DEFINED (FTN_EXPORTS)
CDEC$     ATTRIBUTES DLLEXPORT::CSCALM
CDEC$ ENDIF
      INTEGER*4 M, N, LDA
      COMPLEX*8 S
      COMPLEX*8 A(LDA*N)
      INTEGER*4 I

      IF (M .EQ. LDA) THEN
          CALL CSCAL (M * N, S, A, 1)
      ELSE
          DO 10 I = 0, N-1
              CALL CSCAL (M, S, A(I*LDA+1), 1)
10        CONTINUE
      END IF
      RETURN
      END !SUBROUTINE CSCALM

      SUBROUTINE ZSCALM (M, N, S, A, LDA) 
CDEC$ IF DEFINED (FTN_EXPORTS)
CDEC$     ATTRIBUTES DLLEXPORT::ZSCALM
CDEC$ ENDIF
      INTEGER*4  M, N, LDA
      COMPLEX*16 S
      COMPLEX*16 A(LDA*N)
      INTEGER*4  I

      IF (M .EQ. LDA) THEN
          CALL ZSCAL (M * N, S, A, 1)
      ELSE
          DO 10 I = 0, N-1
              CALL ZSCAL (M, S, A(I*LDA+1), 1)
10        CONTINUE
      END IF
      RETURN
      END !SUBROUTINE ZSCALM
