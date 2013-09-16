C     Matrix copy routines
C
C     Copyright (C), Sergei Nikolaev, 1992-2006, http://cvmlib.com
C
C     Input/Output parameters:
C
C     M   - rows (int)(input)
C     N   - columns (int)(input)
C     A   - source matrix (real)(input)
C     B   - destination matrix (real)(output)
C     LDA - leading dimesion of A (int)(input)
C     LDB - leading dimesion of B (int)(input)

      SUBROUTINE SCOPYM (M, N, A, LDA, B, LDB) 
CDEC$ IF DEFINED (FTN_EXPORTS)
CDEC$     ATTRIBUTES DLLEXPORT::SCOPYM
CDEC$ ENDIF
      INTEGER*4 M, N, LDA, LDB
      REAL*4    A(LDA*N), B(LDB*N)
      INTEGER*4 I

      IF (M .EQ. LDA .AND. M .EQ. LDB) THEN
          CALL SCOPY (M * N, A, 1, B, 1)
      ELSE
          DO 10 I = 0, N-1
              CALL SCOPY (M, A(I*LDA+1), 1, B(I*LDB+1), 1)
10        CONTINUE
      END IF
      RETURN
      END !SUBROUTINE SCOPYM

      SUBROUTINE DCOPYM (M, N, A, LDA, B, LDB) 
CDEC$ IF DEFINED (FTN_EXPORTS)
CDEC$     ATTRIBUTES DLLEXPORT::DCOPYM
CDEC$ ENDIF
      INTEGER*4 M, N, LDA, LDB
      REAL*8    A(LDA*N), B(LDB*N)
      INTEGER*4 I

      IF (M .EQ. LDA .AND. M .EQ. LDB) THEN
          CALL DCOPY (M * N, A, 1, B, 1)
      ELSE
          DO 20 I = 0, N-1
              CALL DCOPY (M, A(I*LDA+1), 1, B(I*LDB+1), 1)
20        CONTINUE
      END IF
      RETURN
      END !SUBROUTINE DCOPYM

      SUBROUTINE CCOPYM (M, N, A, LDA, B, LDB) 
CDEC$ IF DEFINED (FTN_EXPORTS)
CDEC$     ATTRIBUTES DLLEXPORT::CCOPYM
CDEC$ ENDIF
      INTEGER*4 M, N, LDA, LDB
      COMPLEX*8 A(LDA*N), B(LDB*N)
      INTEGER*4 I

      IF (M .EQ. LDA .AND. M .EQ. LDB) THEN
          CALL CCOPY (M * N, A, 1, B, 1)
      ELSE
          DO 30 I = 0, N-1
              CALL CCOPY (M, A(I*LDA+1), 1, B(I*LDB+1), 1)
30        CONTINUE
      END IF
      RETURN
      END !SUBROUTINE CCOPYM

      SUBROUTINE ZCOPYM (M, N, A, LDA, B, LDB) 
CDEC$ IF DEFINED (FTN_EXPORTS)
CDEC$     ATTRIBUTES DLLEXPORT::ZCOPYM
CDEC$ ENDIF
      INTEGER*4  M, N, LDA, LDB
      COMPLEX*16 A(LDA*N), B(LDB*N)
      INTEGER*4  I

      IF (M .EQ. LDA .AND. M .EQ. LDB) THEN
          CALL ZCOPY (M * N, A, 1, B, 1)
      ELSE
          DO 40 I = 0, N-1
              CALL ZCOPY (M, A(I*LDA+1), 1, B(I*LDB+1), 1)
40        CONTINUE
      END IF
      RETURN
      END !SUBROUTINE ZCOPYM

