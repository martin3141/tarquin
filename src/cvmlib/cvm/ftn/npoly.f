C     Working array size calculator
C
C     Copyright (C), Sergei Nikolaev, 1992-2006, http://cvmlib.com
C

      INTEGER*4 FUNCTION NPOLY (M, N)
CDEC$ IF DEFINED (FTN_EXPORTS)
CDEC$     ATTRIBUTES DLLEXPORT::NPOLY
CDEC$ ENDIF
      INTEGER*4 M, N
      INTEGER*4 FLOOR, CEILING

      NPOLY = 0
      IF (M .GT. 0 .AND. N .GT. 1) THEN
          NPOLY  = (FLOOR (DFLOAT (N - 1) /
     1              DFLOAT (CEILING (DSQRT (DFLOAT (N - 1)))))
     2                                               + 2) * M * M
      ENDIF
      RETURN
      END !FUNCTION NPOLY
