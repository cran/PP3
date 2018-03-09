C
C        FORTRAN 77 SUBROUTINE
C        Part of the 3D Projection Pursuit Software Suite
C
C        This program is COPYRIGHT (C) 1990 G P Nason 
C
C
C        AUTHOR:        G P Nason
C
C        TITLE:        P3
C
C        CREATED:        18/3/90
C        LAST MODIFIED:        18/3/90
C        
C
C        Description
C        -----------
C
C        Calculates the 3-dimensional moment projection index from the
C        trivariate k-statistics KSTAT(,,). (Calculated from STOK).
C
C
C        Variable Description
C        --------------------
C
C        P3INDX        -        The projection index.
C
C        KSTAT        -        Array of precalculated k-statistics. This is
C                        a 3-dimensional array. Each dimension ranging from
C                        0:SSIZE
C
C        SSIZE        -        FORTRAN dimension of KSTAT. Recommended to be 4.
C
        SUBROUTINE P3(P3INDX, KSTAT, SSIZE)
        
        INTEGER SSIZE
        DOUBLEPRECISION P3INDX
        DOUBLEPRECISION KSTAT(0:SSIZE, 0:SSIZE, 0:SSIZE)

        DOUBLEPRECISION THIRD,FOURTH

C
C        Calculate contribution from third order k-statistics.
C
        THIRD = KSTAT(3,0,0)**2 + 3.0*KSTAT(2,1,0)**2 +
     1                3.0*KSTAT(2,0,1)**2 + 3.0*KSTAT(1,2,0)**2 +
     2                6.0*KSTAT(1,1,1)**2 + 3.0*KSTAT(1,0,2)**2 +
     3                KSTAT(0,3,0)**2 + 3.0*KSTAT(0,2,1)**2 +
     4                3.0*KSTAT(0,1,2)**2 + KSTAT(0,0,3)**2
C
C        Calculate contribution from fourth order k-statistics.
C
        FOURTH = KSTAT(4,0,0)**2 + 4.0*KSTAT(3,1,0)**2 +
     1                4.0*KSTAT(3,0,1)**2 + 6.0*KSTAT(2,2,0)**2 +
     2                12.0*KSTAT(2,1,1)**2 + 6.0*KSTAT(2,0,2)**2 +
     3                4.0*KSTAT(1,3,0)**2 + 12.0*KSTAT(1,2,1)**2 +
     4                12.0*KSTAT(1,1,2)**2 + 4.0*KSTAT(1,0,3)**2 +
     5                KSTAT(0,4,0)**2 + 4.0*KSTAT(0,3,1)**2 +
     6                6.0*KSTAT(0,2,2)**2 + 4.0*KSTAT(0,1,3)**2 +
     7                KSTAT(0,0,4)**2
C
C        Put them together
C
        P3INDX = THIRD + 0.25*FOURTH
        RETURN
        END
