C
C        FORTRAN 77 SUBROUTINE
C        Part of the 3D Projection Pursuit Software Suite
C
C        This program is COPYRIGHT (C) 1990 G P Nason.
C
C
C        AUTHOR:        G P Nason
C
C        TITLE:        CDP3DX
C
C        CREATED:        16/3/90
C        LAST MODIFIED:        04/4/90
C
C        Checked gpn        Visual check compare to notes only.
C        
C
C        Description
C        -----------
C
C        Computes the derivative of the projection index with respect to
C        each of the projection directions.
C
C
C        Variable Description
C        --------------------
C
C        DP3DX        -        Array containing derivatives of the projection
C                        index with respect to the 3k projection direction
C                        components.
C                        The array is REAL and of dimension 3xKMAX. The first
C                        array index refers to which projection vector you are
C                        interested in. So 1=A, 2=B and 3=C. The second array
C                        index refers to the actual component of that vector.
C
C                        So the derivative of the projection index with respect
C                        to the third component of vector B (i.e. B3) is stored
C                        in DP3DX(2,3).
C
C       KMAX    -       One of the FORTRAN dimensions of arrays DP3DX.
C
C       K       -       The dimensionality of the data set.
C
C        KSTAT        -        Array of necessary k-statistics. This is a three
C                        dimensional array, each dimension ranging from
C                        0:SSIZE.
C
C        DKDX        -        Differentials of all necessary k-statistics with
C                        respect to the 3k projection direction components.
C                        This is a five dimensional array. The first three
C                        dimensions range from 0:SSIZE (0:4) and refer to the
C                        k-statistic in question. So DKDX(1,2,0,x,y) refer to
C                        the differential of K(1,2,0) with respect to something.
C
C                        The last two dimensions refer to what you're
C                        differentiating with respect to in a manner similar
C                        to the two dimensions of DP3DX. The last two dimensions
C                        of DKDX are 3xKMAX. The fourth dimension refers to
C                        which projection vector that you're interested in
C                        (again A=1,B=2, and C=3). The fifth dimension refers to
C                        the actual component of that vector.
C
C                        For example the differential  of K(1,1,2) with respect
C                        to C7 (the 7th component of vector C) should be
C                        provided in DKDX(1,1,2,3,7).
C
C       SSIZE   -       Recommended to be 4. It is the FORTRAN dimension of
C                        some of the dimensions of the array DKDX.
C
C       n.b. K must be <= KMAX.
C       Note: The array DKDX has some indices from ZERO to SSIZE.
C
        SUBROUTINE CDP3DX(DP3DX, KMAX, K, KSTAT, DKDX, SSIZE)

        INTEGER KMAX,K,SSIZE
        DOUBLEPRECISION DP3DX(3,KMAX)
        DOUBLEPRECISION KSTAT(0:SSIZE,0:SSIZE,0:SSIZE)
        DOUBLEPRECISION DKDX(0:SSIZE,0:SSIZE,0:SSIZE,3,KMAX)

        INTEGER VECTOR,COMPNT
C
C        To split up a long formula. We split the derivative calculation
C        into the sum of derivatives of third order and fourth order k-statistics
C
        DOUBLEPRECISION THIRD
        DOUBLEPRECISION FOURTH

C
C        We must calculate the derivative for each vector (VECTOR)(A,B,C) and
C        for each component (COMPNT) of each vector.

        DO 20 VECTOR=1,3
            DO 10 COMPNT=1,K

C
C                Using third order k-statistics.
C
                THIRD =     KSTAT(3,0,0)*DKDX(3,0,0,VECTOR,COMPNT) +
     1                  3.0*KSTAT(2,1,0)*DKDX(2,1,0,VECTOR,COMPNT) +
     2                  3.0*KSTAT(2,0,1)*DKDX(2,0,1,VECTOR,COMPNT) +
     3                  3.0*KSTAT(1,2,0)*DKDX(1,2,0,VECTOR,COMPNT) +
     4                  6.0*KSTAT(1,1,1)*DKDX(1,1,1,VECTOR,COMPNT) +
     5                  3.0*KSTAT(1,0,2)*DKDX(1,0,2,VECTOR,COMPNT) +
     6                      KSTAT(0,3,0)*DKDX(0,3,0,VECTOR,COMPNT) +
     7                  3.0*KSTAT(0,2,1)*DKDX(0,2,1,VECTOR,COMPNT) +
     8                  3.0*KSTAT(0,1,2)*DKDX(0,1,2,VECTOR,COMPNT) +
     9                      KSTAT(0,0,3)*DKDX(0,0,3,VECTOR,COMPNT) 

C
C                Using fourth order k-statistics.
C
                FOURTH =     KSTAT(4,0,0)*DKDX(4,0,0,VECTOR,COMPNT) +
     1                   4.0*KSTAT(3,1,0)*DKDX(3,1,0,VECTOR,COMPNT) +
     2                   4.0*KSTAT(3,0,1)*DKDX(3,0,1,VECTOR,COMPNT) +
     3                   6.0*KSTAT(2,2,0)*DKDX(2,2,0,VECTOR,COMPNT) +
     4                  12.0*KSTAT(2,1,1)*DKDX(2,1,1,VECTOR,COMPNT) +
     5                   6.0*KSTAT(2,0,2)*DKDX(2,0,2,VECTOR,COMPNT) +
     6                   4.0*KSTAT(1,3,0)*DKDX(1,3,0,VECTOR,COMPNT) +
     7                  12.0*KSTAT(1,2,1)*DKDX(1,2,1,VECTOR,COMPNT) +
     8                  12.0*KSTAT(1,1,2)*DKDX(1,1,2,VECTOR,COMPNT) +
     9                   4.0*KSTAT(1,0,3)*DKDX(1,0,3,VECTOR,COMPNT) +
     +                       KSTAT(0,4,0)*DKDX(0,4,0,VECTOR,COMPNT) +
     1                   4.0*KSTAT(0,3,1)*DKDX(0,3,1,VECTOR,COMPNT) +
     2                   6.0*KSTAT(0,2,2)*DKDX(0,2,2,VECTOR,COMPNT) +
     3                   4.0*KSTAT(0,1,3)*DKDX(0,1,3,VECTOR,COMPNT) +
     4                       KSTAT(0,0,4)*DKDX(0,0,4,VECTOR,COMPNT)
C123456789012345678901234567890123456789012345678901234567890123456789012

C
C                Now put them together
C
                DP3DX(VECTOR,COMPNT) = 2.0*THIRD + 0.5*FOURTH
 10            CONTINUE
 20        CONTINUE
        RETURN
        END
