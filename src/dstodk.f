C
C        FORTRAN 77 SUBROUTINE
C        Part of the 3D Projection Pursuit Software Suite
C
C        This program is COPYRIGHT (C) 1990 G P Nason 
C
C
C        AUTHOR:        G P Nason
C
C        TITLE:        DSTODK
C
C        CREATED:        18/3/90
C        LAST MODIFIED:        04/4/90
C
C        040490        gpn        Added n[p] multiplication factor
C        
C
C        Description
C        -----------
C
C        Calculates the derivatives of trivariate k-statistics KSTAT(,,) from
C        the derivatives of power sums S(,,).
C
C
C        Variable Description
C        --------------------
C
C        DKDX        -        Array that will contain the derivatives of the
C                        trivariate k-statistics with
C                       respect to the 3k projection direction components.
C                       This is a five dimensional array. The first three
C                       dimensions range from 0:SSIZE (0:4) and refer to the
C                       k-statistic in question. So DKDX(1,2,0,x,y) refer to
C                       the differential of K(1,2,0) with respect to something.
C
C                       The last two dimensions refer to what you're
C                       differentiating with respect to in a manner similar
C                       to the two dimensions of DP3DX. The last two dimensions
C                       of DKDX are 3xKMAX. The fourth dimension refers to
C                       which projection vector that you're interested in
C                       (again A=1,B=2, and C=3). The fifth dimension refers to
C                       the actual component of that vector.
C
C                       For example the differential  of K(1,1,2) with respect
C                       to C7 (the 7th component of vector C) should be
C                       provided in DKDX(1,1,2,3,7).
C
C        SSIZE        -        Actual FORTRAN dimension of KSTAT. Recommended value=4
C
C       DSA     -       Derivatives of S with respect to components of A.
C                       This is a four dimensional array. The first three
C                       dimensions index the element of S we are differentiating
C                       and the last dimension refers to which element of A
C                       we are differentiating by. So the dimension of this
C                       array has to be SSIZE*SSIZE*SSIZE*KMAX.
C
C       DSB     -       As for DSA except for differentiating with respect to
C                       components of B.
C
C       DSC     -       As for DSA except for differentiating with respect to
C                       components of C.
C
C        N        -        The number of data-points.
C
C        KMAX        -        Fourth dimension size of DSA,DSB & DSC. Fifth
C                        dimension size of DKDX.
C
C        K        -        Dimensionality of data set.
C
        SUBROUTINE DSTODK(DKDX, SSIZE, DSA, DSB, DSC, N, KMAX, K)

        INTEGER SSIZE,KMAX,K
        DOUBLEPRECISION DKDX(0:SSIZE,0:SSIZE,0:SSIZE,3,KMAX)
        DOUBLEPRECISION DSA(0:SSIZE,0:SSIZE,0:SSIZE,KMAX)
        DOUBLEPRECISION DSB(0:SSIZE,0:SSIZE,0:SSIZE,KMAX)
        DOUBLEPRECISION DSC(0:SSIZE,0:SSIZE,0:SSIZE,KMAX)
        INTEGER N

        INTEGER VECTOR
        INTEGER COMPNT

        DOUBLEPRECISION CONS3
        DOUBLEPRECISION CONS4

C
C        Calculate 3rd order multiplication factor
C
        CONS3 = DBLE(N)/DBLE((N-1)*(N-2))
C
C        Calculate 4th order multiplication factor
C
        CONS4 = DBLE(N*(N+1))/DBLE((N-1)*(N-2)*(N-3))

C
C        Calculate derivatives of k-statistics from derivs of power sums
C

C
C        For derivatives wrt components of A
C
        VECTOR=1

        DO 10 COMPNT=1,K
C
C        Third order k-statistics
C
        DKDX(1,1,1,VECTOR,COMPNT) = CONS3*DSA(1,1,1,COMPNT)
        DKDX(2,1,0,VECTOR,COMPNT) = CONS3*DSA(2,1,0,COMPNT)
        DKDX(2,0,1,VECTOR,COMPNT) = CONS3*DSA(2,0,1,COMPNT)
        DKDX(1,2,0,VECTOR,COMPNT) = CONS3*DSA(1,2,0,COMPNT)
        DKDX(1,0,2,VECTOR,COMPNT) = CONS3*DSA(1,0,2,COMPNT)
        DKDX(0,2,1,VECTOR,COMPNT) = CONS3*DSA(0,2,1,COMPNT)
        DKDX(0,1,2,VECTOR,COMPNT) = CONS3*DSA(0,1,2,COMPNT)
        DKDX(3,0,0,VECTOR,COMPNT) = CONS3*DSA(3,0,0,COMPNT)
        DKDX(0,3,0,VECTOR,COMPNT) = CONS3*DSA(0,3,0,COMPNT)
        DKDX(0,0,3,VECTOR,COMPNT) = CONS3*DSA(0,0,3,COMPNT)
C
C        Fourth order k-statistics
C
        DKDX(2,1,1,VECTOR,COMPNT) = CONS4*DSA(2,1,1,COMPNT)
        DKDX(1,2,1,VECTOR,COMPNT) = CONS4*DSA(1,2,1,COMPNT)
        DKDX(1,1,2,VECTOR,COMPNT) = CONS4*DSA(1,1,2,COMPNT)
        DKDX(3,1,0,VECTOR,COMPNT) = CONS4*DSA(3,1,0,COMPNT)
        DKDX(3,0,1,VECTOR,COMPNT) = CONS4*DSA(3,0,1,COMPNT)
        DKDX(1,3,0,VECTOR,COMPNT) = CONS4*DSA(1,3,0,COMPNT)
        DKDX(1,0,3,VECTOR,COMPNT) = CONS4*DSA(1,0,3,COMPNT)
        DKDX(0,3,1,VECTOR,COMPNT) = CONS4*DSA(0,3,1,COMPNT)
        DKDX(0,1,3,VECTOR,COMPNT) = CONS4*DSA(0,1,3,COMPNT)
C
        DKDX(4,0,0,VECTOR,COMPNT) = CONS4*DSA(4,0,0,COMPNT)
        DKDX(0,4,0,VECTOR,COMPNT) = CONS4*DSA(0,4,0,COMPNT)
        DKDX(0,0,4,VECTOR,COMPNT) = CONS4*DSA(0,0,4,COMPNT)
C
        DKDX(0,2,2,VECTOR,COMPNT) = CONS4*DSA(0,2,2,COMPNT)
        DKDX(2,0,2,VECTOR,COMPNT) = CONS4*DSA(2,0,2,COMPNT)
        DKDX(2,2,0,VECTOR,COMPNT) = CONS4*DSA(2,2,0,COMPNT)

 10        CONTINUE

C
C        For derivatives wrt components of B
C
        VECTOR=2

        DO 20 COMPNT=1,K
C
C        Third order k-statistics
C
        DKDX(1,1,1,VECTOR,COMPNT) = CONS3*DSB(1,1,1,COMPNT)
        DKDX(2,1,0,VECTOR,COMPNT) = CONS3*DSB(2,1,0,COMPNT)
        DKDX(2,0,1,VECTOR,COMPNT) = CONS3*DSB(2,0,1,COMPNT)
        DKDX(1,2,0,VECTOR,COMPNT) = CONS3*DSB(1,2,0,COMPNT)
        DKDX(1,0,2,VECTOR,COMPNT) = CONS3*DSB(1,0,2,COMPNT)
        DKDX(0,2,1,VECTOR,COMPNT) = CONS3*DSB(0,2,1,COMPNT)
        DKDX(0,1,2,VECTOR,COMPNT) = CONS3*DSB(0,1,2,COMPNT)
        DKDX(3,0,0,VECTOR,COMPNT) = CONS3*DSB(3,0,0,COMPNT)
        DKDX(0,3,0,VECTOR,COMPNT) = CONS3*DSB(0,3,0,COMPNT)
        DKDX(0,0,3,VECTOR,COMPNT) = CONS3*DSB(0,0,3,COMPNT)
C
C        Fourth order k-statistics
C
        DKDX(2,1,1,VECTOR,COMPNT) = CONS4*DSB(2,1,1,COMPNT)
        DKDX(1,2,1,VECTOR,COMPNT) = CONS4*DSB(1,2,1,COMPNT)
        DKDX(1,1,2,VECTOR,COMPNT) = CONS4*DSB(1,1,2,COMPNT)
        DKDX(3,1,0,VECTOR,COMPNT) = CONS4*DSB(3,1,0,COMPNT)
        DKDX(3,0,1,VECTOR,COMPNT) = CONS4*DSB(3,0,1,COMPNT)
        DKDX(1,3,0,VECTOR,COMPNT) = CONS4*DSB(1,3,0,COMPNT)
        DKDX(1,0,3,VECTOR,COMPNT) = CONS4*DSB(1,0,3,COMPNT)
        DKDX(0,3,1,VECTOR,COMPNT) = CONS4*DSB(0,3,1,COMPNT)
        DKDX(0,1,3,VECTOR,COMPNT) = CONS4*DSB(0,1,3,COMPNT)
C
        DKDX(4,0,0,VECTOR,COMPNT) = CONS4*DSB(4,0,0,COMPNT)
        DKDX(0,4,0,VECTOR,COMPNT) = CONS4*DSB(0,4,0,COMPNT)
        DKDX(0,0,4,VECTOR,COMPNT) = CONS4*DSB(0,0,4,COMPNT)
C
        DKDX(0,2,2,VECTOR,COMPNT) = CONS4*DSB(0,2,2,COMPNT)
        DKDX(2,0,2,VECTOR,COMPNT) = CONS4*DSB(2,0,2,COMPNT)
        DKDX(2,2,0,VECTOR,COMPNT) = CONS4*DSB(2,2,0,COMPNT)

 20        CONTINUE

C
C        For derivatives wrt components of C
C
        VECTOR=3

        DO 30 COMPNT=1,K
C
C        Third order k-statistics
C
        DKDX(1,1,1,VECTOR,COMPNT) = CONS3*DSC(1,1,1,COMPNT)
        DKDX(2,1,0,VECTOR,COMPNT) = CONS3*DSC(2,1,0,COMPNT)
        DKDX(2,0,1,VECTOR,COMPNT) = CONS3*DSC(2,0,1,COMPNT)
        DKDX(1,2,0,VECTOR,COMPNT) = CONS3*DSC(1,2,0,COMPNT)
        DKDX(1,0,2,VECTOR,COMPNT) = CONS3*DSC(1,0,2,COMPNT)
        DKDX(0,2,1,VECTOR,COMPNT) = CONS3*DSC(0,2,1,COMPNT)
        DKDX(0,1,2,VECTOR,COMPNT) = CONS3*DSC(0,1,2,COMPNT)
        DKDX(3,0,0,VECTOR,COMPNT) = CONS3*DSC(3,0,0,COMPNT)
        DKDX(0,3,0,VECTOR,COMPNT) = CONS3*DSC(0,3,0,COMPNT)
        DKDX(0,0,3,VECTOR,COMPNT) = CONS3*DSC(0,0,3,COMPNT)
C
C        Fourth order k-statistics
C
        DKDX(2,1,1,VECTOR,COMPNT) = CONS4*DSC(2,1,1,COMPNT)
        DKDX(1,2,1,VECTOR,COMPNT) = CONS4*DSC(1,2,1,COMPNT)
        DKDX(1,1,2,VECTOR,COMPNT) = CONS4*DSC(1,1,2,COMPNT)
        DKDX(3,1,0,VECTOR,COMPNT) = CONS4*DSC(3,1,0,COMPNT)
        DKDX(3,0,1,VECTOR,COMPNT) = CONS4*DSC(3,0,1,COMPNT)
        DKDX(1,3,0,VECTOR,COMPNT) = CONS4*DSC(1,3,0,COMPNT)
        DKDX(1,0,3,VECTOR,COMPNT) = CONS4*DSC(1,0,3,COMPNT)
        DKDX(0,3,1,VECTOR,COMPNT) = CONS4*DSC(0,3,1,COMPNT)
        DKDX(0,1,3,VECTOR,COMPNT) = CONS4*DSC(0,1,3,COMPNT)
C
        DKDX(4,0,0,VECTOR,COMPNT) = CONS4*DSC(4,0,0,COMPNT)
        DKDX(0,4,0,VECTOR,COMPNT) = CONS4*DSC(0,4,0,COMPNT)
        DKDX(0,0,4,VECTOR,COMPNT) = CONS4*DSC(0,0,4,COMPNT)
C
        DKDX(0,2,2,VECTOR,COMPNT) = CONS4*DSC(0,2,2,COMPNT)
        DKDX(2,0,2,VECTOR,COMPNT) = CONS4*DSC(2,0,2,COMPNT)
        DKDX(2,2,0,VECTOR,COMPNT) = CONS4*DSC(2,2,0,COMPNT)

 30        CONTINUE

        RETURN
        END
