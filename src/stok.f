C
C        FORTRAN 77 SUBROUTINE
C        Part of the 3D Projection Pursuit Software Suite
C
C        This program is COPYRIGHT (C) 1990 G P Nason 
C
C
C        AUTHOR:        G P Nason
C
C        TITLE:        STOK
C
C        CREATED:        18/3/90
C        LAST MODIFIED:        04/4/90
C
C        040490 gpn        Added the constant n[p] for all the kstatistics.
C        
C
C        Description
C        -----------
C
C        Calculates the trivariate k-statistics KSTAT(,,) from the power
C        sums S(,,).
C
C
C        Variable Description
C        --------------------
C
C        KSTAT        -        An array that contains the necessary trivariate
C                        k-statistics on exit. It is a three-dimensional
C                        REAL array, with each dimension ranging from
C                        0:SSIZE.
C
C        SSIZE        -        Actual FORTRAN dimension of KSTAT. Recommended value=4
C
C        S        -        Supplied trivariate power sums (calculated by CALCS).
C                        This is also a 3-dimensional REAL array with each
C                        dimension ranging 0:SSIZE.
C
C        N        -        The number of data-points.
C
        SUBROUTINE STOK(KSTAT, SSIZE, S, N)

        INTEGER SSIZE
        DOUBLEPRECISION KSTAT(0:SSIZE,0:SSIZE,0:SSIZE)
        DOUBLEPRECISION S(0:SSIZE,0:SSIZE,0:SSIZE)
        INTEGER N
        DOUBLEPRECISION CONS3
        DOUBLEPRECISION CONS4

C
C        Multiplication constant for third order sums
C
        CONS3 = DBLE(N)/DBLE((N-1)*(N-2))
C
C        Multiplication constant for fourth order sums
C
        CONS4 = 1.0/DBLE((N-1)*(N-2)*(N-3))

C
C        Calculate k-statistics from power sums
C
C
C        Third order k-statistics
C
        KSTAT(1,1,1) = CONS3*S(1,1,1)
        KSTAT(2,1,0) = CONS3*S(2,1,0)
        KSTAT(2,0,1) = CONS3*S(2,0,1)
        KSTAT(1,2,0) = CONS3*S(1,2,0)
        KSTAT(1,0,2) = CONS3*S(1,0,2)
        KSTAT(0,2,1) = CONS3*S(0,2,1)
        KSTAT(0,1,2) = CONS3*S(0,1,2)
        KSTAT(3,0,0) = CONS3*S(3,0,0)
        KSTAT(0,3,0) = CONS3*S(0,3,0)
        KSTAT(0,0,3) = CONS3*S(0,0,3)
C
C        Fourth order k-statistics
C
        KSTAT(2,1,1) = CONS4*N*(N+1)*S(2,1,1)
        KSTAT(1,2,1) = CONS4*N*(N+1)*S(1,2,1)
        KSTAT(1,1,2) = CONS4*N*(N+1)*S(1,1,2)
        KSTAT(3,1,0) = CONS4*N*(N+1)*S(3,1,0)
        KSTAT(3,0,1) = CONS4*N*(N+1)*S(3,0,1)
        KSTAT(1,3,0) = CONS4*N*(N+1)*S(1,3,0)
        KSTAT(1,0,3) = CONS4*N*(N+1)*S(1,0,3)
        KSTAT(0,3,1) = CONS4*N*(N+1)*S(0,3,1)
        KSTAT(0,1,3) = CONS4*N*(N+1)*S(0,1,3)
C
        KSTAT(4,0,0) = CONS4*(N*(N+1)*S(4,0,0) + 3.0*(1-N))
        KSTAT(0,4,0) = CONS4*(N*(N+1)*S(0,4,0) + 3.0*(1-N))
        KSTAT(0,0,4) = CONS4*(N*(N+1)*S(0,0,4) + 3.0*(1-N))
C
        KSTAT(0,2,2) = CONS4*(N*(N+1)*S(0,2,2) + (1-N))
        KSTAT(2,0,2) = CONS4*(N*(N+1)*S(2,0,2) + (1-N))
        KSTAT(2,2,0) = CONS4*(N*(N+1)*S(2,2,0) + (1-N))

        RETURN
        END
