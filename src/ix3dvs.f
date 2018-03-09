C
C        FORTRAN 77 SUBROUTINE
C        Part of the 3D Projection Pursuit Software Suite
C
C        This program is COPYRIGHT (C) 1990 G P Nason 
C
C
C        AUTHOR:        G P Nason
C
C        TITLE:        IX3DVS
C
C        CREATED:        18/3/90
C        LAST MODIFIED:        25/4/90
C
C        060490        gpn        Change format of printing of orthogonalised vectors
C                        into a REDUCE type output.
C
C        250490        gpn        Add printing of projection index and it's derivatives
C                        if the TEXT variable is 1.
C        
C
C        Description
C        -----------
C
C        Calculates the projection index and all of it's derivatives given
C        a projection plane defined by A,B and C. You must also supply
C        precalculated product moments T(,,) and U(,,,).
C
C
C        Variable Description
C        --------------------
C
C        P3INDX        -        REAL. Returned value of the projection index
C
C        DP3DX        -        REAL array. Dimension 3xKMAX. Returned values of the
C                        derivatives of the projection index with respect to
C                        each of the direction vectors (first dimension) and
C                        each of the components of these vectors (second
C                        dimension).
C
C        KMAX        -        INTEGER. FORTRAN dimension. Puts an upper limit on
C                        the dimensionality of the data set.
C
C        K        -        INTEGER. The dimensionality of the data set.
C
C        N        -        INTEGER. The number of points in the data set.
C
C        KSTAT        -        REAL array of workspace. A 3-dimensional array, each
C                        dimension ranging from 0 to SSIZE. On exit contains the
C                        current values of the projected k-statistics.
C
C        DKDX        -        REAL array of workspace. A 5-dimensional array.
C                        The first 3 dimensions ranging from 0:SSIZE. The fourth
C                        ranges from 1:3 and the fifth from 1:KMAX. On exit
C                        contains the current values of the derivatives of the
C                        projected trivariate k-statistics with respect to
C                        each of the components of the projection directions.
C
C        S        -        REAL array of workspace. A 3-dimensional array,
C                        each dimension ranging from 0:SSIZE. On exit contains
C                        the current values of the projected power-sums.
C
C        DSA,DSB,DSC        REAL arrays of workspace. Each a 4-dimensional array
C                        The first 3 dimensions range from 0:SSIZE and the
C                        fourth 1:KMAX. On exit contains the derivatives of
C                        the projected power sums with respect to the components
C                        of the projection directions.
C
C        SSIZE        -        INTEGER. One of the dimensions of many arrays. We
C                        calculate things like K(0,1,2), S(0,0,4) etc. The most
C                        any of these indices will be is 4, the least is 0. So
C                        this variable MUST be 4 !
C
C        A,B,C        -        REAL arrays. Supply the three projection directions
C                        in these three arrays, Each contains KMAX elements.
C
C        T        -        REAL array. A three-dimensional array, each dimension
C                        ranges from 1:KMAX. This is the third order product
C                        moment tensor and should be supplied by the user
C                        and calculated by MOMENT then MUST be passed through
C                        FILMOM. Unchanged on exit
C
C        U        -        REAL array. A four-dimensional array, each dimension
C                        ranges from 1:KMAX. This is the fourth order product
C                        moment tensor and should be supplied by the user
C                        and calculated by MOMENT then MUST be passed through
C                        FILMOM. Unchanged on exit
C
C        BCT-CCCU        DOUBLEPRECISION arrays. These are all one-dimensional
C                        arrays, the dimension ranging 1:KMAX. On exit they
C                        contain really nothing of interest.
C
C        D,E,F                REAL arrays. Each are one-dimensional arrays, the
C                        dimension ranging 1:KMAX. On exit they contain
C                        orthonormalised versions of A,B,C.
C
C        TEXT        -        INTEGER. A flag. TEXT=1 means that information
C                        messages shall be printed. TEXT=0 means that they
C                        will not.
C
        SUBROUTINE IX3DVS(P3INDX, DP3DX, KMAX, K, N, KSTAT, DKDX, S,
     1                        DSA,DSB,DSC, SSIZE, A, B, C, T, U,
     2                  BCT,ACT,ABT,AAT,BBT,CCT,
     3                  ABCU,AACU,AABU,BBCU,ABBU,BCCU,ACCU,AAAU,
     4                        BBBU,CCCU,
     5                        D,E,F,TEXT,ERROR)


        INTEGER KMAX,K,N,SSIZE,TEXT,ERROR
        INTEGER DUMMY

C        Returned information
C
        DOUBLEPRECISION P3INDX
        DOUBLEPRECISION DP3DX(3,KMAX)
C
C        Supplied information
C
        DOUBLEPRECISION A(KMAX)
        DOUBLEPRECISION B(KMAX)
        DOUBLEPRECISION C(KMAX)
        DOUBLEPRECISION T(KMAX,KMAX,KMAX)
        DOUBLEPRECISION U(KMAX,KMAX,KMAX,KMAX)
C
C        Workspace . . .
C
C           ... that may hold useful information
C        
        DOUBLEPRECISION KSTAT(0:SSIZE, 0:SSIZE, 0:SSIZE)
        DOUBLEPRECISION DKDX(0:SSIZE, 0:SSIZE, 0:SSIZE,3,KMAX)
        DOUBLEPRECISION S(0:SSIZE, 0:SSIZE, 0:SSIZE)
        DOUBLEPRECISION DSA(0:SSIZE, 0:SSIZE, 0:SSIZE, KMAX)
        DOUBLEPRECISION DSB(0:SSIZE, 0:SSIZE, 0:SSIZE, KMAX)
        DOUBLEPRECISION DSC(0:SSIZE, 0:SSIZE, 0:SSIZE, KMAX)
        DOUBLEPRECISION D(KMAX),E(KMAX),F(KMAX)
C
C           ... probably not useful information
C
        DOUBLEPRECISION BCT(KMAX),ACT(KMAX),ABT(KMAX)
        DOUBLEPRECISION AAT(KMAX),BBT(KMAX),CCT(KMAX)
        DOUBLEPRECISION ABCU(KMAX), AACU(KMAX), AABU(KMAX), BBCU(KMAX)
        DOUBLEPRECISION ABBU(KMAX), BCCU(KMAX), ACCU(KMAX), AAAU(KMAX)
        DOUBLEPRECISION BBBU(KMAX), CCCU(KMAX)

C
C       Use TEXT to stop compilers whining
C
        DUMMY = TEXT


C
C        Orthonormalise the vectors A,B,C and store in D,E,F
C
        CALL GRAMSC(A,B,C,KMAX,K,D,E,F,ERROR)

        IF (ERROR.NE.0) THEN
                RETURN
        ENDIF

C
C        Calculate the S power-sums and put in S
C
        CALL CALCS(S, D, E, F, KMAX, K, SSIZE, T, U)
C
C        Calculate the derivatives of the S power-sums, and put in DSA,DSB,DSC
C
        CALL DERIVS(DSA, DSB, DSC, D, E, F,KMAX,K, SSIZE, S, T, U,
     1                  BCT,ACT,ABT,AAT,BBT,CCT,
     2                  ABCU,AACU,AABU,BBCU,ABBU,BCCU,ACCU,AAAU,
     3                        BBBU,CCCU)
C
C        Calculate k-statistics from S power sums, and store in KSTAT
C
        CALL STOK(KSTAT, SSIZE, S, N)
C
C        Calculate derivatives of k-statistics from derivatives of S power sums
C        and put into DKDX
C
        CALL DSTODK(DKDX, SSIZE, DSA, DSB, DSC, N, KMAX, K)

C
C        Calculate the projection index, and put into P3INDX
C
        CALL P3(P3INDX, KSTAT, SSIZE)
C
C        Calculate the derivatives, and put into DP3DX
C
        CALL CDP3DX(DP3DX, KMAX, K, KSTAT, DKDX, SSIZE)

        RETURN
        END
