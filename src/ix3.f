C
C        FORTRAN 77 SUBROUTINE
C        Part of the 3D Projection Pursuit Software Suite
C
C        This program is COPYRIGHT (C) 1990 G P Nason 
C
C
C        AUTHOR:        G P Nason
C
C        TITLE:        IX3
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
C        Calculates the projection index only, based exactly on IX3DVS under
C        a projection plane defined by A,B and C. You must also supply
C        precalculated product moments T(,,) and U(,,,).
C
C
C        Variable Description
C        --------------------
C
C        P3INDX        -        REAL. Returned value of the projection index
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
C        S        -        REAL array of workspace. A 3-dimensional array,
C                        each dimension ranging from 0:SSIZE. On exit contains
C                        the current values of the projected power-sums.
C
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
C        D,E,F                REAL arrays. Each are one-dimensional arrays, the
C                        dimension ranging 1:KMAX. On exit they contain
C                        orthonormalised versions of A,B,C.
C
C        TEXT        -        INTEGER. A flag. TEXT=1 means that information
C                        messages shall be printed. TEXT=0 means that they
C                        will not.
C
        SUBROUTINE IX3(P3INDX, KMAX, K, N, KSTAT, S,
     1                        SSIZE, A, B, C, T, U, D,E,F,TEXT,ERROR)


        INTEGER KMAX,K,N,SSIZE,TEXT,ERROR
        INTEGER DUMMY

C
C        Returned information
C
        DOUBLEPRECISION P3INDX
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
        DOUBLEPRECISION S(0:SSIZE, 0:SSIZE, 0:SSIZE)
        DOUBLEPRECISION D(KMAX),E(KMAX),F(KMAX)

C
C       Use text to stop compilers whining
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
C        Calculate k-statistics from S power sums, and store in KSTAT
C
        CALL STOK(KSTAT, SSIZE, S, N)
C
C        Calculate the projection index, and put into P3INDX
C
        CALL P3(P3INDX, KSTAT, SSIZE)

        RETURN
        END
