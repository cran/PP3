C
C        FORTRAN 77 SUBROUTINE
C        Part of the 3D Projection Pursuit Software Suite
C
C        This program is COPYRIGHT (C) 1990 G P Nason 
C
C
C        AUTHOR:        G P Nason
C
C        TITLE:        MOMENT
C
C        CREATED:        25/1/90
C        LAST MODIFIED:        25/1/90
C        
C
C        Description
C        -----------
C
C        This subroutine calculates the third and fourth product moments
C        for a set of data X. X is a K by N data matrix.
C
C        X is subscripted as follows X(R,I) which refers to the ith case
C        on the rth variable.
C
C        Two arrays must be supplied. These are T and U.
C
C        Array T is the third moment array. It is a 3dimensional array
C        with K variables on each dimension. Array U is the fourth moment
C        array. It is a four dimensional array with K variables on each
C        dimension.
C
C        T and U are obtained as follows:
C
C                  N
C                   _
C        T(p,q,r) = \ X(p,i)*X(q,i)*X(r,i)
C                   /
C                   _
C                  i=1
C
C                     N
C                     _
C        U(p,q,r,s) = \ X(p,i)*X(q,i)*X(r,i)*X(s,i)
C                     /
C                     _
C                    i=1
C
C
C        Variable Description
C        --------------------
C
C        X        Data matrix of true dimension (KMAX,NMAX) but only
C                (K,N) will actually be used.
C
C        NMAX        Maximum dimension of X (max. number of cases).
C
C        KMAX        Maximum dimension of X (max. number of variables).
C
C        N        Number of cases actually used.
C
C        K        Number of variables actually used.
C
C        T        Third order product moment array of dimension (KMAX,KMAX,KMAX)
C
C        U        Fourth order product moment array of dim. (KMAX,KMAX,KMAX,KMAX)        
C
        SUBROUTINE MOMENT(X, KMAX, NMAX, K, N, T, U)

C
C        Argument declarations
C
        INTEGER NMAX,KMAX,N,K

        DOUBLEPRECISION X(KMAX,NMAX)
        DOUBLEPRECISION T(KMAX,KMAX,KMAX)
        DOUBLEPRECISION U(KMAX,KMAX,KMAX,KMAX)

C
C        Local program variables
C

        INTEGER P,Q,R,S,I
        DOUBLEPRECISION TEMP

C
C        Routine starts here
C

        DO 50 P = 1,K
          DO 40 Q = P,K
            DO 30 R = Q,K
C
C                Initialise T
C
                TEMP = 0.0D0
C
C                Create the T moments
C
                DO 20 I = 1,N
                        TEMP=TEMP+DBLE(X(P,I))*DBLE(X(Q,I))*DBLE(X(R,I))
 20                CONTINUE

                T(P,Q,R) = DBLE(TEMP)

                DO 10 S = R,K
C
C                        Initialise U
C
                        TEMP = 0.0D0 
C
C                        Create the U moments
C
                        DO 5 I = 1,N
                          TEMP=TEMP+DBLE(X(P,I))*DBLE(X(Q,I))
     1                                *DBLE(X(R,I))*DBLE(X(S,I))
 5                        CONTINUE
                        U(P,Q,R,S) = DBLE(TEMP)
 10                CONTINUE
 30            CONTINUE
 40          CONTINUE
 50        CONTINUE

        RETURN
        END
