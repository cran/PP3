C
C        FORTRAN 77 SUBROUTINE
C        Part of the 3D Projection Pursuit Software Suite
C
C        This program is COPYRIGHT (C) 1990 G P Nason.
C
C
C        AUTHOR:        G P Nason
C
C        TITLE:        CALCS
C
C        CREATED:        31/1/90
C        LAST MODIFIED:        31/1/90
C        
C
C        Description
C        -----------
C
C        This program calculates the power sums S(,,) from the product
C        moments T(,,) and U(,,,). It incorporates the current direction
C        information which is specified by three k-tuple vectors A(),B()
C        and C(). These three define the projection space. Note that
C        they should always be orthonormal.
C
C        Example.
C        --------
C
C        To calculate S(2,0,1) the program works out the sum:
C
C                   K
C                   _
C        S(2,0,1) = \  A(M)*A(N)*C(P)*T(M,N,P)
C                   /
C                   _
C                M=N=P=1
C
C        To calculate S(2,1,1) the program works out the sum:
C
C                   K
C                   _
C        S(2,1,1) = \  A(M)*A(N)*B(P)*C(Q)*T(M,N,P,Q)
C                   /
C                   _
C                M=N=P=Q=1
C
C
C
C        Variable Description
C        --------------------
C
C        S        -        An array which will contain the required power
C                        sums on exit. It is of dimension SSIZE**3.
C                        This is only provisional and SSIZE should be
C                        4, since the largest subscript that we will
C                        require is 4, in S(4,0,0) etc.
C
C        A,B,C        -        The 3 projection orthonormal projection vectors
C                        that define the projection space. The FORTRAN
C                        arrays are of dimension KMAX, although only K
C                        of these values will be used.
C
C        KMAX        -        The FORTRAN dimension of arrays A,B,C
C        
C        K        -        The dimensionality of the data set.
C
C        SSIZE        -        The FORTRAN dimension of S. This should be 4.
C
C        T        -        third order product moment tensor
C
C        U        -        fourth order product moment tensor
C
C        n.b. K must be <= KMAX.
C        n.b. The arrays T,U must be full-up (i.e. not just calculated
C        with MOMENT.)
C

        SUBROUTINE CALCS(S, A, B, C, KMAX, K, SSIZE, T, U)

        INTEGER KMAX,K,SSIZE
        DOUBLEPRECISION S(0:SSIZE,0:SSIZE,0:SSIZE)
        DOUBLEPRECISION A(KMAX),B(KMAX),C(KMAX)
        DOUBLEPRECISION T(KMAX,KMAX,KMAX), U(KMAX,KMAX,KMAX,KMAX)

        INTEGER M,N,P,Q

C
C        Initialise the S array
C

C
C        Zero everything
C
        DO 3 M=0,SSIZE
          DO 2 N=0,SSIZE
            DO 1 P=0,SSIZE
                S(M,N,P) = 0.0
 1            CONTINUE
 2          CONTINUE
 3        CONTINUE

C
C        Third order - must be calculated
C
        S(1,1,1) = 0.0
        S(2,1,0) = 0.0
        S(2,0,1) = 0.0
        S(1,2,0) = 0.0
        S(1,0,2) = 0.0
        S(0,2,1) = 0.0
        S(0,1,2) = 0.0
        S(3,0,0) = 0.0
        S(0,3,0) = 0.0
        S(0,0,3) = 0.0
C
C        Fourth order - must be calculated
C
        S(2,1,1) = 0.0
        S(1,2,1) = 0.0
        S(1,1,2) = 0.0
        S(3,1,0) = 0.0
        S(3,0,1) = 0.0
        S(1,3,0) = 0.0
        S(1,0,3) = 0.0
        S(0,3,1) = 0.0
        S(0,1,3) = 0.0
        S(4,0,0) = 0.0
        S(0,4,0) = 0.0
        S(0,0,4) = 0.0
        S(0,2,2) = 0.0
        S(2,0,2) = 0.0
        S(2,2,0) = 0.0

C
C        Calculate the power sums
C

        DO 40 M=1,K
          DO 30 N=1,K
            DO 20 P=1,K
C
C                Third order
C

                S(1,1,1) = S(1,1,1) + A(M)*B(N)*C(P)*T(M,N,P)

                S(2,1,0) = S(2,1,0) + A(M)*A(N)*B(P)*T(M,N,P)

                S(2,0,1) = S(2,0,1) + A(M)*A(N)*C(P)*T(M,N,P)

                S(1,2,0) = S(1,2,0) + A(M)*B(N)*B(P)*T(M,N,P)

                S(1,0,2) = S(1,0,2) + A(M)*C(N)*C(P)*T(M,N,P)

                S(0,2,1) = S(0,2,1) + B(M)*B(N)*C(P)*T(M,N,P)

                S(0,1,2) = S(0,1,2) + B(M)*C(N)*C(P)*T(M,N,P)

                S(3,0,0) = S(3,0,0) + A(M)*A(N)*A(P)*T(M,N,P)

                S(0,3,0) = S(0,3,0) + B(M)*B(N)*B(P)*T(M,N,P)

                S(0,0,3) = S(0,0,3) + C(M)*C(N)*C(P)*T(M,N,P)
C
C                Fourth order
C
                DO 10 Q=1,K

                  S(2,1,1) = S(2,1,1) + A(M)*A(N)*B(P)*C(Q)*U(M,N,P,Q)

                  S(1,2,1) = S(1,2,1) + A(M)*B(N)*B(P)*C(Q)*U(M,N,P,Q)

                  S(1,1,2) = S(1,1,2) + A(M)*B(N)*C(P)*C(Q)*U(M,N,P,Q)

                  S(3,1,0) = S(3,1,0) + A(M)*A(N)*A(P)*B(Q)*U(M,N,P,Q)

                  S(3,0,1) = S(3,0,1) + A(M)*A(N)*A(P)*C(Q)*U(M,N,P,Q)

                  S(1,3,0) = S(1,3,0) + A(M)*B(N)*B(P)*B(Q)*U(M,N,P,Q)

                  S(1,0,3) = S(1,0,3) + A(M)*C(N)*C(P)*C(Q)*U(M,N,P,Q)

                  S(0,3,1) = S(0,3,1) + B(M)*B(N)*B(P)*C(Q)*U(M,N,P,Q)

                  S(0,1,3) = S(0,1,3) + B(M)*C(N)*C(P)*C(Q)*U(M,N,P,Q)

                  S(4,0,0) = S(4,0,0) + A(M)*A(N)*A(P)*A(Q)*U(M,N,P,Q)

                  S(0,4,0) = S(0,4,0) + B(M)*B(N)*B(P)*B(Q)*U(M,N,P,Q)

                  S(0,0,4) = S(0,0,4) + C(M)*C(N)*C(P)*C(Q)*U(M,N,P,Q)

                  S(0,2,2) = S(0,2,2) + B(M)*B(N)*C(P)*C(Q)*U(M,N,P,Q)

                  S(2,0,2) = S(2,0,2) + A(M)*A(N)*C(P)*C(Q)*U(M,N,P,Q)

                  S(2,2,0) = S(2,2,0) + A(M)*A(N)*B(P)*B(Q)*U(M,N,P,Q)
 10                CONTINUE
 20            CONTINUE
 30          CONTINUE
 40        CONTINUE

        RETURN
        END
