C
C       FORTRAN 77 SUBROUTINE
C       Part of the 3D Projection Pursuit Software Suite
C
C       This program is COPYRIGHT (C) 1990 G P Nason.
C
C
C       AUTHOR: G P Nason
C
C       TITLE:  DERIVS
C
C       CREATED:        15/3/90
C       LAST MODIFIED:  10/4/90
C       
C        050490        gpn        Correction to sum for DSA(1,2,0,R) made. Correction
C                        was made to the multiplier of -A(R) of that sum.
C                        Corrected S(2,1,0) to the correct S(1,2,0)
C
C        100490        gpn        Correction made to DSA(0,1,2). Multiplier of S(1,0,2)
C                        changed from A(R) to B(R).
C
C        100490        gpn        Correction made to DSB(1,1,2,R). The B(R) term was
C                        multiplied by 2.0. In fact the C(R) term should be
C                        mulitplied by 2.0.
C
C        100490        gpn        Correction made to DSB(0,3,1,R). The mulitplier for
C                        C(R) was S(1,3,0) it should be (and now is) S(0,4,0)
C
C       Description
C       -----------
C
C       This program calculates the derivatives of the power sums S(,,)
C       with respect to the projection direction components A,B and C.
C       These vectors have to be an orthonormal set. The derivatives are
C       placed in the array DSA(,,,), DSB(,,,) and DSC(,,,).
C
C       Note: The arrays S,DSA,DSB, and DSC have indices from ZERO to SSIZE.
C
C
C       Variable Description
C       --------------------
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
C       A,B,C   -       The 3 orthonormal projection vectors that define the
C                       projection space. The FORTRAN arrays are of dimension
C                       KMAX although only K of these values will be used.
C
C       KMAX    -       The FORTRAN dimension of arrays A,B,C.
C
C       K       -       The dimensionality of the data set.
C
C       SSIZE   -       The FORTRAN dimension of S. This should be 4.
C
C       S       -       An array which contains the necessary power sums.
C                       It is of dimension SSIZE**3. SSIZE should be 4 as
C                       the largest subscript that we will require will be
C                       in expressions such as S(4,0,0).
C
C       T       -       third order product moment tensor
C
C       U       -       fourth order product moment tensor
C
C       BCT-CCT -       Workspace. These 6 arrays are all DOUBLEPRECISION
C                       of dimension KMAX.
C
C       ABCU-CCCU       Workspace. These 10 arrays are all DOUBLEPRECISION
C                       of dimension KMAX.
C
C       n.b. K must be <= KMAX.
C       n.b. The arrays T,U must be full-up (i.e. not just calculated
C       with MOMENT.)
C

        SUBROUTINE DERIVS(DSA, DSB, DSC, A, B, C,KMAX,K, SSIZE, S, T, U,
     1                  BCT,ACT,ABT,AAT,BBT,CCT,
     2                  ABCU,AACU,AABU,BBCU,ABBU,
     3                        BCCU,ACCU,AAAU,BBBU,CCCU)

        INTEGER KMAX,K,SSIZE
        DOUBLEPRECISION DSA(0:SSIZE, 0:SSIZE, 0:SSIZE, KMAX)
        DOUBLEPRECISION DSB(0:SSIZE, 0:SSIZE, 0:SSIZE, KMAX)
        DOUBLEPRECISION DSC(0:SSIZE, 0:SSIZE, 0:SSIZE, KMAX)
        DOUBLEPRECISION A(KMAX),B(KMAX),C(KMAX)
        DOUBLEPRECISION S(0:SSIZE, 0:SSIZE, 0:SSIZE)
        DOUBLEPRECISION T(KMAX,KMAX,KMAX), U(KMAX,KMAX,KMAX,KMAX)

        INTEGER M,N,P,R

C
C       Next set of variables used for calculating the third order
C       derivatives.
C
        DOUBLEPRECISION BCT(KMAX),ACT(KMAX),ABT(KMAX)
        DOUBLEPRECISION AAT(KMAX),BBT(KMAX),CCT(KMAX)
C
C       Next set of variables used for calculating the fourth order
C       derivatives.
C
        DOUBLEPRECISION ABCU(KMAX), AACU(KMAX), AABU(KMAX), BBCU(KMAX)
        DOUBLEPRECISION ABBU(KMAX), BCCU(KMAX), ACCU(KMAX), AAAU(KMAX)
        DOUBLEPRECISION BBBU(KMAX), CCCU(KMAX)
C
C        Initialisation of DSA,DSB and DSC
C
        DO 10 R=1,K
          DO 8  M=0,SSIZE
            DO 6 N=0,SSIZE
              DO 4 P=0,SSIZE
                DSA(M,N,P,R) = 0.0
                DSB(M,N,P,R) = 0.0
                DSC(M,N,P,R) = 0.0
 4              CONTINUE
 6            CONTINUE
 8          CONTINUE
 10        CONTINUE

C
C       We do sums like bm cn Tmnr in one loop
C

C
C       Rth derivative
C
        DO 40 R = 1,K
C
C               Initialisations
C
                BCT(R) = 0.0D0
                ACT(R) = 0.0D0
                ABT(R) = 0.0D0
                AAT(R) = 0.0D0
                BBT(R) = 0.0D0
                CCT(R) = 0.0D0
C
C               Calculations
C
                DO 30 M = 1,K
                        DO 20 N = 1,K

        BCT(R) = BCT(R) + DBLE(B(M))*DBLE(C(N))*DBLE(T(M,N,R))
        ACT(R) = ACT(R) + DBLE(A(M))*DBLE(C(N))*DBLE(T(M,N,R))
        ABT(R) = ABT(R) + DBLE(A(M))*DBLE(B(N))*DBLE(T(M,N,R))
        AAT(R) = AAT(R) + DBLE(A(M))*DBLE(A(N))*DBLE(T(M,N,R))
        BBT(R) = BBT(R) + DBLE(B(M))*DBLE(B(N))*DBLE(T(M,N,R))
        CCT(R) = CCT(R) + DBLE(C(M))*DBLE(C(N))*DBLE(T(M,N,R))

 20                     CONTINUE
 30             CONTINUE
 40     CONTINUE

C
C       Now calculate the derivatives
C
        DO 50 R=1,K

C
C               For S(1,1,1)
C
                DSA(1,1,1,R) = BCT(R) - A(R)*S(1,1,1) - B(R)*S(2,0,1)
     1                  - C(R)*S(2,1,0)

                DSB(1,1,1,R) = ACT(R) - A(R)*S(2,0,1) - B(R)*S(1,1,1)
     1                  - C(R)*S(1,2,0)

                DSC(1,1,1,R) = ABT(R) - A(R)*S(2,1,0) - B(R)*S(1,2,0)
     1                  - C(R)*S(1,1,1)

C
C               For S(2,1,0)
C
                DSA(2,1,0,R) = 2.0 * (ABT(R) - A(R)*S(2,1,0))
     1                  - B(R)*S(3,0,0)

                DSB(2,1,0,R) = AAT(R) - A(R)*S(3,0,0) - B(R)*S(2,1,0)

C               DSC(2,1,0,R) is always zero (& has already been set).

C
C               For S(2,0,1)
C
                DSA(2,0,1,R) = 2.0 * (ACT(R) - A(R)*S(2,0,1))
     1                  - C(R)*S(3,0,0)

                DSB(2,0,1,R) = - C(R)*S(2,1,0)

                DSC(2,0,1,R) = AAT(R) - A(R)*S(3,0,0) - B(R)*S(2,1,0)
     1                  - C(R)*S(2,0,1)

C
C               For S(1,2,0)
C
                DSA(1,2,0,R) = BBT(R) - A(R)*S(1,2,0)
     1                  - 2.0*B(R)*S(2,1,0)

                DSB(1,2,0,R) = 2.0*(ABT(R) - A(R)*S(2,1,0)
     1                  - B(R)*S(1,2,0))

C               DSC(1,2,0,R) is always zero (& has already been set).

C
C               For S(1,0,2)
C
                DSA(1,0,2,R) = CCT(R) - A(R)*S(1,0,2)
     1                  - 2.0*C(R)*S(2,0,1)

                DSB(1,0,2,R) = -2.0*C(R)*S(1,1,1)

                DSC(1,0,2,R) = 2.0*(ACT(R)- A(R)*S(2,0,1)
     1                  - B(R)*S(1,1,1) - C(R)*S(1,0,2))
C
C               For S(0,2,1)
C
                DSA(0,2,1,R) = -2.0*B(R)*S(1,1,1) - C(R)*S(1,2,0)

                DSB(0,2,1,R) = 2.0*(BCT(R) - A(R)*S(1,1,1)
     1                  - B(R)*S(0,2,1)) - C(R)*S(0,3,0)

                DSC(0,2,1,R) = BBT(R) - A(R)*S(1,2,0) - B(R)*S(0,3,0)
     1                  - C(R)*S(0,2,1)
C
C               For S(0,1,2)
C
                DSA(0,1,2,R) = -B(R)*S(1,0,2) - 2.0*C(R)*S(1,1,1)

                DSB(0,1,2,R) = CCT(R) - A(R)*S(1,0,2)
     1                  - B(R)*S(0,1,2) - 2.0*C(R)*S(0,2,1)

                DSC(0,1,2,R) = 2.0*(BCT(R) - A(R)*S(1,1,1)
     1                  - B(R)*S(0,2,1) - C(R)*S(0,1,2))
C
C               For S(3,0,0)
C
                DSA(3,0,0,R) = 3.0*(AAT(R) - A(R)*S(3,0,0))

C               DSB(3,0,0,R) is always zero (& has already been set).

C               DSC(3,0,0,R) is always zero (& has already been set).
C
C               For S(0,3,0)
C
                DSA(0,3,0,R) = -3.0*B(R)*S(1,2,0)

                DSB(0,3,0,R) = 3.0*(BBT(R) - A(R)*S(1,2,0)
     1                  - B(R)*S(0,3,0))

C               DSC(0,3,0,R) is always zero (& has already been set).
C
C               For S(0,0,3)
C
                DSA(0,0,3,R) = -3.0*C(R)*S(1,0,2)

                DSB(0,0,3,R) = -3.0*C(R)*S(0,1,2)

                DSC(0,0,3,R) = 3.0*(CCT(R) - A(R)*S(1,0,2)
     1                  - B(R)*S(0,1,2) - C(R)*S(0,0,3))
 50     CONTINUE
C
C       We do sums like am bn cp Umnpr in one loop
C

C
C       Rth derivative
C
        DO 90 R = 1,K
C
C         Initialisations
C
          ABCU(R) = 0.0D0
          AACU(R) = 0.0D0
          AABU(R) = 0.0D0
          BBCU(R) = 0.0D0
          ABBU(R) = 0.0D0
          BCCU(R) = 0.0D0
          ACCU(R) = 0.0D0
          AAAU(R) = 0.0D0
          BBBU(R) = 0.0D0
          CCCU(R) = 0.0D0
        
C
C         Calculations
C
          DO 80 M = 1,K
            DO 70 N = 1,K
              DO 60 P = 1,K

              ABCU(R)=ABCU(R)+DBLE(A(M))*DBLE(B(N))*
     1                  DBLE(C(P))*DBLE(U(M,N,P,R))

              AACU(R)=AACU(R)+DBLE(A(M))*DBLE(A(N))*
     1                  DBLE(C(P))*DBLE(U(M,N,P,R))

              AABU(R)=AABU(R)+DBLE(A(M))*DBLE(A(N))*
     1                  DBLE(B(P))*DBLE(U(M,N,P,R))

              BBCU(R)=BBCU(R)+DBLE(B(M))*DBLE(B(N))*
     1                  DBLE(C(P))*DBLE(U(M,N,P,R))

              ABBU(R)=ABBU(R)+DBLE(A(M))*DBLE(B(N))*
     1                  DBLE(B(P))*DBLE(U(M,N,P,R))

              BCCU(R)=BCCU(R)+DBLE(B(M))*DBLE(C(N))*
     1                  DBLE(C(P))*DBLE(U(M,N,P,R))

              ACCU(R)=ACCU(R)+DBLE(A(M))*DBLE(C(N))*
     1                  DBLE(C(P))*DBLE(U(M,N,P,R))

              AAAU(R)=AAAU(R)+DBLE(A(M))*DBLE(A(N))*
     1                  DBLE(A(P))*DBLE(U(M,N,P,R))

              BBBU(R)=BBBU(R)+DBLE(B(M))*DBLE(B(N))*
     1                  DBLE(B(P))*DBLE(U(M,N,P,R))

              CCCU(R)=CCCU(R)+DBLE(C(M))*DBLE(C(N))*
     1                  DBLE(C(P))*DBLE(U(M,N,P,R))
 60           CONTINUE
 70         CONTINUE
 80       CONTINUE
 90     CONTINUE
C
C       Now calculate the derivatives
C
        DO 100 R=1,K
C
C               For S(2,1,1)
C
                DSA(2,1,1,R) = 2.0*(ABCU(R) - A(R)*S(2,1,1))
     1                  - B(R)*S(3,0,1) - C(R)*S(3,1,0)

                DSB(2,1,1,R) = AACU(R) - A(R)*S(3,0,1)
     1                  - B(R)*S(2,1,1) - C(R)*S(2,2,0)

                DSC(2,1,1,R) = AABU(R) - A(R)*S(3,1,0)
     1                  - B(R)*S(2,2,0) - C(R)*S(2,1,1)
C
C               For S(1,2,1)
C
                DSA(1,2,1,R) = BBCU(R) - A(R)*S(1,2,1)
     1                  - 2.0*B(R)*S(2,1,1) -C(R)*S(2,2,0)

                DSB(1,2,1,R) = 2.0*(ABCU(R) - A(R)*S(2,1,1)
     1                  - B(R)*S(1,2,1)) -C(R)*S(1,3,0)

                DSC(1,2,1,R) = ABBU(R) - A(R)*S(2,2,0)
     1                  - B(R)*S(1,3,0) -C(R)*S(1,2,1)
C
C               For S(1,1,2)
C
                DSA(1,1,2,R) = BCCU(R) - A(R)*S(1,1,2)
     1                  - B(R)*S(2,0,2) -2.0*C(R)*S(2,1,1)

                DSB(1,1,2,R) = ACCU(R) - A(R)*S(2,0,2)
     1                  - B(R)*S(1,1,2) -2.0*C(R)*S(1,2,1)

                DSC(1,1,2,R) = 2.0*(ABCU(R) - A(R)*S(2,1,1)
     1                  - B(R)*S(1,2,1) - C(R)*S(1,1,2))
C
C               For S(3,1,0)
C
                DSA(3,1,0,R) = 3.0*(AABU(R) - A(R)*S(3,1,0))
     1                  - B(R)*S(4,0,0)

                DSB(3,1,0,R) = AAAU(R) - A(R)*S(4,0,0) - B(R)*S(3,1,0)

C               DSC(3,1,0,R) is always zero (& has already been set).
C
C               For S(3,0,1)
C
                DSA(3,0,1,R) = 3.0*(AACU(R) - A(R)*S(3,0,1))
     1                  - C(R)*S(4,0,0)

                DSB(3,0,1,R) = -C(R)*S(3,1,0)

                DSC(3,0,1,R) = AAAU(R) - A(R)*S(4,0,0)
     1                  - B(R)*S(3,1,0) -C(R)*S(3,0,1)
C
C               For S(1,3,0)
C
                DSA(1,3,0,R) = BBBU(R) - A(R)*S(1,3,0) -
     1                  3.0*B(R)*S(2,2,0)

                DSB(1,3,0,R) = 3.0*(ABBU(R) - A(R)*S(2,2,0)
     1                  - B(R)*S(1,3,0))

C               DSC(1,3,0,R) is always zero (& has already been set).
C
C               For S(1,0,3)
C
                DSA(1,0,3,R) = CCCU(R) - A(R)*S(1,0,3)
     1                  - 3.0*C(R)*S(2,0,2)

                DSB(1,0,3,R) = -3.0*C(R)*S(1,1,2)

                DSC(1,0,3,R) = 3.0*(ACCU(R) - A(R)*S(2,0,2)
     1                  - B(R)*S(1,1,2) -C(R)*S(1,0,3))
C
C               For S(0,3,1)
C
                DSA(0,3,1,R) = -3.0*B(R)*S(1,2,1) - C(R)*S(1,3,0)

                DSB(0,3,1,R) = 3.0*(BBCU(R) - A(R)*S(1,2,1)
     1                  - B(R)*S(0,3,1)) -C(R)*S(0,4,0)

                DSC(0,3,1,R) = BBBU(R) - A(R)*S(1,3,0)
     1                  - B(R)*S(0,4,0) -C(R)*S(0,3,1)
C
C               For S(0,1,3)
C
                DSA(0,1,3,R) = -B(R)*S(1,0,3) - 3.0*C(R)*S(1,1,2)

                DSB(0,1,3,R) = CCCU(R) - A(R)*S(1,0,3) - B(R)*S(0,1,3)
     1                  -3.0*C(R)*S(0,2,2)

                DSC(0,1,3,R) = 3.0*(BCCU(R) - A(R)*S(1,1,2)
     1                  - B(R)*S(0,2,2) -C(R)*S(0,1,3))
C
C               For S(4,0,0)
C
                DSA(4,0,0,R) = 4.0*(AAAU(R) - A(R)*S(4,0,0))

C               DSB(4,0,0,R) is always zero (& has already been set).

C               DSC(4,0,0,R) is always zero (& has already been set).
C
C               For S(0,4,0)
C
                DSA(0,4,0,R) = -4.0*B(R)*S(1,3,0)

                DSB(0,4,0,R) = 4.0*(BBBU(R) - A(R)*S(1,3,0)
     1                  - B(R)*S(0,4,0))

C               DSC(0,4,0,R) is always zero (& has already been set).
C
C               For S(0,0,4)
C
                DSA(0,0,4,R) = -4.0*C(R)*S(1,0,3)

                DSB(0,0,4,R) = -4.0*C(R)*S(0,1,3)

                DSC(0,0,4,R) = 4.0*(CCCU(R) - A(R)*S(1,0,3)
     1                  - B(R)*S(0,1,3) -C(R)*S(0,0,4))
C
C               For S(0,2,2)
C
                DSA(0,2,2,R) = -2.0*(B(R)*S(1,1,2) + C(R)*S(1,2,1))

                DSB(0,2,2,R) = 2.0*(BCCU(R) - A(R)*S(1,1,2)
     1                  - B(R)*S(0,2,2) -C(R)*S(0,3,1))

                DSC(0,2,2,R) = 2.0*(BBCU(R) - A(R)*S(1,2,1)
     1                  - B(R)*S(0,3,1) -C(R)*S(0,2,2))
C
C               For S(2,0,2)
C
                DSA(2,0,2,R) = 2.0*(ACCU(R) - A(R)*S(2,0,2)
     1                  - C(R)*S(3,0,1))

                DSB(2,0,2,R) = -2.0*C(R)*S(2,1,1)

                DSC(2,0,2,R) = 2.0*(AACU(R) - A(R)*S(3,0,1)
     1                  - B(R)*S(2,1,1) -C(R)*S(2,0,2))
C
C               For S(2,2,0)
C
                DSA(2,2,0,R) = 2.0*(ABBU(R) - A(R)*S(2,2,0)
     1                  - B(R)*S(3,1,0))

                DSB(2,2,0,R) = 2.0*(AABU(R) - A(R)*S(3,1,0)
     1                  - B(R)*S(2,2,0))

C               DSC(2,2,0,R) is always zero (& has already been set).
 100    CONTINUE
        RETURN
        END
