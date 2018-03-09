C        G P Nason        061089
C
C        SULMAT
C        Originally called MULMAT but this collided with an S routine
C        of the same name.
C
C        Multiplies two matrices together C=AB
C
C        Arguments:
C                A        -        first matrix
C                KMAX        -        max rows of A
C                K        -        actual rows of A to multiply
C                LMAX        -        max cols of A
C                L        -        actual cols of A to multiply
C                B        -        post multiplier matrix
C                MMAX        -        max rows of B
C                M        -        actual rows of B to multiply
C                NMAX        -        max cols of B
C                N        -        actual cols of B to multiply
C                C        -        product matrix
C                CROWS        -        max rows of C
C                CCOLS        -        max cols of C
C
C        Variables
C                H        -        integer counter
C                I        -        integer counter
C                J        -        integer counter
C
C
        SUBROUTINE SULMAT(A,KMAX,K,LMAX,L,B,MMAX,M,NMAX,N,C,CROWS,CCOLS,
     1           ERROR)
        INTEGER KMAX,K,LMAX,L,MMAX,M,NMAX,N,CROWS,CCOLS
        DOUBLEPRECISION A(KMAX,LMAX), B(MMAX,NMAX), C(CROWS, CCOLS)
        INTEGER H,I,J
        INTEGER ERROR

C
C        Check to see sizes do not exceed maximum
C
        IF ((K.GT.KMAX).OR.(L.GT.LMAX).OR.(M.GT.MMAX).OR.(N.GT.NMAX)
     +        .OR.(K.GT.CROWS).OR.(N.GT.CCOLS)) THEN
                CALL INTPR("MULMAT: Illegal array size",
     1                        26,K,0)
                ERROR = 17
                RETURN
C
C        Check that multiplication is conformable
C
        ELSE IF (L.NE.M) THEN
                CALL INTPR("MULMAT: Multiplication not conformable",
     1                        39,K,0)
                ERROR  = 18
                RETURN
        ELSE
C
C                Perform multiplication
C
                DO 30 I=1,K
                        DO 20 J=1,N
                                C(I,J)=0.0
                                DO 10 H=1,L
                                        C(I,J)=C(I,J)+ A(I,H) * B(H,J)
 10                                CONTINUE
 20                        CONTINUE
 30                CONTINUE
        END IF
        END
