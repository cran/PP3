C        G P Nason                051089
C
C        SCAMAT
C
C        Scales a matrix by a (real) scalar quantity
C
C        Arguments:
C                X        -        matrix to scale
C                KMAX        -        max rows of matrix
C                NMAX        -        max cols of matrix
C                K        -        number of rows of matrix X
C                N        -        number of columns of matrix X
C                R        -        scalar multiplier
C
C        Variables:
C                I        -        counter
C                J        -        counter
C
        SUBROUTINE SCAMAT(X, KMAX, NMAX, K, N, R)
        INTEGER KMAX, NMAX, K, N
        DOUBLEPRECISION X(KMAX, NMAX), R

        INTEGER I,J
C
C        Check dimensions
C
        IF ((K.GT.KMAX).OR.(N.GT.NMAX)) THEN
                CALL INTPR("SCAMAT: Illegal array size",27,K,0)
        ELSE
C
C                Scale the matrix
C
                DO 20 I=1,K
                        DO 10 J=1,N
                                X(I,J) = R * X(I,J)
 10                        CONTINUE
 20                CONTINUE
        END IF
        END
