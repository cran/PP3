C        G P Nason                061089
C
C        ADDMAT
C
C        Add two matrices: C=A+B
C
C        Arguments:
C                A        -        matrix to add
C                KMAX        -        matrix dimension - rows
C                NMAX        -        matrix dimension - cols
C                K        -        rows
C                N        -        cols
C                B        -        matrix to add
C                IMAX        -        matrix dimension - rows
C                JMAX        -        matrix dimension - cols
C                I        -        rows
C                J        -        cols
C                C        -        result matrix
C                CROWS        -        matrix dimension
C                CCOLS        -        matrix dimension

        SUBROUTINE ADDMAT(A,KMAX,NMAX,K,N,B,IMAX,JMAX,I,J,C,CROWS,CCOLS,
     1                ERROR)

        INTEGER KMAX,NMAX,IMAX,JMAX,K,N,I,J,CROWS,CCOLS
        DOUBLEPRECISION A(KMAX,NMAX),B(IMAX,JMAX),C(CROWS,CCOLS)
        INTEGER ERROR

        INTEGER X,Y

C
C        Check array sizes
C
      IF ((K.GT.KMAX).OR.(N.GT.NMAX).OR.(I.GT.IMAX).OR.(J.GT.JMAX)) THEN
                CALL INTPR('ADDMAT: Illegal array size',26,K,0)
                ERROR = 9
                RETURN
        ELSE IF (K.NE.I) THEN
                CALL INTPR('ADDMAT: Row dimensions not the same',
     1                        35,K,0)
                ERROR = 10
                RETURN
        ELSE IF (N.NE.J) THEN
                CALL INTPR('ADDMAT: Col dimensions not the same',
     1                        35,K,0)
                ERROR = 11
                RETURN
        ELSE IF (K.GT.CROWS) THEN
                CALL INTPR('ADDMAT: Result rows not big enough',
     1                        35,K,0)
                ERROR = 12
                RETURN
        ELSE IF (N.GT.CCOLS) THEN
                CALL INTPR('ADDMAT: Result cols not big enough',
     1                        35,K,0)
                ERROR = 13
                RETURN
        ELSE
C
C                Finally do the addition
C
                DO 20 X=1,K
                        DO 10 Y=1,N
                                C(X,Y) = A(X,Y) + B(X,Y)
 10                        CONTINUE
 20                CONTINUE
        ENDIF
        END
