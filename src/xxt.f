C        G P Nason        181089
C
C        XXT
C
C        Forms matrix product XX'
C
C
        SUBROUTINE XXT(X, MAXROW, MAXCOL, K, N, ANS, ERROR)

        INTEGER MAXROW,MAXCOL,K,N
        DOUBLEPRECISION X(MAXROW,MAXCOL), ANS(MAXROW,MAXROW)
        INTEGER ERROR

        INTEGER I,J,L

        IF ((K.GT.MAXROW).OR.(N.GT.MAXCOL)) THEN
                CALL INTPR(
     1            'XXT: Cannot multiply elements past array bounds',
     2                  45,K,0)
                ERROR =24
                RETURN
        ELSE
                CONTINUE
        ENDIF

        DO 30 I=1,K
                DO 20 J=1,K
                        ANS(I,J) = 0.0

                        DO 10 L=1,N
                                ANS(I,J) = ANS(I,J) + (X(I,L) * X(J,L))
 10                        CONTINUE
 20                CONTINUE
 30        CONTINUE

        END
