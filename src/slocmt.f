
C
C        SLOCMT        -         Slow matrix centering
C

        SUBROUTINE SLOCMT(DATA,KMAX,NMAX,K,N,XCEN)
        INTEGER KMAX,NMAX,K,N
        DOUBLEPRECISION DATA(KMAX,NMAX)
        DOUBLEPRECISION XCEN(KMAX,NMAX)

        INTEGER I,J
        DOUBLEPRECISION TOTAL,MEAN


C
C        For each variable in turn
C

        DO 10 I=1,K
C
C                Work out mean of variable
C
                TOTAL = 0.0

                DO 20 J=1,N

                        TOTAL = TOTAL + DATA(I,J)
 20                CONTINUE

                MEAN = TOTAL/N

C
C                Now subtract mean from each element
C
                DO 30 J=1,N
                        XCEN(I,J) = DATA(I,J) - MEAN
 30                CONTINUE
 10        CONTINUE
        RETURN
        END
