C
C        FORTRAN 77 SUBROUTINE
C        Part of the 3D Projection Pursuit Software Suite
C
C        This program is COPYRIGHT (C) 1990 G P Nason 
C
C
C        AUTHOR:        G P Nason
C
C        TITLE:        FILMOM
C
C        CREATED:        31/1/90
C        LAST MODIFIED:        31/1/90
C        
C
C        Description
C        -----------
C
C        The MOMENT subroutine calculates the third and fourth order
C        product moment tensors T and U. However these arrays are totally
C        symmetric in all their arguments and so only the "upper triangular"
C        parts of these arrays are actually calculated.
C
C        Later calculations may use T and U in an unpredictable way and
C        so the arrays should contain all the values.
C
C
C        Variable Description
C        --------------------
C

        SUBROUTINE FILMOM(T, U, KMAX, K)

        INTEGER KMAX,K
        DOUBLEPRECISION T(KMAX,KMAX,KMAX),U(KMAX,KMAX,KMAX,KMAX)

        INTEGER P,Q,R,S,I,J,M,N

        DO 40 P=1,K
          DO 30 Q=1,K
            DO 20 R=1,K

                CALL SORT3(P,Q,R,I,J,M)

                T(P,Q,R) = T(I,J,M)

C                PRINT *,'T(',P,',',Q,',',R,')=',T(P,Q,R)

                DO 10 S=1,K

                        CALL SORT4(P,Q,R,S,I,J,M,N)

                        U(P,Q,R,S) = U(I,J,M,N)

C                        PRINT *,'U(',P,',',Q,',',R,',',S,')=',U(P,Q,R,S)
 10                CONTINUE
 20            CONTINUE
 30          CONTINUE
 40        CONTINUE
        RETURN
        END
