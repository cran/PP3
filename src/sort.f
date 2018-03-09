C
C        FORTRAN 77 SUBROUTINE
C        Part of the 3D Projection Pursuit Software Suite
C
C        This program is COPYRIGHT (C) 1990 G P Nason 
C
C
C        AUTHOR:        G P Nason
C
C        TITLE: SORT
C
C        CREATED:        31/1/90
C        LAST MODIFIED:        31/1/90
C        
C
C        Description
C        -----------
C
C        Sorts numbers into increasing order
C
C
C        Variable Description
C        --------------------
C        P,Q,R,S        -        Variables to be sorted
C
C        I,J,K,L        -        Sorted variables
C
C        Example
C        -------
C
C        1,4,3
C                CALL SORT3...
C        1,3,4
C
C        1,4,3,2
C                CALL SORT4...
C        1,2,3,4
C

        SUBROUTINE SORT3(P,Q,R,I,J,K)

        INTEGER P,Q,R,I,J,K
        INTEGER TEMP

        I=P
        J=Q
        K=R

        IF (I.GT.K) THEN
                TEMP = K
                K = I
                I = TEMP
        ENDIF

        IF (J.GT.K) THEN
                TEMP = K
                K = J
                J = TEMP
        ENDIF

        IF (I.GT.J) THEN
                TEMP = J
                J = I
                I = TEMP
        ENDIF
        RETURN
        END

        SUBROUTINE SORT4(P,Q,R,S,I,J,K,L)

        INTEGER P,Q,R,S,I,J,K,L
        INTEGER TEMP

        L=S
        CALL SORT3(P,Q,R,I,J,K)

        IF (L.GT.K) THEN
                RETURN
        ELSE IF (L.GT.J) THEN
                TEMP = K
                K = L
                L = TEMP
                RETURN
        ELSE IF (L.GT.I) THEN
                TEMP = L
                L = K
                K = J
                J = TEMP
                RETURN
        ELSE
                TEMP = L
                L = K
                K = J
                J = I
                I = TEMP
                RETURN
        ENDIF
        END
