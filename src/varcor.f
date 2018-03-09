C
C        G P Nason                051089
C
C        VARCOR
C
C        Converts variance matrix to correlation matrix
C
C        Arguments:
C                VAR        -        variance matrix
C                VKMAX        -        max rows of VAR
C                VNMAX        -        max cols of VAR
C                VK        -        columns of VAR to transform
C                VN        -        rows of VAR to transform
C                COR        -        correlation matrix
C                CKMAX        -        max rows of COR
C                CNMAX        -        max cols of COR
C                CK        -        cols of COR to transform
C                CN        -        rows of COR to transform
C
C        Variables:
C                I        -        counter
C                J        -        counter

        SUBROUTINE VARCOR(VAR, VKMAX, VNMAX, VK, VN, COR, CKMAX, CNMAX,
     1                CK, CN,ERROR)
        INTEGER VKMAX, VNMAX, VK, VN, CKMAX, CNMAX, CK, CN
        DOUBLEPRECISION VAR(VKMAX, VNMAX), COR(CKMAX, CNMAX)
        INTEGER ERROR
        INTEGER I,J

C
C        Check array sizes
C
        IF ((VK.GT.VKMAX).OR.(VN.GT.VNMAX).OR.(CK.GT.CKMAX).OR.
     +                (CN.GT.CNMAX)) THEN
                CALL INTPR("VARCOR: Illegal array size",24,VK,0)
                ERROR = 20
                RETURN
        ELSE IF (VK.NE.VN) THEN
                CALL INTPR("VARCOR: Variance matrix not square",
     1                        33,VK,0)
                ERROR = 21
                RETURN
        ELSE IF (CK.NE.CN) THEN
                CALL INTPR("VARCOR: Correlation matrix not square",
     1                        36,VK,0)
                ERROR = 22
                RETURN
        ELSE IF (VK.NE.CK) THEN
                CALL INTPR("VARCOR: Var. & cor. mtx not same order",
     1                        37,VK,0)
                ERROR = 23
                RETURN
        ELSE
C
C                Can produce correlation matrix
C
                DO 20 I=1,CK
                        DO 10 J= I+1,CK
                            COR(I,J) = VAR(I,J)/ SQRT(VAR(I,I)*VAR(J,J))
                            COR(J,I) = COR(I,J)
 10                        CONTINUE
                        COR(I,I) = 1.0
 20                CONTINUE
        END IF
        END
