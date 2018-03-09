C
C        FORTRAN 77 SUBROUTINE
C        Part of the 3D Projection Pursuit Software Suite
C
C        This program is COPYRIGHT (C) 1990 G P Nason 
C
C
C        AUTHOR:        G P Nason
C
C        TITLE:        GRAMSC
C
C        CREATED:        11/3/90
C        LAST MODIFIED:        11/3/90
C        
C
C        Description
C        -----------
C
C        This subroutine is supplied with three linearly independent vectors
C        A,B and C. A set of orthonormal vectors is constructed from these
C        by the Gram-Schmidt orthonormalisation procedure and deposited in
C        the vectors D,E and F.
C
C
C        Variable Description
C        --------------------
C
C        A,B,C        -        Three linearly independent vectors
C
C        KMAX        -        Maximum array size of A,B,C
C
C        K        -        Actual array size of A,B,C
C
C        D,E,F        -        Returned orthonormal vectors.
C
        SUBROUTINE GRAMSC(A,B,C,KMAX,K,D,E,F,ERROR)
        INTEGER KMAX,K,ERROR
        DOUBLEPRECISION A(KMAX),B(KMAX),C(KMAX)
        DOUBLEPRECISION D(KMAX),E(KMAX),F(KMAX)

        DOUBLEPRECISION MODA,MODE,MODF
        DOUBLEPRECISION BTD,CTE,CTD

        INTEGER I
        DOUBLEPRECISION TOL
        PARAMETER(TOL=1.0D-6)


C
C        Calculate the modulus of A
C

        MODA = 0.0D0
        DO 10 I=1,K
                MODA = MODA + DBLE(A(I))*DBLE(A(I))
 10        CONTINUE
        MODA = SQRT(MODA)

C
C        Check that A is reasonable
C
        IF (ABS(MODA).LT.TOL) THEN
                CALL DBLEPR(
     1                  'GRAMSC: Vector A is too small, less than ',
     2                  41,TOL,1)
                ERROR = 99
                RETURN
        ENDIF

C
C        Now construct the first vector D
C        And at the same time make BTD and CTD
C

        BTD = 0.0D0
        CTD = 0.0D0
        DO 20 I=1,K
                D(I) = DBLE(A(I))/MODA
                BTD = BTD + DBLE(B(I))*DBLE(D(I))
                CTD = CTD + DBLE(C(I))*DBLE(D(I)) 
 20        CONTINUE

C
C        Now construct the second vector E
C        And calculate MODE at the same time.
C

        MODE = 0.0D0
        DO 30 I=1,K
                E(I) = B(I) - BTD*DBLE(D(I))
                MODE = MODE + DBLE(E(I))*DBLE(E(I))
 30        CONTINUE
        MODE=SQRT(MODE)
C
C        Check that E is reasonable
C
        IF (ABS(MODE).LT.TOL) THEN
                CALL DBLEPR(
     1                  'GRAMSC: Vector E is too small, less than ',
     2                  41,TOL,1)
                ERROR = 98
                RETURN
        ENDIF


C
C        Now scale the vector E so it has unit modulus
C        Also calculate CTE on the way.
C

        CTE = 0.0D0
        DO 40 I=1,K
                E(I) = DBLE(E(I))/MODE
                CTE = CTE + DBLE(C(I))*DBLE(E(I)) 
 40        CONTINUE

C
C        Now make the vector F
C        And calculate MODF at the same time
C

        MODF = 0.0D0
        DO 50 I=1,K
                F(I) = C(I) - CTE*DBLE(E(I)) - CTD*DBLE(D(I))
                MODF = MODF + DBLE(F(I))*DBLE(F(I))
 50        CONTINUE
        MODF = SQRT(MODF)
C
C        Check that F is reasonable
C
        IF (ABS(MODF).LT.TOL) THEN
                CALL DBLEPR(
     1                  'GRAMSC: Vector F is too small, less than ',
     2                  41,TOL,1)
                ERROR = 97
                RETURN
        ENDIF

C
C        Now scale F
C
        DO 60 I=1,K
                F(I) = DBLE(F(I))/MODF
 60        CONTINUE

C
C        Do ON check
C
        RETURN
        END
