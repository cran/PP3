C
C        FORTRAN 77 SUBROUTINE
C
C        This program is COPYRIGHT (C) 1991 G P Nason 
C
C
C        AUTHOR:        G P Nason
C
C        TITLE:        SPHCOR
C
C        CREATED:        08/02/91
C        LAST MODIFIED:        25/02/91
C
C        Modification History
C        --------------------
C
C        250291 gpn      Changed routine so that TRNMAT is postmultiplied
C                        by the inverse standard deviations rather than
C                        XCEN to be premultiplied by them. This means that
C                        TRNMAT can now be regarded as the transformation
C                        matrix from standard space to sphered space.
C
C        060318	gpn	 Replaced SEIGEN algorithm by call to DSYEVR which
C                        is LAPACK eigenvalue routine in R 
C        
C
C        Description
C        -----------
C
C        Spheres a multivariate data set. 
C
C        First the data set is centred to form XCEN. This is achieved by
C        routine SLOCMT.
C
C        Then the variance matrix VAR of the data is formed. This is achieved
C        by computing the outer product (XCEN*XCEN') which is KxK matrix. The
C        matrix VAR is then computed by scaling each of the outer product
C        matrix's elements by (1.0/DBLE(N)) (n.b. NOT N-1).
C
C        The correlation matrix  COR is the formed by the subroutine VARCOR.
C
C        To decorrelate the system an eigendecomposition of COR is performed
C        by DSYEVR. The inverse square root matrix of R is computed and this
C        is put into TRNMAT. The complete transformation matrix TRNMAT is
C        formed by multiplying the column i of TRNMAT by VARRii for each i.
C        Where VARii is a diagonal matrix with the iith element being
C        the reciprocal of the standard deviation of ith variable. The INVMAT
C        matrix is computed simultaneously as the inverse of TRNMAT.
C
C        The sphered data is then computed by the combination
C
C                DATMAT = (TRNMAT)*XCEN
C
C        The DATMAT matrix will have zero mean and identity variance.
C        TRNMAT represented the transformation from variable to sphered space
C        and INVMAT the inverse.
C        
C
C        Variable Description
C        --------------------
C
C        DATMAT    - DBLE array of dimension MAXROWxMAXCOL. Contains the
C                    data to be sphered. Contains the sphered data on
C                    exit (input/output).
C
C        MAXROW    - INTEGER. The maximum number of rows of the data matrix
C                    (variables). (input)
C
C        MAXCOL    - INTEGER. The maximum number of columns of the data
C                    matrix (observations). (input)
C
C        K         - INTEGER. Actual number of variables (input)
C
C        N         - INTEGER. Actual number of variables (input)
C
C        XCEN      - DBLE array of dimension MAXROWxMAXCOL. Contains the
C                    centred data matrix. (OUTPUT)
C
C        VAR       - DBLE array of dimension MAXROWxMAXROW. Contains the
C                    outer product matrix XCENxXCEN'. (OUTPUT)
C
C        COR       - DBLE array of dimension MAXROWxMAXROW. Contains the
C                    correlation matrix of the original data on OUTPUT.
C
C        EVALUE    - DBLE array of dimension MAXROW. Contains the eigenvalues
C                    of the inverse square root of the eigenvalues of the
C                    correlation matrix as computed by DSYEVR. (OUTPUT).
C
C        EV        - DBLE array of dimension MAXROW. Workspace
C
C        EET       - DBLE array of dimension MAXROWxMAXROW. Workspace
C
C        TRNMAT    - DBLE array of dimension MAXROWxMAXROW. The KxK matrix
C                    which transforms coordinates in variable space to
C                    sphered coordinates. (OUTPUT)
C
C        INVMAT    - DBLE array of dimension MAXROWxMAXROW. Inverse
C                    of TRNMAT. (OUTPUT)
C
C        ERROR     - INTEGER. Error code. 0=no error, 0<>error, see below
C
C        DOMXRW    - INTEGER. Length of ISUPPZ array (see DSYEVR help)
C
C        ISUPPZ	   - INTEGER. Array of length DOMXRW, (see DSYEVR help)
C
C        LWORK     - INTEGER. Length of WORK workspace array
C
C        WORK      - DBLE array of length LWORK (workspace)
C
C        LIWORK    - INTEGER. Length of IWORK workspace array
C
C        IWORK     - INTEGER. Array of length LIWORK (workspace)
C
C        INFO      - INTEGER. Error/information code (OUTPUT)
C
C        LAVAL     - DBLE. Array of length MAXROW, eigenvalues of
C                    correlation matrix of variables (OUTPUT)
C
C        LAVEC     - DBLE. Array of dimension MAXROWxMAXROW containing
C                    eigenvectors of correlation of variables (OUTPUT) 
C
C
C        Error codes.
C        ------------
C
C        1        -  K is larger than MAXROW.
C
C        2        -  N is larger than MAXCOL.
C
C        3        -  N is smaller than 2, we can only perform sphering
C                    on data sets with N.GE.2
C
C        4        -  No longer used
C
C        5        -  No longer used
C
C        6        -  One of the eigenvalues was less than TOL. This
C                    means that the variance matrix of the data is
C                    nearly singular. (Or worse still there are negative
C                    eigenvalues from the supposedly symmetric outer
C                    product matrix.)
C
C        10000+D  -  Error code D-10000 was returned by DSYEVR
C
C        Local Variables:
C
C        TOL      -  tolerance for real number comparisons
C        I,J        -        counters
C                RDUMMY,IDUMMY         -       not used
C                ZERO           - PARAMETER set to zero
C                ABSTOL - set to ZERO
C                M      - number of eigenvalues found


        SUBROUTINE SPHCOR(DATMAT, MAXROW, MAXCOL, K, N,
     1          XCEN, VAR, COR, EVALUE, EV, EET, TRNMAT, INVMAT,
     2          ERROR, DOMXRW, ISUPPZ, LWORK, WORK, LIWORK, IWORK,
     3          INFO, LAVAL, LAVEC)

        INTEGER MAXROW,MAXCOL,K,N,ERROR,IDUMMY,M
        INTEGER DOMXRW,LWORK,LIWORK
        INTEGER ISUPPZ(DOMXRW)
        DOUBLEPRECISION WORK(LWORK) 
        INTEGER IWORK(LIWORK)
        INTEGER INFO
        DOUBLEPRECISION LAVAL(MAXROW)
        DOUBLEPRECISION LAVEC(MAXROW, MAXROW)
        DOUBLEPRECISION DATMAT(MAXROW,MAXCOL),RDUMMY,ABSTOL,ZERO
        PARAMETER(ZERO=0.0D+0)

C
C        Tolerance for eigenvalues to be declared small or negative.
C
        DOUBLEPRECISION TOL
        PARAMETER (TOL=1e-20)

C
C        Contains the centred data matrix
C
        DOUBLEPRECISION XCEN(MAXROW, MAXCOL)
C
C        Contains the variance matrix and covariance matrix
C
        DOUBLEPRECISION VAR(MAXROW, MAXROW)        
        DOUBLEPRECISION COR(MAXROW, MAXROW)        
C
C        Eigenvalues and vectors of COR
C
        DOUBLEPRECISION EVALUE(MAXROW)
        DOUBLEPRECISION EV(MAXROW)
        DOUBLEPRECISION EET(MAXROW,MAXROW)
        DOUBLEPRECISION TRNMAT(MAXROW,MAXROW)
        DOUBLEPRECISION INVMAT(MAXROW,MAXROW)
C
C        Temporary variable to hold jth standard deviation (j=1,...,K)
C
        DOUBLEPRECISION SDJ

C        Workspace
        INTEGER I,J

C
C        Executable code begins here
C
        RDUMMY = 0.0D0
        IDUMMY = 0

C
C        Check that number of rows and columns are within limits
C
        IF (K.GT.MAXROW) THEN
                ERROR = 1
                RETURN
        ELSE IF (N.GT.MAXCOL) THEN
                ERROR = 2
                RETURN
        ENDIF

C
C        Check that we have enough observations to perform sphering
C

        IF (N.LT.2) THEN
                ERROR = 3
                RETURN
        ENDIF

C
C        Form centred data matrix XCEN
C
        CALL SLOCMT(DATMAT, MAXROW, MAXCOL, K, N, XCEN)

C
C        Now form VAR = XX'/N (first XX')
C
        CALL XXT(XCEN,MAXROW,MAXCOL,K,N,VAR, ERROR)

        IF (ERROR.NE.0) THEN
                RETURN
        ENDIF
C
C        Now form S = A/N (variance matrix)
C
        CALL SCAMAT(VAR, MAXROW, MAXROW, K, K, 1.0/DBLE(N))
C
C        Form correlation matrix C
C
        CALL VARCOR(VAR, MAXROW, MAXROW, K, K,
     1                    COR, MAXROW, MAXROW, K, K, ERROR)

        IF (ERROR.NE.0) THEN
                RETURN
        ENDIF
C
C        Do spectral decomposition of COR by LAPACK 
C
        ABSTOL = ZERO

        CALL DSYEVR('V', 'A', 'U', MAXROW, COR, MAXROW, RDUMMY, RDUMMY,
     1          IDUMMY, IDUMMY, ABSTOL, M, LAVAL, LAVEC, MAXROW, 
     2          ISUPPZ, WORK, LWORK, IWORK, LIWORK, INFO)
C
C       DSYEVR ERROR returned
C
        IF (INFO.NE.0) THEN
C
C                       Add 10000 onto ERROR code
C
                        ERROR = 10000+INFO
                        RETURN
        ENDIF


C
C       Transfer REAL->DBLE
C
C       Modification: Without SEIGEN use LAVAL, but need to reverse order
C
        DO 120 I=1,K
                EVALUE(I) = LAVAL(K-I+1)
 120        CONTINUE

C
C        Check that the eigenvalues are not too small
C
        DO 10 I=1,K
                IF (EVALUE(I).LT.TOL) THEN
                        ERROR = 6
                        RETURN
                ENDIF
 10        CONTINUE
C
C        Initialise TRNMAT and INVMAT
C
        DO 12 I=1,K
                DO 14 J=1,K
                                TRNMAT(I,J) = 0.0
                                INVMAT(I,J) = 0.0
 14                CONTINUE
 12        CONTINUE

C
C       Modification: Without SEIGEN use LAVEC, but need to reverse order
C       Now form sqrt of matrix
C
        DO 30 I=1,K
C
C                For the ith eigenvalue
C
C                form e
C                        make up eigenvector
                DO 20 J=1,K
                        EV(J) = DBLE(LAVEC(J,K-I+1))
 20                CONTINUE
C
C                form ee'
C
                CALL XXT(EV,MAXROW,1,K,1,EET,ERROR)

                IF (ERROR.NE.0) THEN
                        RETURN
                ENDIF
C
C                scale-by-reciprocal of eigenvalue
C
                CALL SCAMAT(EET, MAXROW, MAXROW, K, K,
     1                        SQRT(1.0/EVALUE(I)))
C
C                add-on-to-total
C
                CALL ADDMAT(TRNMAT, MAXROW, MAXROW, K, K, EET, MAXROW,
     1                    MAXROW, K, K, TRNMAT, MAXROW, MAXROW, ERROR)

                IF (ERROR.NE.0) THEN
                        RETURN
                ENDIF
C
C                scale EET for the inverse of TRNMAT
C
                CALL SCAMAT(EET, MAXROW, MAXROW, K, K, EVALUE(I))
C
C                Add onto INVMAT
C
                CALL ADDMAT(INVMAT, MAXROW, MAXROW, K, K, EET, MAXROW,
     1                    MAXROW, K, K, INVMAT, MAXROW, MAXROW, ERROR)

                IF (ERROR.NE.0) THEN
                        RETURN
                ENDIF

 30        CONTINUE
C
C        Multiply columns of TRNMAT by the sd inverses
C        and rows of INVMAT by the sd's 
C
C        Column J row I for TRNMAT
C        Row J Column I for INVMAT
        DO 60 J=1,K
C
C                jth variables standard deviation
C
                SDJ = SQRT(VAR(J,J))

                DO 50 I=1,K
                        TRNMAT(I,J) = TRNMAT(I,J)/SDJ
                        INVMAT(J,I) = INVMAT(J,I)*SDJ
 50                CONTINUE
 60        CONTINUE
 
C
C        Work out scaled sphered data matrix
C
        CALL SULMAT(TRNMAT,MAXROW,K,MAXROW,K,XCEN,MAXROW,K,MAXCOL,N,
     +                DATMAT,MAXROW,MAXCOL, ERROR)

        IF (ERROR.NE.0) THEN
                RETURN
        ENDIF


        ERROR = 0
        RETURN
        END
