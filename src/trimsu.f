C
C        FORTRAN 77 SUBROUTINE
C
C        This program is COPYRIGHT (C) 1991 G P Nason 
C
C
C        AUTHOR:        G P Nason
C
C        TITLE:        TRIMSU
C
C        CREATED:        13/03/91
C        LAST MODIFIED:        13/03/91
C        
C
C        Description
C        -----------
C
C        Trimming subroutine. Since projection pursuit with the moment index
C        find many projections with outliers in we need to modify some data
C        sets to remove the outlier altogether, or modify the outlier to make
C        it less of a pain.
C
C        We adapt the down-weighting methods suggested by John Tukey in
C        Jones and Sibsons' 1987 JRSS A paper `What is Projection Pursuit'.
C
C        This routine can carry out:
C
C                a) outlier analysis. The data set is not altered, however
C                possible outliers are indicated and logged;
C
C                b) outlier removal. Outliers are removed from the data set;
C
C                c) log trimming, change the outlier's distance from the origin
C                from r to 1+ln(r);
C
C                d) square root trimming, change the outlier's distance from the
C                origin from r to 3-2/(sqrt(r)).
C
C        b) is an obvious EDA thing to do; (c) and (d) were suggested by Tukey.
C
C        Whatever this routine does it is logged to fort.LOGFIL, where LOGFIL
C        is a defined parameter in this program. Generally the operation of
C        the routine is to alter the data set in DATMAT and deposit the results
C        into DATCP. However if option (a) is chosen then DATCP will contain
C        identical entries to DATMAT.
C
C
C        Variable Description
C        --------------------
C
C        DATMAT        -        DBLE array of dimension KMAXDxNMAXD. The original
C                        data matrix (input).
C
C        KMAXD        -        INTEGER. First FORTRAN dimension of DATMAT and DATCP
C                        arrays, the maximum number of variables. (input)
C
C        NMAXD        -        INTEGER. Second FORTRAN dimension of DATMAT and DATCP
C                        arrays, the maximum number of observations (input).
C
C        K        -        INTEGER. Actual number of variables (input)
C
C        N        -        INTEGER. Actual number of observations (input)
C
C        LIMIT        -        DBLE. The distance over which an observation is
C                        considered an outlier. (input).
C
C        ACTION        -        INTEGER. Type of trimming. You might like to use
C                        the parameters NOACTN,REMOVE,TRMLOG,TRMSQR below to
C                        drive this routine. They are all INTEGERs set to
C                        0,1,2, and 3 respectively. They cause the following
C                        actions to be taken:
C
C                        0        -        No trimming is done. DATCP will contain
C                                        exactly the same entries as DATMAT.
C                                        However those observations that have
C                                        distances from the origin greater than
C                                        LIMIT will have this fact logged in
C                                        fort.LOGFIL.
C
C                        1        -        Outliers are removed. DATCP will
C                                        a trimmed version of DATMAT. Note the
C                                        number of observations of DATCP may
C                                        well be less than DATMAT, and this
C                                        number is stored in COPYN. Removed
C                                        outliers are logged in fort.LOGFIL.
C
C                        2        -        Log trimming is performed. Altered
C                                        outliers are logged in fort.LOGFIL.
C
C                        3        -        Square root trimming is performed.
C                                        Altered outliers are logged in
C                                        fort.LOGFIL.
C
C        DATCP        -        DBLE array of dimension KMAXDxNMAXD. Contains the
C                        trimmed (depending on action) dataset. (output)
C
C        COPYN        -        INTEGER. The actual number of observations in DATCP.
C                        This is always LE N since we may remove some
C                        observations.        

        SUBROUTINE TRIMSU(DATMAT,KMAXD,NMAXD,K,N,LIMIT,ACTION,
     1                DATCP, COPYN, ERROR)

        INTEGER KMAXD,NMAXD,K,N,ACTION,COPYN,ERROR

        DOUBLEPRECISION LIMIT

        INTEGER NOACTN,REMOVE,TRMLOG,TRMSQR
        PARAMETER(NOACTN=0,REMOVE=1,TRMLOG=2,TRMSQR=3)

        DOUBLEPRECISION DATMAT(KMAXD, NMAXD)
        DOUBLEPRECISION DATCP(KMAXD,NMAXD)
C
C        Local variables
C
        INTEGER I,J
        DOUBLEPRECISION MODSQ
        DOUBLEPRECISION MODU
C
C        Executable code starts here
C

        IF ((ACTION.LT.0).OR.(ACTION.GT.3)) THEN
                ERROR = 96
                RETURN
        ENDIF
C
C        Initialise COPYN
C
        COPYN = 0
C
C        For each observation ...
C
C        OPEN(9,FILE='trimfile')

        DO 50 I=1,N
C
C                Find the length of the observation over the variable space
C                (MODU).
C
                MODSQ = 0.0D0

                DO 10 J=1,K

                        MODSQ = MODSQ +        DATMAT(J,I)*DATMAT(J,I)
 10                CONTINUE

                MODU = SQRT(MODSQ)
C
C                Now on the basis of MODU, the length of this observation
C                vector, decide whether to take any action. First print out
C                to trimfile what action is taken.
C

                IF (MODU.GT.LIMIT) THEN
                        CONTINUE
C                        IF (ACTION.EQ.TRMSQR) THEN
C                        WRITE(9,*),'Observation ',I,': Modulus: ',
C     1                                MODU,' Square root trimmed '
C                        ELSE IF (ACTION.EQ.TRMLOG) THEN
C                        WRITE(9,*),'Observation ',I,': Modulus: ',
C     1                                MODU,' Log trimmed '
C                        ELSE IF (ACTION.EQ.REMOVE) THEN
C                        WRITE(9,*),'Observation ',I,': Modulus: ',
C     1                                MODU,' REMOVED '
C                        ELSE
C                        WRITE(9,*),'Observation ',I,': Modulus: ',
C     1                          MODU,' NOT REMOVED '
C                        ENDIF

                ENDIF
C
C                Now take the appropriate trimming action if necessary
C
                IF ((MODU.LE.LIMIT).OR.(ACTION.EQ.NOACTN)) THEN
C
C                        Do nothing to the observation, copy it exactly from
C                        DATMAT to DATCP.
C
                        COPYN = COPYN + 1

                        DO 20 J=1,K
                                DATCP(J,COPYN) = DATMAT(J,I)
 20                        CONTINUE
                ELSE
C
C                        MODU is greater than the limit so we may have to
C                        do something
C
                        IF (ACTION.EQ.REMOVE) THEN
C
C                                Don't do anything, i.e. the observation is
C                                not copied across, so it is effectively
C                                removed.
C
                                CONTINUE
                        ELSE IF (ACTION.EQ.TRMLOG) THEN
C
C                                Log trimming
C
                                COPYN = COPYN + 1
                                DO 30 J=1,K
                                        DATCP(J,COPYN) =
     1                                                (1.0+LOG(MODU))*
     2                                                DATMAT(J,I)/MODU
 30                                CONTINUE

                        ELSE IF (ACTION.EQ.TRMSQR) THEN
C
C                                Square root trimming
C
                                COPYN = COPYN + 1
                                DO 40 J=1,K
                                        DATCP(J,COPYN)=
     1                                            (3.0-2.0/SQRT(MODU))*
     2                                            DATMAT(J,I)/MODU
 40                                CONTINUE
                        ENDIF

                ENDIF

 50        CONTINUE
C        CLOSE(9)

        RETURN
        END
