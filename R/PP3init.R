PP3init <-
function (xm, action=0, limit=3.0) 
{

# xm is a data matrix
# The number of rows corresponds to the number of variables
# The number of cols corresponds to the number of observations/cases/individuals
#
# (This might be the transpose of R's usual multivariate data matrix convention)

# Note: ncol and maxrow are different deliberately (and maxcol/n)
maxrow <- k <- nrow(xm)
maxcol <- n <- ncol(xm)

datcp <- xcen <- matrix(0, nrow=maxrow, ncol=maxcol)

lavec <- invmat2 <- trnmat2 <- rcor <- invmat <- trnmat <- eet <- evectr <- corm <- varm <- matrix(0, nrow=maxrow, ncol=maxrow)

laval <- rvalue <- wrkspc <- ev <- evalue <- rep(0, maxrow)

ansT <- array(0, dim=c(maxrow, maxrow, maxrow))
ansU <- array(0, dim=c(maxrow, maxrow, maxrow, maxrow))

domxrw <- 2*maxrow
lwork <- 26*maxrow
liwork <- 10*maxrow

isuppz <- rep(0, domxrw)
work <- rep(0, lwork)
iwork <- rep(0, liwork)

copyn <- 0	# Number of observations after outliers removed, if any
error <- 0
info <- 0

#
# Sphere the data
#
answer.sphcor <- .Fortran(C_sphcor,
		   X = as.double(xm),
		   MAXROW = as.integer(maxrow),
		   MAXCOL = as.integer(maxcol),
		   K = as.integer(k),
		   N = as.integer(n),
		   XCEN = as.double(xcen),
		   VAR = as.double(varm),
		   COR = as.double(corm),
		   EVALUE = as.double(evalue),
		   EV = as.double(ev),
		   EET = as.double(eet),
		   TRNMAT = as.double(trnmat),
		   INVMAT = as.double(invmat),
		   ERROR = as.integer(error),
		   DOMXRW = as.integer(domxrw),
		   ISUPPZ = as.integer(isuppz),
		   LWORK = as.integer(lwork),
		   WORK = as.double(work),
		   LIWORK = as.integer(liwork),
		   IWORK = as.integer(iwork),
		   INFO = as.integer(info),
		   LAVAL = as.double(laval),
		   LAVEC = as.double(lavec))

if (answer.sphcor$ERROR != 0)	{
	stop(paste("FORTRAN sphcor1 error code was ", answer.sphcor$ERROR))
	}

if (answer.sphcor$INFO != 0)	{
	stop(paste("FORTRAN DSYEVR info code was ", answer.sphcor$INFO))
	}

#
# Do trimming if necessary
#
answer.trimsu <- .Fortran(C_trimsu,
		SPHDAT = as.double(answer.sphcor$X),
		K = as.integer(k),
		N = as.integer(n),
		K = as.integer(k),
		N = as.integer(n),
		LIMIT = as.double(limit),
		ACTION= as.integer(action),
		DATCP = as.double(datcp),
		COPYN = as.integer(copyn),
		ERROR = as.integer(error))

if (answer.trimsu$ERROR != 0)	{
	stop(paste("FORTRAN trimsu error code was ", answer.trimsu$ERROR))
	}

#
# Resphere
#
answer.sphcor2 <- .Fortran(C_sphcor,
		   X = as.double(answer.trimsu$DATCP),
		   MAXROW = as.integer(maxrow),
		   MAXCOL = as.integer(maxcol),
		   K = as.integer(k),
		   N = as.integer(answer.trimsu$COPYN),
		   XCEN = as.double(xcen),
		   VAR = as.double(varm),
		   COR = as.double(corm),
		   EVALUE = as.double(evalue),
		   EV = as.double(ev),
		   EET = as.double(eet),
		   TRNMAT2 = as.double(trnmat2),
		   INVMAT2 = as.double(invmat2),
		   ERROR = as.integer(error),
		   DOMXRW = as.integer(domxrw),
		   ISUPPZ = as.integer(isuppz),
		   LWORK = as.integer(lwork),
		   WORK = as.double(work),
		   LIWORK = as.integer(liwork),
		   IWORK = as.integer(iwork),
		   INFO = as.integer(info),
		   LAVAL = as.double(laval),
		   LAVEC = as.double(lavec))

if (answer.sphcor2$ERROR != 0)	{
	stop(paste("FORTRAN sphcor2 error code was ", answer.sphcor2$ERROR))
	}

if (answer.sphcor2$INFO != 0)	{
	stop(paste("FORTRAN DSYEVR info2 code was ", answer.sphcor$INFO))
	}



#
# Work out product moments
#
answer.moment <- .Fortran(C_moment,
		   DATCP = as.double(answer.sphcor2$X),
		   MAXROW = as.integer(maxrow),
		   MAXCOL = as.integer(maxcol),
		   K = as.integer(k),
		   N = as.integer(answer.trimsu$COPYN),
		   ansT = as.double(ansT),
		   ansU = as.double(ansU))


ansT <- answer.moment$ansT
ansU <- answer.moment$ansU
#
# Fill up both of ansT and ansU arrays from their triangular form to full
#
answer.filmom <- .Fortran(C_filmom,
		ansT = as.double(ansT),
		ansU = as.double(ansU),
		MAXROW = as.integer(maxrow),
		K = as.integer(k))

ansT <- answer.filmom$ansT
ansU <- answer.filmom$ansU


ans <- list(COPYN=answer.trimsu$COPYN, ansT=ansT, ansU=ansU, answer.sphcor=answer.sphcor, answer.sphcor2=answer.sphcor2)
return(ans)




}
