PP3ix3dvsFromTU <-
function (the.init, avec, bvec, cvec, maxrow, k, maxcol, n, text, type="value")
{

if (length(avec) != maxrow)
	stop(paste("avec should be a vector containing ", maxrow, "elements\n"))

if (length(bvec) != maxrow)
	stop(paste("bvec should be a vector containing ", maxrow, "elements\n"))

if (length(cvec) != maxrow)
	stop(paste("cvec should be a vector containing ", maxrow, "elements\n"))


p3indx <- 0
dp3dx <- matrix(0, nrow=3, ncol=maxrow)
ssize <- 4
s.wrkspc <- kstat <- array(0, dim=rep(ssize+1, 3))
dkdx <- array(0, dim=c(ssize+1, ssize+1, ssize+1, 3, maxrow))
dsa <- dsb <- dsc <- array(0, dim=c(ssize+1, ssize+1, ssize+1, maxrow))

dvec <- evec <- fvec <- rep(0, maxrow)
bct <- act <- abt <- aat <- bbt <- cct <- rep(0, maxrow)
abcu <- aacu <- aabu <- bbcu <- abbu <- bccu <- accu <- aaau <- bbbu <- cccu <- rep(0, maxrow)

error <- 0


answer.ix3dvs <- .Fortran(C_ix3dvs,
	P3INDX = as.double(p3indx),
	DP3DX = as.double(dp3dx),
	KMAX = as.integer(maxrow),
	K = as.integer(k),
	N = as.integer(the.init$COPYN),
	KSTAT = as.double(kstat),
	DKDX = as.double(dkdx),
	S = as.double(s.wrkspc),
	DSA = as.double(dsa),
	DSB = as.double(dsb),
	DSC = as.double(dsc),
	SSIZE = as.integer(ssize),
	A = as.double(avec),
	B = as.double(bvec),
	C = as.double(cvec),
	ansT = as.double(the.init$ansT),
	ansU = as.double(the.init$ansU),
	BCT = as.double(bct),
	ACT = as.double(act),
	ABT = as.double(abt),
	AAT = as.double(aat),
	BBT = as.double(bbt),
	CCT = as.double(cct),
	ABCU = as.double(abcu),
	AACU = as.double(aacu),
	AABU = as.double(aabu),
	BBCU = as.double(bbcu),
	ABBU = as.double(abbu),
	BCCU = as.double(bccu),
	ACCU = as.double(accu),
	AAAU = as.double(aaau),
	BBBU = as.double(bbbu),
	CCCU = as.double(cccu),
	D = as.double(dvec),
	E = as.double(evec),
	F = as.double(fvec),
	TEXT = as.integer(text),
	ERROR = as.integer(error))

if (type=="value")
	answer <- answer.ix3dvs$P3INDX
else if (type=="deriv")	{
	answer <- answer.ix3dvs$DP3DX
	}
else
	stop("Unknown return type")


return(answer)
}
