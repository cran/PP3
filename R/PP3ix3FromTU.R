PP3ix3FromTU <-
function (the.init, avec, bvec, cvec, maxrow, k, maxcol, n, text)
{

if (length(avec) != maxrow)
	stop(paste("avec should be a vector containing ", maxrow, "elements\n"))

if (length(bvec) != maxrow)
	stop(paste("bvec should be a vector containing ", maxrow, "elements\n"))

if (length(cvec) != maxrow)
	stop(paste("cvec should be a vector containing ", maxrow, "elements\n"))


p3indx <- 0
ssize <- 4
s.wrkspc <- kstat <- array(0, dim=rep(ssize+1, 3))

dvec <- evec <- fvec <- rep(0, maxrow)

error <- 0


answer.ix3 <- .Fortran(C_ix3,
	P3INDX = as.double(p3indx),
	KMAX = as.integer(maxrow),
	K = as.integer(k),
	N = as.integer(the.init$COPYN),
	KSTAT = as.double(kstat),
	S = as.double(s.wrkspc),
	SSIZE = as.integer(ssize),
	A = as.double(avec),
	B = as.double(bvec),
	C = as.double(cvec),
	ansT = as.double(the.init$ansT),
	ansU = as.double(the.init$ansU),
	D = as.double(dvec),
	E = as.double(evec),
	F = as.double(fvec),
	TEXT = as.integer(text),
	ERROR = as.integer(error))

if (answer.ix3$ERROR != 0)
	stop(paste("ix3 FORTRAN error code was: ", answer.ix3$ERROR, "\n"))

answer <- answer.ix3$P3INDX

return(answer)
}
