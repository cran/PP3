PP3many <-
function (xm, nrandstarts=100, lapplyfn=lapply, action=0, limit=3.0, text=0)
{

# xm is a data matrix
# The number of rows corresponds to the number of variables
# The number of cols corresponds to the number of observations/cases/individuals
#
# (This might be the transpose of R's usual multivariate data matrix convention)

# Note: ncol and maxrow are different deliberately (and maxcol/n)
maxrow <- k <- nrow(xm)
maxcol <- n <- ncol(xm)

if (k >= n)
	stop("Number of rows (variables) should not exceed number of columns (cases/observation). Might want to supply transpose of data matrix? See help")

origvarnames <- dimnames(xm)[[1]]



#
# Do initialization, this is sphering, trimming and moment computation
#
my.init <- PP3init(xm=xm, action=action, limit=limit)

#
# Set up list of nrandstarts random projection directions
init.par <- matrix(runif(nrandstarts*3*k), nrow=nrandstarts, ncol=3*k) 
init.par.l <- split(init.par, rep(1:nrow(init.par), each=ncol(init.par)))

#
# Generate pseudo p-values
#
the.pseudp <- function(init.par, fn, the.init, maxrow, k, maxcol, n, text){
	answer <- fn(Pvec=init.par, the.init=the.init,
		maxrow=maxrow, k=k, maxcol=maxcol, n=n, text=text) 
	return(answer)
	}

the.pseudp.list <- lapplyfn(init.par.l, the.pseudp, fn=PP3fastIX3,
	the.init=my.init, maxrow=maxrow, k=k, maxcol=maxcol, n=n, text=text)

the.pseudp.vals <- unlist(the.pseudp.list)


the.ofn <- function(init.par, fn, gr, the.init, maxrow, k, maxcol, n, text,
		    control, method)	{
	answer <- optim(par=init.par, fn=fn, gr=gr,
		the.init=the.init, maxrow=maxrow, k=k, maxcol=maxcol, n=n,
		text=text, control=control, method=method)
	return(answer)
	}

out.list <- lapplyfn(init.par.l, the.ofn, fn=PP3fastIX3, gr=PP3slowDF3,
		the.init=my.init, maxrow=maxrow, k=k, maxcol=maxcol, n=n, text=text,
		control=list(fnscale=-1), method="CG")

getindex <- function(litem)
	return(litem$value)

makepdata <- function(litem, xm)
	return(matrix(litem$par, nrow=3) %*% xm)

ix3 <- unlist(lapplyfn(out.list, getindex))

pdata.list <- lapplyfn(out.list, makepdata, xm=xm)

ans.list <- list(ix3=ix3, info=out.list, pdata.list=pdata.list, pseudp.vals=the.pseudp.vals, origvarnames=origvarnames)

class(ans.list) <- "PP3"

return(ans.list)


}
