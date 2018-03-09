getPP3index <-
function(x, number){

ix3 <- x$ix3
n.ix3 <- length(ix3)

if (missing(number))
	return(ix3)

else	{
	if (number < 1)
		stop("number argument has to be bigger than 0")
	else if (number > n.ix3)
		stop("number argument exceeds number of projection indices")

	return(ix3[number])

	}


}
