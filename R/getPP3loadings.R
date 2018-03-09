getPP3loadings <-
function(PP3, number){

if (missing(number))
	stop("You have to specify which set of projection directions you want by specifying a number argument")

if (number < 1)
	stop("Number has to be bigger than 0\n")
else if (number > length(PP3$ix3))
	stop("Number has to be smaller than ", length(PP3$ix3), "\n")

m <- matrix(PP3$info[[number]]$par, nrow=3)
dimnames(m)[[2]] <- PP3$origvarnames
dimnames(m)[[1]] <- c("Comp. 1", "Comp. 2", "Comp. 3")

return(m)
}
