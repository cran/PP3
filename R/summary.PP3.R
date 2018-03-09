summary.PP3 <-
function (object, ...) 
{
	cat("Summary statistics of projection index\n")
	print(summary.default(object$ix3))
	cat("Summary statistics of pseudo-indices\n")
	print(summary(object$pseudp.vals))
}
