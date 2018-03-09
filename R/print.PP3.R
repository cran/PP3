print.PP3 <-
function (x, ...) 
{
    cat("Class 'PP3' : Three-dimensional Projection Pursuit Object:\n")
    cat("       ~~~  : List with", length(x), "components with names\n")
    cat("             ", names(x), "\n\n")
    cat("Number of random start(s): ", length(x$ix3), "\n")
    the.maxes <- which(x$ix3 == max(x$ix3))
    cat("Maximum projection index is ", max(x$ix3), " achieved by ", length(the.maxes), " random start(s).\n")
    cat("(Partial) list of those starts achieving max are: ", the.maxes[1:min(5, length(the.maxes))])
    if (length(the.maxes) <= 5)
	    cat("\n")
    else
	    cat("...\n")
    cat("\nsummary(.):\n----------\n")
    summary.PP3(x, ...)
}
