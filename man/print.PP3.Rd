\name{print.PP3}
\alias{print.PP3}
\title{Print information about a PP3 object.
}
\description{Print information about a PP3 object.
}
\usage{
\method{print}{PP3}(x, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{PP3 object to print
}
  \item{\dots}{Other arguments (not used)
}
}
\details{Prints information about a PP3 object. Primarily, the names
	of the components, the number of random starts that were used,
	the biggest (maximised) projection index that was found, and how
	many random starts achieved the max, and the index numbers of
	those. Then \code{\link{summary.PP3}} is applied.
}
\value{No specific value
}
\references{
Friedman, J.H. and Tukey, J.W. (1974) A projection pursuit algorithm
	for exploratory data analysis. \emph{IEEE Trans. Comput.},
	\strong{23}, 881-890.

Jones, M.C. and Sibson, R. (1987) What is projection pursuit? (with discussion)
	\emph{J. R. Statist. Soc. A}, \strong{150}, 1-36.

Nason, G. P. (1995) Three-dimensional projection pursuit.
	\emph{J. R. Statist. Soc. C}, \strong{44}, 411-430.

Nason, G. P. (2001) Robust projection indices.
	\emph{J. R. Statist. Soc. B}, \strong{63}, 551-567.
}
\author{G. P. Nason
}

\seealso{\code{\link{PP3many}}, \code{\link{summary.PP3}}
}
\examples{
#
# The flea beetle data
#
data(beetle)
#
# Run projection pursuit with 10 random starts (usually MUCH more than this,
# but this example will be run on installation and testing and hence I
# want to minimize computational load. A more reasonable value is 1000)
#
beetle.PP3 <- PP3many(t(beetle), nrandstarts=10)
#
# Output from summary
#
print(beetle.PP3)
#Class 'PP3' : Three-dimensional Projection Pursuit Object:
#       ~~~  : List with 5 components with names
#              ix3 info pdata.list pseudp.vals origvarnames 
#
#Number of random start(s):  10 
#Maximum projection index is  20.39497  achieved by  1  random start(s).
#(Partial) list of those starts achieving max are:  9
#
#summary(.):
#----------
#Summary statistics of projection index
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  13.84   15.36   17.50   17.30   19.08   20.39 
#Summary statistics of pseudo p-values
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  11.14   11.78   12.77   13.59   15.31   17.63 
}
\keyword{multivariate}
