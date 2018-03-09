\name{summary.PP3}
\alias{summary.PP3}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Print summary information about a PP3 object.
}
\description{Print summary statistics about a PP3 object.
}
\usage{
\method{summary}{PP3}(object, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{PP3 object
}
  \item{\dots}{Other arguments (which aren't used)
}
}
\details{This applies the usual summary default function (which calculates
	summary statistics on a vector of values) to two vectors.
	The first is to the vector of maximised projection indices;
	the intention is so one can see what kinds of values the large ones
	take. The second application is to the pseudo projection indices,
	those computed on random directions without optimisation. Essentially,
	real projection indices that are larger than the maximum pseudo
	indices might be interesting and worth looking at with, e.g.
	the \code{\link{plot.PP3}} function.
}
\value{Nothing explicit is returned
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

\seealso{\code{\link{PP3many}}, \code{\link{plot.PP3}}
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
summary(beetle.PP3)
#Summary statistics of projection index
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  13.84   15.36   17.50   17.30   19.08   20.39 
#Summary statistics of pseudo p-values
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  11.14   11.78   12.77   13.59   15.31   17.63 
}
\keyword{multivariate}
