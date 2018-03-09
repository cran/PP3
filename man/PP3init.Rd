\name{PP3init}
\alias{PP3init}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Initialize projection pursuit code
}
\description{Initialize projection pursuit code. This function only
	need be executed once. Then subsequent calls to projection
	index calculation can reuse the results of this initialization
	many times.
}
\usage{
PP3init(xm, action = 0, limit = 3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xm}{Data matrix containing K rows (variables) and N columns
	(observations).
}
  \item{action}{Possible outlier action. See help to \code{\link{PP3many}}
	for description.
}
  \item{limit}{Possible outlier action. See help to \code{\link{PP3many}}
	for description.
}
}
\details{This function spheres the data and calculates third- and
	fourth-order moment quantities, which can be used for
	subsequent index calculation.
}
\value{A list with the following components.

\item{COPYN}{An integer less than or equal to the number of columns of
	the input matrix. The number of points that were kept after outlier
	processing (for index computation only).}

\item{ansT}{The three-dimensional T summary statistic array}
\item{ansU}{The four-dimensional U summary statistic array}
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

\seealso{\code{\link{PP3many}}
}
\examples{
#
# Not designed for direct user use, but here is an example
#
#
# The flea beetle data
#
data(beetle)
#
# Initialise the PP3 system for this data
#
tmp <- PP3init(t(beetle))
#
# This object contains the ansT and ansU third- and fourth-order tensors
#
}
\keyword{multivariate}
