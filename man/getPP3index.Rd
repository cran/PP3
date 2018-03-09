\name{getPP3index}
\alias{getPP3index}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Extract projection index values from a PP3 object.
}
\description{Extract project index value vector from a PP3 object.
}
\usage{
getPP3index(x, number)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{PP3 object to get indices from.
}
  \item{number}{Identification number of solution. Optional
}
}
\details{If \code{number} is missing then the whole projection index
	vector is returned. If \code{number} is specified then the
	projection index associated with that identification number
	is returned.
}
\value{Either a vector of indices, or a single index, depending on whether
	\code{number} is not specified or it is, respectively.
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

\seealso{\code{\link{getPP3loadings}}, \code{\link{getPP3projdata}},
	\code{\link{PP3many}}
}
\examples{
#
# See help for PP3many, this contains an example of using getPP3index.
#
}
\keyword{multivariate}
