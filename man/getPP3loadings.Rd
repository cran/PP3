\name{getPP3loadings}
\alias{getPP3loadings}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Extract projection loadings from a PP3 object.
}
\description{Extract projection loadings (directions) from a PP3 object.
}
\usage{
getPP3loadings(PP3, number)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{PP3}{The PP3 object.
}
  \item{number}{which projection solution you want to get.
}
}
\details{Gets projection directions associated with a particular solution
from the \code{nrandstarts} that \code{\link{PP3many}} produces.
}
\value{A matrix of dimensional 3xK, where K was the original dimensionality
	of the input data set. Note: they are not an orthogonal set,
	but the back transform of an orthogonal set (back transformed
	from the inverse sphering transform).
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

\seealso{\code{\link{getPP3index}}, \code{\link{PP3many}}
}
\examples{
#
# See example for getPP3index in the help for PP3many, this function
# works similarly but returns loadings rather than indices.
#
}
\keyword{multivariate}
