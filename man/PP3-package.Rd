\name{PP3-package}
\alias{PP3-package}
\alias{PP3}
\docType{package}
\title{
\packageTitle{PP3}
}
\description{
\packageDescription{PP3}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{PP3}
\packageIndices{PP3}
The main routine is \code{\link{PP3many}}. This package carries
out three-dimensional projection pursuit on multivariate data set.
It can be thought of as an alternative to principal components analysis,
where interesting views of the data are presented in a three-dimensional
projected space. This package can directly produce a true three-dimensional
solution and, not, a combination of, e.g. three one-dimensional views.
This permits the elucidation of more complex structures.
The three-dimensional solution can be used to produce colour pixel
values enabling interesting contrast display of colour images from multispectal
ones.
}
\author{
\packageAuthor{PP3}

Maintainer: \packageMaintainer{PP3}
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
\keyword{ package }
\keyword{ multivariate }
\seealso{\code{\link{PP3many}}
}
\examples{
#
# See extended example in PP3many
#
}
