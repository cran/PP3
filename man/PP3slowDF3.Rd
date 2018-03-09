\name{PP3slowDF3}
\alias{PP3slowDF3}
\alias{PP3fastIX3}

\title{Compute the projection index or its derivative.
}
\description{Computes the projection index or its derivative with
	respect to the input projection directions.
	Function is a simple wrapper for call to \code{\link{PP3ix3FromTU}}
	or \code{\link{PP3ix3dvsFromTU}}.
}
\usage{
PP3fastIX3(Pvec, the.init, maxrow, k, maxcol, n, text)
PP3slowDF3(Pvec, the.init, maxrow, k, maxcol, n, text)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Pvec}{The projection direction. Here, a three-dimensional
	  	projection direction (matrix with three columns) is
		stacked into a single vector.
}
  \item{the.init}{Projection index initialization info. From the
	  	function \code{\link{PP3init}}
}
  \item{maxrow}{Maximum number of rows (usually equal to \code{k}) or
	  	variables.
}
  \item{k}{Actual number of rows/variables.
}
  \item{maxcol}{Maximum number of observations (usually equal to \code{n}).
}
  \item{n}{Number of observations.
}
  \item{text}{Integer. If set to 1 then the FORTRAN code prints out
	  	information messages. If set to 0, then it doesn't.
}
}
\details{\code{PP3fastIX3} computes the index only, and
	\code{PP3slowDF3} computes the 
	derivatives of the projection index with respect
		to the current projection direction (or, rather the
		Gram-Schmidt orthonormalised version). The word `slow'
		does not mean slow, but refers to the fact that this
		routine also computes the projection index, but slowly
		because the derivatives are also being computed.
}
\value{\code{PP3fastIX3} computes the projection index with respect
	to the input projection direction.
	\code{PP3slowDF3} computes a numeric vector, of the same length as \code{Pvec} containing
		the derivative of the projection index with respect to
		every entry of \code{Pvec}.
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

\seealso{\code{\link{PP3ix3dvsFromTU}}, \code{\link{PP3many}}
}
\examples{
#
# Not designed for simple user use
#
# Since these functions are simple wrappers for PP3ix3FromTU and
# PP3ix3dvsFromTU, please consult their help functions. All these
# functions do is take a single projection vector and then split it into
# three to provide three separate projection vectors for the called
# functions.
#
}
\keyword{multivariate}
