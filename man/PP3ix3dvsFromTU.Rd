\name{PP3ix3dvsFromTU}
\alias{PP3ix3dvsFromTU}
\alias{PP3ix3FromTU}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute the projection index or its derivatives from the T and U statistics.
}
\description{Compute the projection index or its derivatives from the T and U statistics.
}
\usage{
PP3ix3FromTU(the.init, avec, bvec, cvec, maxrow, k, maxcol, n, text)
PP3ix3dvsFromTU(the.init, avec, bvec, cvec, maxrow, k, maxcol, n, text, type = "value")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{the.init}{Initialization information, such as from \code{\link{PP3init}}.
}
  \item{avec}{First projection vector of length K
}
  \item{bvec}{Second projection vector, of length K
}
  \item{cvec}{Third projection vector, of length K
}
  \item{maxrow}{Maximum number of variables
}
  \item{k}{Actual number of variables (usually \code{maxrow} and \code{k}
	are the same.)
}
  \item{maxcol}{Maximum number of observations
}
  \item{n}{Actual number of observations. (usually \code{maxcol} and \code{n}
	are the same.)
}
  \item{text}{Integer. If equal to one then information messages from
	the FORTRAN code are produced. If equal to zero, then they are not.
}
  \item{type}{What type of information to return. If set to \code{"deriv"}
	then the derivative vector is returned. If set of \code{"value"}
	then the projection index is returned.
}
}
\details{The \code{T} and \code{U} third- and fourth-order moment related
	statistics are computed using \code{\link{PP3init}}. The
	projection index and its derivative are computed using these.
	This function can return either statistic. If you just want the
	index then you can supply the \code{"value"} argument, but
	it might be quicker to use \code{\link{PP3ix3FromTU}} which is
	called by the user-friendly \code{\link{PP3fastIX3}}.
}
\value{If \code{type} is set to \code{"value"} then a single value
	of the projection index is returned. If \code{type} is set to
	\code{"deriv"} then a vector, of length 3xK, with the derivative
	with respect to \code{avec}, \code{bvec} and \code{cvec} is
	returned.
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

\seealso{\code{\link{PP3init}}, \code{\link{PP3fastIX3}}, \code{\link{PP3slowDF3}}
}
\examples{
#
# Not for direct user use, but here is an example
#
#
# Load flea beetle data
#
data(beetle)
#
# Initialize T and U tensors
#
b.init <- PP3init(t(beetle))
#
# Get number of cases and dimensions
#
b.n <- nrow(beetle)
b.k <- ncol(beetle)
fortran.messages <- 0
#
# Select arbitrary projection vectors
#
b.pva <- c(1, rep(0, b.k-1))
b.pvb <- c(0, 1, rep(0, b.k-2))
b.pvc <- c(0, 0, 1, rep(0, b.k-3))
#
# Now compute the projection index for this data for this direction
#
answer <- PP3ix3FromTU(the.init=b.init, avec=b.pva, bvec=b.pvb,
	cvec=b.pvc, maxrow=b.k, k=b.k, maxcol=b.n, n=b.n,
	text=fortran.messages)
#
# Print out answer
#
answer
# [1] 13.49793
#
# Now compute the projection index derivatives for this data for this 
# direction
#
answer <- PP3ix3dvsFromTU(the.init=b.init, avec=b.pva, bvec=b.pvb,
	cvec=b.pvc, maxrow=b.k, k=b.k, maxcol=b.n, n=b.n,
	text=fortran.messages, type="deriv")
#
# Print out answer
#
answer
# [1]  0.000000e+00  0.000000e+00  0.000000e+00 -1.283695e-15  0.000000e+00
# [6]  0.000000e+00  0.000000e+00 -4.649059e-16  0.000000e+00  2.910680e+00
#[11] -4.941646e+00 -1.232917e+00  1.057721e-01 -4.608611e+00 -8.286708e-01
#[16] -1.602697e-01  1.724654e+00 -2.029220e+00
#
# The answer is a vector of length 3xb.k = 18. The first b.k=6 entries
# correspond to the derivative wrt b.pva, the next b.k=6 entries to
# b.pvb, and the last b.k=6 entries to b.pvc.
}
\keyword{multivariate}
