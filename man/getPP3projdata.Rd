\name{getPP3projdata}
\alias{getPP3projdata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Extract projected data from a PP3 object.
}
\description{Extract projected data solution from a PP3 object.
}
\usage{
getPP3projdata(PP3, number)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{PP3}{The PP3 object you wish to get information from
}
  \item{number}{The solution number you want to know about
}
}
\details{The projected data associated with a projection solution
	(out of the \code{nrandstarts} that were executed by
	\code{\link{PP3many}})) can be extracted using this function.
}
\value{The three-dimensional projected data associated with projection
	solution \code{number}. This is a matrix of three rows
	(one for each projected direction) versus N columns, where
	N is the number original cases.
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
\author{G. P. Nason}

\seealso{\code{\link{PP3many}}
}
\examples{
#
# See example for getPP3index in the help for PP3many. This function
# works similarly except that projected observations are returned.
}
\keyword{multivariate}
