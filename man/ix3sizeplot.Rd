\name{ix3sizeplot}
\alias{ix3sizeplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot a histogram of projection index values
}
\description{Plot a histogram of projection index values from a PP3 object,
        and superimpose control set of random projection indices for
        comparison.
}
\usage{
ix3sizeplot(PP3manyobj, main = "Projection Index Value Histogram", nbig = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{PP3manyobj}{PP3 object to plot
}
  \item{main}{Main title of plot
}
  \item{nbig}{The number of big projection indices, and their identification
	numbers to plot.
}
}
\details{See the help page for \code{\link{plot.PP3}} which calls
	this function.
}
\value{No specific return value, graphics output.
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

\seealso{\code{\link{plot.PP3}}
}
\examples{
#
# See example in PP3many, which calls plot.PP3, which calls this function
#
}
\keyword{multivariate}
