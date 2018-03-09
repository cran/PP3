\name{pdataplot}
\alias{pdataplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plots the projected data for a particular random start.
}
\description{This function should be called via the
	\code{\link{plot.PP3}} function with a \code{number} argument specified.
}
\usage{
pdataplot(PP3manyobj, number, colvec = 1, lab = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{PP3manyobj}{The PP3 object you want to plot projected data for.
}
  \item{number}{The identification number that identifies which of
	the random start projection solutions you want to plot,
}
  \item{colvec}{Optional vector of colour indices to plot for each observation.
}
  \item{lab}{Labels, one for each observation. Optional.
}
  \item{\dots}{Other arguments to the plot function.
}
}
\details{Described in more detail in the help for \code{\link{plot.PP3}}
	function when the \code{number} argument is specified.
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
\author{
Friedman, J.H. and Tukey, J.W. (1974) A projection pursuit algorithm
	for exploratory data analysis. \emph{IEEE Trans. Comput.},
	\strong{23}, 881-890.

Jones, M.C. and Sibson, R. (1987) What is projection pursuit? (with discussion)
	\emph{J. R. Statist. Soc. A}, \strong{150}, 1-36.

Nason, G. P. (1995) Three-dimensional projection pursuit.
	\emph{J. R. Statist. Soc. C}, \strong{44}, 411-430.

Nason, G. P. (2001) Robust projection indices.
	\emph{J. R. Statist. Soc. B}, \strong{63}, 551-567.
%%  ~~who you are~~
}
\seealso{\code{\link{plot.PP3}}, \code{\link{PP3many}}
}
\examples{
#
# See example in help for PP3many
#
}
\keyword{multivariate}
