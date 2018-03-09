\name{plot.PP3}
\alias{plot.PP3}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Two types of plot for a PP3 object.  }
\description{Either plot projected data (if number is specified) or
	histogram of projection indices over all random starts
}
\usage{
\method{plot}{PP3}(x, number, main, nbig = 10, colvec = 1, lab = NULL, ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{PP3 object that you wish to plot.
}
  \item{number}{If missing, then the plot information concerning all
	of the project indices over all random starts 
	(using function \code{\link{ix3sizeplot}}) or, if specified,
	plots the projected data associated with the \code{number}
	random start outcome (using function
	\code{\link{pdataplot}}.
}
  \item{main}{Argument for main title.
}
  \item{nbig}{If plotting information on all indices, this argument
	controls the number of the largest projection indices and
	their identification number.
}
  \item{colvec}{If plotting the projected data, you can supply a vector
	of colour indices to colour each observation in the plot
}
  \item{lab}{As for \code{colvec}, but permits a text label for each
	observation.
}
  \item{\dots}{Other arguments to supply to the plot.
}
}
\details{This function can produce two types of plot. The first type, if
	the \code{number} argument is not supplied, produces a histogram
	of projection index values (in black) and superimposes the
	pseudo-index value density estimate, and median, upper quartile,
	0.9 quantile and maximum values as vertical dotted lines.

	If the \code{number} argument is specified, it should be one of
	the projection random start numbers and then the projected data
	associated with that maximised projection index is plotted.
	The points in the plot can be coloured (using \code{colvec}) or
	a text label supplied (using \code{lab}) argument.
}
\value{No specific value.}
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
\seealso{\code{\link{ix3sizeplot}}, \code{\link{pdataplot}},
	\code{\link{PP3many}}
}
\examples{
#
# See example in help for \code{\link{PP3many}}
#
}
\keyword{multivariate}
