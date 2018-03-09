\name{PP3many}
\alias{PP3many}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Main function to carry out three-dimensional projection pursuit.
}
\description{
Given a multivariate data set this function applies
        exploratory projection pursuit to find an interesting
        three-dimensional projection of the input data.
}
\usage{
PP3many(xm, nrandstarts = 100, lapplyfn = lapply, action = 0, limit = 3, text = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xm}{Data matrix. Note: ROWS correspond to variables, and so the
	input matrix here is (probably) the transpose of what you might expect
	(and is accepted by principal components functions). This might
	change in the future.
}
  \item{nrandstarts}{Number of random starts.
}
  \item{lapplyfn}{By default this is lapply. If you can use the parallel
	package then you can replace this with mclapply to get a speed up.
	The latter is faster if you know how to
	use the parallel package and you often need to set the
	\code{options(mc.cores=x)} option to a value of x giving the 
	number of codes you wish to use.
}
  \item{action}{Integer taking the value 0, 1, 2 or 3. This controls
	how multivariate outliers are treated. If 0, then no action is
	taken, the multivariate set is supplied unmolested to the
	projection pursuit routine. If 1, then outliers are removed.
	If 2,3 , then outliers are moved towards the centre according
	to log or square root trimming (the precise formulae are described
	in the trimsu FORTRAN code, and stem from Jones and Sibson (1987).
	
	The outliers are only removed, or trimmed, for the purposes of
	index calculation. They are retained in the final solution and
	projected according to the projection direction decided by the
	projection pursuit algorithm computed on the non-outlier (majority)
	portion of the data.
}
  \item{limit}{This argument is only used if \code{action} is 1, 2 or 3
	and outliers need to be dealt with. If observations are greater
	in distance than \code{limit} then they are considered to be outliers.
}
  \item{text}{Integer. Has no effect if set to zero. If set to one the
	FORTRAN code will print out informational messages.
}
}
\details{
Exploratory projection pursuit is a method introduced by
Friedman and Tukey (1974) which projects multivariate K-dimensional
data down onto a L-dimensional subspace, using a projection A (an
LxK-dimensional matrix). A projection index, I, is devised to measure
how interesting the projection is and hence usually I=I(A). A numerical
optimiser is then employed to find the projection A that maximises
the index of interestingness. 

A goal of exploratory projection pursuit is to find projections that
highlight interesting structure and clustering. Often, this method
can discover interesting structures that existing methods miss, as it
uses a different measure of interestingness.

As with many other optimisation problems of this sort the numerical
optimiser does not find the global maximum, but a local maximum.
In any case, the global maximum might not correspond to the only interesting
projection direction, so it is worth exploring many projections that result
in large projection indices. One way of obtaining these is to run
projection pursuit from several random starting directions. Although this
does not guarantee to find all interesting local optima, it will find many
of them. 

This code implements the three-dimensional version of the moment projection
index proposed by Jones and Sibson, (1987).
Description of the three-dimensional version can be found in
Nason (1995). This version first centres the data (removes its mean),
and then spheres it (transforms its covariance matrix to the identity).
The aim of sphering is to ensure that any structure discovered by
projection pursuit is not related to anything that could be found by
principal components analysis (because principal components "looks for
structure" contained in the covariance matrix, and there will not be
any if the covariance matrix is the identity). The moment index is
an approximation to the entropy index, which looks for departures in
the empirical distribution of the projected data from standard normality.
The heuristic is that anything that is far from normal is `interesting'.
Hence, clustered, bi-, tri- or multimodal is not normal and hence interesting
and so clustered data is deemed to be interesting by the moment index.
The moment index has been criticised for being sensitive to outliers.
More discussion of that, and designs of robust indices can be found
in Nason (2001). However, recent numerical experimentation shows that
the moment index indeed does find projections with outliers, but these 
can be easily and quickly discounted; and the method does routinely
find interesting projections, when they are present.
}

\value{A \code{PP3} class object, which is a list with the following components.

\item{ix3}{A vector containing the optimised projection index correspond
	to each of the \code{nrandstarts} random starts.}

\item{info}{A list of \code{nrandstarts} lists. Each item in the
	list contains a list with information obtained during
	the optimisation process. The components of each list
	are \code{par} the final optimisation parameters corresponding
	to the optimal projection direction (if an optimum was found);
	\code{value} the associated optimal projection index;
	\code{counts} the number of evaluations of the projection
	index and, separately, for the gradient;
	\code{convergence} a code indicating the outcome of the optimisation;
	\code{message} a message that describes the convergence outcome.
	The last three components are information produced by the
	R \code{optim} function and are described more fully there.} 

\item{pdata.list}{The data projected according to the optimal projection
	direction resulting from each random start.
	A list containing \code{nrandstarts} matrices, each
	one corresponding to the optimisation resulting from each
	random start. Each matrix has three rows (corresponding to
	the three projection directions) and a number of columns equal
	to the number of multivariate cases (or observations, the
	number of columns of the input matrix). The
	\code{\link{plot.PP3}} function uses this information
	to produce plots of the data projected according projection
	solutions, when its \code{number} argument is specified.	}

\item{pseudp.vals}{It can be hard to evaluate a given optimal projection
	index value and know whether it can be judged as large. One
	way of assessing it is to compare it to a set of projection
	indices computed (without optimisation) on another set of
	random projection directions. It is known that most arbitrary
	projections result in uninteresting projections, so a collection
	of \code{pseudp.vals} acts as a set of projection index control
	values (NOT p-values) to which the real \code{ix3} values can
	be compared to. The \code{\link{plot.PP3}} functoin uses
	the \code{pseudp.vals} to draw a red density estimate of the
	control values (and also plots median, upper quartile, 0.9
	quantile, maximum) along with a histogram of the optimised
	\code{ix3} values, and one can easily see what should be considered
	to be large such values.
	}

\item{origvarnames}{The names of the original variables (names of rows)
	from the original data matrix. This is taken from the
	\code{dimnames} component of the data matrix.}
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

\note{Nearly all of PP3 was written by Guy Nason for his PhD (1992). 
	Nason is grateful too, and worked under the guidance of the late,
	great, Robin Sibson.
	Robin also wrote a carefully written eigendecomposition subroutine,
	which is part of this package.
}

\seealso{\code{\link{getPP3index}},\code{\link{getPP3loadings}},
	\code{\link{getPP3projdata}},\code{\link{plot.PP3}},
	\code{\link{print.PP3}}, \code{\link{summary.PP3}}
}
\examples{
#
# The flea beetle data
#
data(beetle)
#
# Run projection pursuit with 100 random starts (normally, you'd use MANY
# more random starts, e.g. 1000 or more. Here, we keep the number small to
# help CRAN
#
#
# N.b. I am going to set.seed here, so results match what you might see
# when trying THESE functions, but, in general, you can ignore set.seed
# or set it to your favourite value
#
set.seed(1)

beetle.PP3 <- PP3many(t(beetle), nrandstarts=100)
#
# Look at the output
#
beetle.PP3
#Class 'PP3' : Three-dimensional Projection Pursuit Object:
#		 ~~~  : List with 5 components with names
#		ix3 info pdata.list pseudp.vals origvarnames 
#
#Number of random start(s):  100 
#Maximum projection index is  22.02255  achieved by  1  random start(s).
#(Partial) list of those starts achieving max are:  90
#
#summary(.):
#----------
#   Summary statistics of projection index
#	 Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#       11.51   14.81   16.18   16.34   17.89   22.02 
#   Summary statistics of pseudo p-values
#	Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#	8.592  10.885  12.361  12.466  13.992  18.210 
#
# The print out shows that 100 random starts were executed and the max
# projection index was 22.02255 and only one of those random starts found
# this (sometimes more than one random start converges to the same maximum).
#
# The index number of the run which found the maximum was 90 (this number
# can be useful later to access the maximum).
#
# The summary gives the summary statistics of the 100 projection index
# values found. The max is the same as above, but also the distribution
# can be discerned.
#
# The distribution of the pseudo-p-values (NOT actual p-values) is presented
# after that. These are the projection indices computed purely on random
# directions, not the optimised versions and so you can think of them as
# null values to compare the earlier optimised values. E.g. the maximum of
# the pseudo-projection indices is 18.21, so any actual optimised projection
# index larger than this might be interesting.
#
# Now produce a plot (using all projection index info on 100 runs):
#
\dontrun{plot(beetle.PP3)}
#
# This produces (a) a histogram of the projection indices (b) a red density
# estimate of the pseudo-projection indices and (c) the median, upper quartile,
# 0.9 quantile and maximum of the pseudos as red dotted vertical lines. The red
# information corresponds to a kind of null, and so projection indices larger
# than these values might be interesting. The plot also produces some text:
#
#Big Projection Indices
#Maximum Psuedo p-value:  18.2103 
#Index Number and associated projection indices
#      90       74        4       87       75       54       13        6 
#22.02255 21.36578 21.29397 20.86531 19.59663 19.42427 19.34596 19.26520 
#      23       60 
#19.22810 19.16459
#
# This is a list of the 10 biggest projection indices and their respective
# identity numbers (which one of the random starts generated it). These
# can be used in the plot function with a number argument to generate further
# information/plots about the projection solution. Note, the number of
# biggest projection indices can be controlled with the nbig argument of
# plot.PP3.
#
# Now suppose we wanted to look at the projection solution 74, which had the
# second-biggest projection index. We can plot the projected data with the
# following command:
#
\dontrun{plot(beetle.PP3, number=74, colvec=dimnames(beetle)[[1]])}
#
# The colvec supplies the group structure so the different species can
# be coloured differently. The label argument permits you to put text
# labels there, per point as well as colours.
#
# You can extract information from the beetle.PP3 object using the
# extractor functions, getPP3index, getPP3projdata and the variable loadings
# using getPP3loading. For example,
#
getPP3index(beetle.PP3, 74)
#      74 
#21.36578 
#
# gets the second largest projection index. The third-largest can be obtained
# by replacing 74 by 4, etc.
}
\keyword{multivariate}
