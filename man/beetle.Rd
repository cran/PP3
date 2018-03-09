\name{beetle}
\alias{beetle}
\docType{data}
\title{The Beetle Data
}
\description{The Lubischew flea beetle data consisting of six
	observations measured on 74 flea beetles on three species.
	The rows dimnames attribute contains the species id.
}
\usage{data(beetle)}
\format{A 74x6 data matrix.
}
\source{
Lubischew, A.A. (1962) On the use of discriminant functions in taxonomy.
	\emph{Biometrics}, \bold{18}, 455-477.
}
\references{
Lubischew, A.A. (1962) On the use of discriminant functions in taxonomy.
	\emph{Biometrics}, \bold{18}, 455-477.
}
\examples{
data(beetle)
#
# Here is the sample mean of this multivariate data set
#
apply(beetle, 2, mean)
# [1] 177.25676 123.95946  50.35135 134.81081  12.98649  95.37838
}
\keyword{datasets}
