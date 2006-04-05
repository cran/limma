\name{weighted.median}
\alias{weighted.median}
\title{Weighted Median}
\description{
Compute a weighted median of a numeric vector.
}
\usage{
weighted.median(x, w, na.rm = FALSE)
}
\arguments{
\item{x}{a numeric vector containing the values whose mean is to be computed.}
\item{w}{a vector of weights the same length as \code{x} giving
the weights to use for each element of \code{x}.}
\item{na.rm}{a logical value indicating whether \code{NA}
values in \code{x} should be stripped before the computation proceeds.}
}
\details{
If \code{w} is missing then all elements of \code{x} are
given the same weight.

Missing values in \code{w} are not handled.

The weighted median is the median of the discrete distribution with
values given by \code{x} and probabilities given by \code{w/sum(w)}.
}

\value{numeric value giving the weighted median}

\seealso{
\code{\link{median}}, \code{\link{weighted.mean}}
}

\examples{
## GPA from Siegel 1994
wt <- c(5,  5,  4,  1)/15
x <- c(3.7,3.3,3.5,2.8)
xm <- weighted.median(x,wt)
}

\keyword{univar}
