\name{normalizeScale}
\alias{normalizeMedianDeviations}
\alias{normalizeMedians}
\title{Normalize Columns of a Matrix to have the Same Scale}
\description{
Performs scale normalization of an M-value matrix or an A-value matrix across a series of arrays.
Users do not normally need to call these functions directly - use \code{normalizeBetweenArrays} instead.
}
\usage{
normalizeMedianDeviations(x)
normalizeMedians(x)
}
\arguments{
  \item{x}{numeric matrix}
}
\value{
\code{normalizeMedianDeviations} produces a numeric matrix of the same size as that input which has been scaled so that each column as the same median-absolute value.

\code{normalizeMedians} produces a numeric matrix which has been scaled so that each column has the same median-value.
}
\details{
If \code{x} is a matrix of log-ratios of expression (M-values) then \code{normalizeMedianDeviations} is very similar to scaling to equalize the median absolute deviation (MAD) as in Yang et al (2001, 2002).
Here the median-absolute value is used for preference to as to not re-center the M-values.

\code{normalizeMedians} is used for A-values of overall expression.
}
\author{Gordon Smyth}
\seealso{
  An overview of LIMMA functions for normalization is given in \link{4.Normalization}.
}
\examples{
M <- cbind(Array1=rnorm(10),Array2=2*rnorm(10))
normalizeMedianDeviations(M)

A <- cbind(Array1=rlnorm(10),Array2=2*rlnorm(10))
normalizeMedians(A)
}
\keyword{array}
