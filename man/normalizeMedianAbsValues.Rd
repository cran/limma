\name{normalizeMedianAbsValues}
\alias{normalizeMedianAbsValues}
\title{Normalize Columns of a Matrix to have the Median Absolute Value}
\description{
Performs scale normalization of an M-value matrix or an A-value matrix across a series of arrays.
Users do not normally need to call these functions directly - use \code{normalizeBetweenArrays} instead.
}
\usage{
normalizeMedianAbsValues(x)
}
\arguments{
  \item{x}{numeric matrix}
}
\value{
A numeric matrix of the same size as that input which has been scaled so that each column as the same median-absolute value.
}
\details{
If \code{x} is a matrix of log-ratios of expression (M-values) then \code{normalizeMedianAbsValues} is very similar to scaling to equalize the median absolute deviation (MAD) as in Yang et al (2001, 2002).
Here the median-absolute value is used for preference to as to not re-center the M-values.

\code{normalizeMedianAbsValues} is also used to scale the A-values when scale-normalization is applied to an \code{MAList} object.
}
\author{Gordon Smyth}
\seealso{
  An overview of LIMMA functions for normalization is given in \link{05.Normalization}.
}
\examples{
M <- cbind(Array1=rnorm(10),Array2=2*rnorm(10))
normalizeMedianAbsValues(M)
}
\keyword{array}
