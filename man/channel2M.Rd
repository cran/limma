\name{designI2M}
\alias{designI2M}
\alias{designI2A}
\title{Convert Individual Channel Design Matrix to M-A Format}
\description{
Convert a design matrix in terms of individual channels to ones in terms of M-values or A-values for two-color microarray data.
}
\usage{
designI2M(design)
designI2A(design)
}
\arguments{
  \item{design}{numeric model matrix with one row for each channel observation, i.e., twice as many rows as arrays}
}

\details{
If \code{design} is a model matrix suitable for modelling individual log-intensities for two color microarray data, then \code{designI2M} computes the corresponding model matrix for modelling M-values (log-ratios) and \code{designI2A} computes the model matrix for modelling A-values (average log-intensities).

Note that the matrices \code{designI2M(design)} or \code{designI2A(design)} may be singular if not all of the coefficients are estimable from the M or A-values.
In that case there will be columns containing entirely zeros.
}

\value{
numeric model matrix with half as many rows as \code{design}
}

\seealso{
\code{\link[stats]{model.matrix}} in the stats package.

An overview of individual channel linear model functions in limma is given by \link{07.SingleChannel}.
}

\author{Gordon Smyth}

\examples{
X <- cbind(1,c(1,1,1,1,0,0,0,0),c(0,0,0,0,1,1,1,1))
designI2M(X)
designI2A(X)
}

\keyword{regression}
