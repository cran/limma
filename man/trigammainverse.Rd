\name{trigammaInverse}
\alias{trigammaInverse}
\title{Inverse Trigamma Function}
\description{
The inverse of the trigamma function.
}
\usage{
trigammaInverse(x)
}
\arguments{
  \item{x}{numeric vector or array}
}
\details{
The function uses Newton's method with a clever starting value to ensure monotonic convergence.
}
\value{
Numeric vector or array \code{y} satisfying \code{trigamma(y)==x}.
}
\author{Gordon Smyth}
\seealso{
   \code{\link[base:Special]{trigamma}}
}
\warning{
This function does not accept a data.frame as argument although the internal function \code{trigamma} does.
}
\examples{
y <- trigammaInverse(5)
trigamma(y)
}
\keyword{math}
