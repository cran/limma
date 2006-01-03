\name{dnormexp}
\alias{dnormexp}
\alias{dnormexp.saddle}
\title{The Normal + Exponential Convolution Distribution}
\description{
Density of the normal + exponential convolution distribution.
This function is called by \code{normexp.fit}.
}
\usage{
dnormexp(x,normmean,normsd,expmean,log=FALSE)
dnormexp.saddle(x,normmean,normsd,expmean,log=FALSE,secondorder=TRUE)
}
\arguments{
  \item{x}{numeric vector of quantiles}
  \item{normmean}{mean of normal distribution}
  \item{normsd}{standard deviation of normal distribution}
  \item{expmean}{mean of exponential distribution}
  \item{log}{logical, if \code{TRUE}, density values \code{d} are given as \code{log(d)}}
  \item{secondorder}{logical, if \code{TRUE} the second order saddle-point approximation is used}
}
\details{
\code{dnormexp} computes the exact density but may fail for some extreme parameter values.
\code{dnormexp.saddle} computes the saddle-point approximation to the density.
}
\value{
numeric vector of density or log-density values.
}
\author{Jeremy Silver and Gordon Smyth}
\seealso{
An overview of background correction functions is given in \code{\link{04.Background}}.
}
\keyword{models}
