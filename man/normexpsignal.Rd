\name{normexp.signal}
\alias{normexp.signal}
\title{Expected Signal Given Observed Foreground Under Normal+Exp Model}
\description{
Adjust foreground intensities for observed background using Normal+Exp Model.
This function is called by \code{backgroundCorrect} and is not normally called directly by the user.
}
\usage{
normexp.signal(mu,sigma,alpha,foreground)
}
\arguments{
  \item{mu}{numeric vector giving expected values for signal}
  \item{sigma}{numeric scalar giving standard deviation of background intensities}
  \item{alpha}{numeric scaler giving exponential decay rate of foreground intensities}
  \item{foreground}{numeric vector of foreground intensities}
}
\details{
In general the vector \code{mu} is computed conditional on background at each spot.
}
\value{
Numeric vector containing adjusted intensities.
}
\author{Gordon Smyth}
\seealso{
An overview of background correction functions is given in \code{\link{04.Background}}.
}
\keyword{models}
