\name{normexp.signal}
\alias{normexp.signal}
\title{Expected Signal Given Observed Foreground Under Normal+Exp Model}
\description{
Adjust foreground intensities for observed background using Normal+Exp Model.
This function is called by \code{backgroundCorrect} and is not normally called directly by the user.
}
\usage{
normexp.signal(par, x)
}
\arguments{
  \item{par}{numeric vector containing the parameters of the Normal+Exp distribution, see \code{normexp.fit}}
  \item{x}{numeric vector of (background corrected) intensities}
}
\details{
In general the vector \code{normmean} is computed conditional on background at each spot.
}
\value{
Numeric vector containing adjusted intensities.
}
\author{Gordon Smyth}
\seealso{
\code{\link{normexp.fit}}

An overview of background correction functions is given in \code{\link{04.Background}}.
}
\examples{# See normexp.fit}
\keyword{models}
