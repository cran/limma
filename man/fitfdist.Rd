\name{fitFDist}
\alias{fitFDist}
\title{Moment Estimation of Scaled F-Distribution}
\description{
Moment estimation of the parameters of a scaled F-distribution given one of the degrees of freedom.
This function is called internally by \code{ebayes} and is not usually called directly by a user.
}
\usage{
fitFDist(x,df1)
}
\arguments{
  \item{x}{numeric vector or array of positive values representing a sample from an F-distribution.}
  \item{df1}{the first degrees of freedom of the F-distribution. May be an integer or a vector of the same length as \code{x}.}
}
\details{
The function estimates \code{scale} and \code{df2} under the assumption that \code{x} is distributed as \code{scale} times an F-distributed random variable on \code{df1} and \code{df2} degrees of freedom.
}
\value{
A list containing the components
  \item{scale}{scale factor for F-distribution}
  \item{df2}{the second degrees of freedom of the F-distribution}
}
\author{Gordon Smyth}
\seealso{
\code{\link{ebayes}}, \code{\link{trigammaInverse}}
}
\keyword{distribution}
