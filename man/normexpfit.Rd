\name{normexp.fit}
\alias{normexp.fit}
\title{Fit Normal+Exp Convolution Model to Observed Intensities}
\description{
Fit normal+exponential convolution model to observed intensities.
The normal part represents the background and the exponential represents the signal intensities.
This function is called by \code{backgroundCorrect} and is not normally called directly by the user.
}
\usage{
normexp.fit(x, trace=FALSE)
}
\arguments{
  \item{x}{numeric vector of (background corrected) intensities}
  \item{trace}{logical, if \code{TRUE}, tracing information on the progress of the optimization is given.}
}
\details{
This function uses maximum likelihood estimation to fit a model to the foreground and background intensities.
The model is a $normal(\mu,\sigma^2)+exponential(\alpha)$ convolution model for the background corrected intensities.

This is essentially the same model which is used by \code{bg.correct.rma} in the affy package.
The difference is that the parameters are estimated by maximum likelihood and that the estimated background is subtracted before applying the model-based background.
}
\value{
A list containing the components
  \item{par}{numeric vector giving estimated values of $\mu$, $\log(\sigma)$ and $\log\alpha$}
  \item{m2loglik}{numeric scalar giving minus twice the log-likelihood}
  \item{convergence}{integer code indicating successful convergence or otherwise of the optimization. See \code{optim}.}
}
\author{Jeremy Silver and Gordon Smyth}
\seealso{
\code{\link{normexp.signal}}

An overview of background correction functions is given in \code{\link{04.Background}}.
}
\examples{
f <- c(2,3,1,10,3,20,5,6)
b <- c(2,2,2,2,2,2,2,2)
out <- normexp.fit(f-b)
normexp.signal(out$par, x=f-b)
}
\keyword{models}
