\name{fit.normexp}
\alias{fit.normexp}
\title{Fit Normal+Exp Convolution Model to Observed Intensities}
\description{
Fit normal+exponential convolution model to observed intensities.
The normal part represents the background and the exponential represents the signal intensities.
This function is called by \code{backgroundCorrect} and is not normally called directly by the user.
}
\usage{
fit.normexp(foreground,background=NULL,background.matrix=NULL,trace=0,beta.start=NULL)
}
\arguments{
  \item{foreground}{numeric vector of foreground intensities}
  \item{background}{optional vector of background intensity values}
  \item{background.matrix}{option design matrix for regression on background values}
  \item{trace}{integer value passed to \code{optim}. If positive then tracing information on the progress of the optimization is given.  Higher values give more information.}
  \item{beta.start}{optional numeric vector giving starting values for the regression coefficients}
}
\details{
Uses Nelder-Mead simplex algorithm to maximize likelihood based on $normal(\mu,\sigma^2)+exponential(\alpha)$ convolution model for the foreground intensities.
The values $\mu$ may depend on any covariates, for example the observed background values.
}
\value{
A list containing the components
  \item{beta}{numeric vector of estimated regression coefficients}
  \item{sigma}{numeric scalar giving estimated value of $\sigma$}
  \item{alpha}{numeric scalar giving estimated value of $\alpha$}
  \item{m2loglik}{numeric scalar giving minus twice the log-likelihood}
  \item{convergence}{integer code indicating successful convergence or otherwise of the optimization. See \code{optim}.}
}
\author{Gordon Smyth}
\seealso{
\code{\link{signal.normexp}}, \code{\link[stats]{optim}}

An overview of normalization and background correction functions is given in \code{\link{4.Normalization}}.
}
\examples{
f <- c(2,3,1,10,3,20,5,6)
b <- c(2,2,2,2,2,2,2,2)
out <- fit.normexp(f,b)
}
\keyword{models}
