\name{normexp}
\alias{normexp.m2loglik}
\alias{normexp.grad}
\title{Normal + Exponential Log-Likelihood}
\description{
Marginal log-likelihood of foreground values for normal + exponential model and its derivatives.
These functions are called by \code{normexp.fit} and are not normally called directly by the user.
}
\usage{
normexp.m2loglik(par,x)
normexp.grad(par,x)
}
\arguments{
  \item{par}{numeric vector of parameters}
  \item{x}{numeric vector of (background corrected) intensities}
}
\details{
The parameter vector \code{par} holds the normal mean, the normal log-standard deviation and the exponential mean.

Computes minus twice the log-likelihood based on the $normal(\mu,\sigma^2)+exponential(\alpha)$ convolution model for the foreground intensities.
The elements of \code{par} are $\mu$, $\log(\sigma)$ and $\log(\alpha)$.
}
\value{
\code{normexp.m2loglik} returns a numeric scalar holding minus-twice the log-likelihood.
\code{normexp.grad} returns a numeric vector holding the derivatives with respect to the elements of \code{par}.
}
\author{Jeremy Silver and Gordon Smyth}
\seealso{
An overview of background correction functions is given in \code{\link{04.Background}}.
}
\keyword{models}
