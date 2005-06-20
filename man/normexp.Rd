\name{normexp}
\alias{normexp.m2loglik}
\alias{normexp.grad}
\title{Normal + Exponential Log-Likelihood}
\description{
Marginal log-likelihood of foreground values for normal + exponential model and its derivatives.
This function is called by \code{backgroundCorrect} and is not normally called directly by the user.
}
\usage{
normexp.m2loglik(theta,foreground,background=0)
normexp.grad(theta,foreground,background=0)
}
\arguments{
  \item{theta}{numeric vector of parameters}
  \item{foreground}{numeric vector of foreground intensities}
  \item{background}{optional vector of background intensity values}
}
\details{
Computes minus twice the log-likelihood based on the $normal(\mu,\sigma^2)+exponential(\alpha)$ convolution model for the foreground intensities.
The last two elements of \code{theta} are $\log(\sigma)$ and $\log(\alpha)$.
The other elements are the vector $\beta$.

The vector $\beta$ is related to $\mu$ through $\mu=X\beta$ if $X=$\code{background.matrix} is non-null.
Otherwise $\mu=\beta_1+\beta_2$\code{background} if \code{background} is non-null.
If both \code{background} and \code{background.matrix} are null, then $\mu=\beta_1$.

\code{normexp.m2loglik} is minus-twice the log-likelihood.
\code{normexp.grad} is the gradient vector.
}
\value{
Numeric scalar giving minus twice the log-likelihood
}
\author{Jeremy Silver and Gordon Smyth}
\seealso{
An overview of normalization and background correction functions is given in \code{\link{05.Normalization}}.
}
\keyword{models}
