\name{m2loglik.normexp}
\alias{m2loglik.normexp}
\title{Normal + Exponential Log-Likelihood}
\description{
Marginal log-likelihood of foreground values for normal + exponential model.
This function is called by \code{backgroundCorrect} and is not normally called directly by the user.
}
\usage{
m2loglik.normexp(theta,foreground,background=NULL,background.matrix=NULL)
}
\arguments{
  \item{theta}{numeric vector of parameters}
  \item{foreground}{numeric vector of foreground intensities}
  \item{background}{optional vector of background intensity values}
  \item{background.matrix}{option design matrix for regression on background values}
}
\details{
Computes the likelihood based on the $normal(\mu,\sigma^2)+exponential(\alpha)$ convolution model for the foreground intensities.
The last two elements of \code{theta} are $\log(\sigma)$ and $\log(\alpha)$.
The other elements are the vector $\beta$.

The vector $\beta$ is related to $\mu$ through $\mu=X\beta$ if $X=$\code{background.matrix} is non-null.
Otherwise $\mu=\beta_1+\beta_2$\code{background} if \code{background} is non-null.
If both \code{background} and \code{background.matrix} are null, then $\mu=\beta_1$.
}
\value{
Numeric scalar giving minus twice the log-likelihood
}
\author{Gordon Smyth}
\seealso{
\code{\link{fit.normexp}}, \code{\link{signal.normexp}}

An overview of normalization and background correction functions is given in \code{\link{4.Normalization}}.
}
\keyword{models}
