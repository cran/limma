\name{m2loglik.normexp}
\alias{m2loglik.normexp}
\alias{sumloglik}
\alias{grsumloglik}
\alias{grloglik}
\title{Normal + Exponential Log-Likelihood}
\description{
Marginal log-likelihood of foreground values for normal + exponential model and its derivatives.
This function is called by \code{backgroundCorrect} and is not normally called directly by the user.
}
\usage{
m2loglik.normexp(theta,foreground,background=NULL,background.matrix=NULL)
sumloglik(theta,foreground,background)
grsumloglik(theta,foreground,background)
}
\arguments{
  \item{theta}{numeric vector of parameters}
  \item{foreground}{numeric vector of foreground intensities}
  \item{background}{optional vector of background intensity values}
  \item{background.matrix}{option design matrix for regression on background values}
}
\details{
Computes minus twice the log-likelihood based on the $normal(\mu,\sigma^2)+exponential(\alpha)$ convolution model for the foreground intensities.
The last two elements of \code{theta} are $\log(\sigma)$ and $\log(\alpha)$.
The other elements are the vector $\beta$.

The vector $\beta$ is related to $\mu$ through $\mu=X\beta$ if $X=$\code{background.matrix} is non-null.
Otherwise $\mu=\beta_1+\beta_2$\code{background} if \code{background} is non-null.
If both \code{background} and \code{background.matrix} are null, then $\mu=\beta_1$.

\code{m2loglik.normexp} is minus-twice the log-likelihood.
\code{sumloglik} is minus the log-likelihood.
\code{grsumloglik} is the gradient vector and \code{grloglik} the gradient components for individual spots.
}
\value{
Numeric scalar giving minus twice the log-likelihood
}
\section{Bugs}{
There are undocumented arguments in the function \code{grloglik}.
These will disappear soon.
}
\author{Jeremy Silver and Gordon Smyth}
\seealso{
\code{\link{fit.normexp}}, \code{\link{signal.normexp}}

An overview of normalization and background correction functions is given in \code{\link{4.Normalization}}.
}
\keyword{models}
