\name{tmixture}
\alias{tmixture.vector}
\alias{tmixture.matrix}
\title{Estimate Scale Factor in Mixture of t-Distributions}
\description{
This function estimates the unscaled standard deviation of the log fold change for differentially expressed genes.
It is called by the function \code{ebayes} and is not intended to be called by users.
}
\usage{
tmixture.vector(tstat,stdev.unscaled,df,proportion,v0.lim=NULL)
tmixture.matrix(tstat,stdev.unscaled,df,proportion,v0.lim=NULL)
}
\arguments{
  \item{tstat}{numeric vector or matrix of t-statistics}
  \item{stdev.unscaled}{numeric matrix conformal with \code{tstatf} containing the unscaled standard deviations for the coefficient estimators}
  \item{df}{numeric vector giving the degrees of freedom associated with \code{tstat}}
  \item{proportion}{assumed proportion of genes which are differentially expressed}
  \item{v0.lim}{numeric vector of length 2, assumed lower and upper limits for the estimated unscaled standard deviation}
}
\details{
The values in each column of \code{tstat} are assumed to follow a mixture of an ordinary t-distribution, with mixing proportion \code{1-proportion}, and \code{(v0+v1)/v1} times a t-distribution, with mixing proportion \code{proportion}.
Here \code{v1=stdev.unscaled^2} and \code{v0} is the value to be estimated.
}
\value{
Numeric vector of length equal to the number of columns of \code{tstat} and \code{stdev.unscaled}.
}
\seealso{
\code{\link{ebayes}}
}
\author{Gordon Smyth}
\keyword{htest}
