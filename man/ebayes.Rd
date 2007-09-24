\name{ebayes}
\alias{ebayes}
\alias{eBayes}
\title{Empirical Bayes Statistics for Differential Expression}
\description{Given a series of related parameter estimates and standard errors, compute moderated t-statistics, moderated F-statistic, and log-odds of differential expression by empirical Bayes shrinkage of the standard errors towards a common value.}
\usage{
ebayes(fit,proportion=0.01,stdev.coef.lim=c(0.1,4))
eBayes(fit,proportion=0.01,stdev.coef.lim=c(0.1,4))
}
\arguments{
  \item{fit}{an \code{MArrayLM} fitted model object produced by \code{lmFit} or \code{contrasts.fit}, or an unclassed list produced by \code{lm.series}, \code{gls.series} or \code{mrlm} containing components \code{coefficients}, \code{stdev.unscaled}, \code{sigma} and \code{df.residual}}
  \item{proportion}{numeric value between 0 and 1, assumed proportion of genes which are differentially expressed}
  \item{stdev.coef.lim}{numeric vector of length 2, assumed lower and upper limits for the standard deviation of log2 fold changes for differentially expressed genes}
}
\value{
\code{ebayes} produces an ordinary list with the following components.
\code{eBayes} adds the following components to \code{fit} to produce an augmented object, usually of class \code{MArrayLM}.
  \item{t}{numeric vector or matrix of moderated t-statistics}
  \item{p.value}{numeric vector of p-values corresponding to the t-statistics}
  \item{s2.prior}{estimated prior value for \code{sigma^2}}
  \item{df.prior}{degrees of freedom associated with \code{s2.prior}}
  \item{s2.post}{vector giving the posterior values for \code{sigma^2}}
  \item{lods}{numeric vector or matrix giving the log-odds of differential expression}
  \item{var.prior}{estimated prior value for the variance of the log2-fold-change for differentially expressed gene}
  \item{F}{numeric vector of moderated F-statistics for testing all contrasts defined by the columns of \code{fit} simultaneously equal to zero}
  \item{F.p.value}{numeric vector giving p-values corresponding to \code{F}}
}
\details{
These functions is used to rank genes in order of evidence for differential expression.
They use an empirical Bayes method to shrink the probe-wise sample variances towards a common value and to augmenting the degrees of freedom for the individual variances (Smyth, 2004).
The functions accept as input argument \code{fit} a fitted model object from the functions \code{lmFit}, \code{lm.series}, \code{mrlm} or \code{gls.series}.
The fitted model object may have been processed by \code{contrasts.fit} before being passed to \code{eBayes} to convert the coefficients of the design matrix into an arbitrary number of contrasts which are to be tested equal to zero.
The columns of \code{fit} define a set of contrasts which are to be tested equal to zero.

The empirical Bayes moderated t-statistics test each individual contrast equal to zero.
For each probe (row), the moderated F-statistic tests whether all the contrasts are zero.
The F-statistic is an overall test computed from the set of t-statistics for that probe.
This is exactly analogous the relationship between t-tests and F-statistics in conventional anova, except that the residual mean squares and residual degrees of freedom have been moderated between probes.

The estimates \code{s2.prior} and \code{df.prior} are computed by \code{fitFDist}.
\code{s2.post} is the weighted average of \code{s2.prior} and \code{sigma^2} with weights proportional to \code{df.prior} and \code{df.residual} respectively.
The \code{lods} is sometimes known as the B-statistic.
The F-statistics \code{F} are computed by \code{classifyTestsF} with \code{fstat.only=TRUE}.

\code{eBayes} doesn't compute ordinary (unmoderated) t-statistics by default, but these can be easily extracted from 
the linear model output, see the example below.

\code{ebayes} is the earlier and leaner function.
\code{eBayes} is intended to have a more object-orientated flavor as it produces objects containing all the necessary components for downstream analysis.
}
\seealso{
\code{\link{squeezeVar}}, \code{\link{fitFDist}}, \code{\link{tmixture.matrix}}.

An overview of linear model functions in limma is given by \link{06.LinearModels}.
}
\author{Gordon Smyth}
\references{
Loennstedt, I. and Speed, T. P. (2002). Replicated microarray data. \emph{Statistica Sinica} \bold{12}, 31-46.

Smyth, G. K. (2004). Linear models and empirical Bayes methods for assessing differential expression in microarray experiments.
\emph{Statistical Applications in Genetics and Molecular Biology}, \bold{3}, No. 1, Article 3. \url{"http://www.bepress.com/sagmb/vol3/iss1/art3"}
}
\examples{
#  See also lmFit examples

#  Simulate gene expression data,
#  6 microarrays and 100 genes with one gene differentially expressed
set.seed(2004); invisible(runif(100))
M <- matrix(rnorm(100*6,sd=0.3),100,6)
M[1,] <- M[1,] + 1
fit <- lmFit(M)

#  Ordinary t-statistic
par(mfrow=c(1,2))
ordinary.t <- fit$coef / fit$stdev.unscaled / fit$sigma
qqt(ordinary.t,df=fit$df.residual,main="Ordinary t")
abline(0,1)

#  Moderated t-statistic
eb <- eBayes(fit)
qqt(eb$t,df=eb$df.prior+eb$df.residual,main="Moderated t")
abline(0,1)
#  Points off the line may be differentially expressed
par(mfrow=c(1,1))
}
\keyword{htest}
