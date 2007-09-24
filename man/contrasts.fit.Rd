\name{contrasts.fit}
\alias{contrasts.fit}
\title{Compute Contrasts from Linear Model Fit}
\description{
Given a linear model fit to microarray data, compute estimated coefficients and standard errors for a given set of contrasts.
}
\usage{
contrasts.fit(fit, contrasts=NULL, coefficients=NULL)
}
\arguments{
  \item{fit}{an \code{MArrayLM} object or a list object produced by the function \code{lm.series} or equivalent. Must contain components \code{coefficients} and \code{stdev.unscaled}.}
  \item{contrasts}{numeric matrix with row corresponding to coefficients in \code{fit} and columns containing contrasts. May be a vector if there is only one contrast.}
  \item{coefficients}{vector indicating which coefficients are to be kept in the revised fit object. An alternative way to specify the \code{contrasts}.}
}
\value{
  An list object of the same class as \code{fit}. This is a list with components components
  \item{coefficients}{numeric matrix containing the estimated coefficients for each contrast for each probe.}
  \item{stdev.unscaled}{numeric matrix conformal with \code{coef} containing the unscaled standard deviations for the coefficient estimators.}
  \item{cor.coefficients}{numeric correlation matrix for the estimated coefficients}
  \item{...}{any other components input in \code{fit}}
}
\details{
This function accepts input from any of the functions \code{\link{lmFit}}, \code{lm.series}, \code{mrlm}, \code{gls.series} or \code{lmscFit}.
The function re-orientates the fitted model object from the coefficients of the original design matrix to any set of contrasts of the original coefficients.
The coefficients, unscaled standard deviations and correlation matrix are re-calculated in terms of the contrasts.

The idea of this function is to fit a full-rank model using \code{lmFit} or equivalent, then use \code{contrasts.fit} to obtain coefficients and standard errors for any number of contrasts of the coefficients of the original model.
Unlike the design matrix input to \code{lmFit}, which normally has one column for each treatment in the experiment, the matrix \code{contrasts} may have any number of columns and these are not required to be linearly independent.
Methods of assessing differential expression, such as \code{eBayes} or \code{classifyTestsF}, can then be applied to fitted model object.

The \code{coefficients} argument provides a simpler way to specify the \code{contrasts} matrix when the desired contrasts are just a subset of the original coefficients.

Warning. For efficiency reasons, this function does not re-factorize the design matrix for each probe.
A consequence is that, if the design matrix is non-orthogonal and the original fit included quality weights or missing values, then the unscaled standard deviations produced by this function are approximate rather than exact.
The approximation is usually acceptable.
The results are always exact if the original fit was a oneway model.
}
\seealso{
An overview of linear model functions in limma is given by \link{06.LinearModels}.
}
\author{Gordon Smyth}
\examples{
#  Simulate gene expression data: 6 microarrays and 100 genes
#  with one gene differentially expressed in first 3 arrays
M <- matrix(rnorm(100*6,sd=0.3),100,6)
M[1,1:3] <- M[1,1:3] + 2
#  Design matrix corresponds to oneway layout, columns are orthogonal
design <- cbind(First3Arrays=c(1,1,1,0,0,0),Last3Arrays=c(0,0,0,1,1,1))
fit <- lmFit(M,design=design)
#  Would like to consider original two estimates plus difference between first 3 and last 3 arrays
contrast.matrix <- cbind(First3=c(1,0),Last3=c(0,1),"Last3-First3"=c(-1,1))
fit2 <- contrasts.fit(fit,contrast.matrix)
fit2 <- eBayes(fit2)
#  Large values of eb$t indicate differential expression
results <- decideTests(fit2,method="global",p=0.3)
summary(results)
vennCounts(results)
}
\keyword{htest}
