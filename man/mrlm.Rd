\name{rlm}
\alias{rlm.series}
\alias{mrlm}
\title{Robust Linear Model for Series of Microarrays}
\description{Fit linear models for each gene to a series of microarrays. Fit is by robust M-estimation.}
\usage{
mrlm(M,design=NULL,ndups=1,spacing=1,weights=NULL,...)
rlm.series(x,...)
}
\arguments{
  \item{M}{numeric matrix containing log-ratio or log-expression values for a series of microarrays, rows correspond to genes and columns to arrays.}
  \item{x}{same as \code{M}}
  \item{design}{the design matrix of the microarray experiment, with rows corresponding to arrays and columns to comparisons to be estimated. The number of rows must match the number of columns of \code{M}. Defaults to the unit vector meaning that the arrays are treated as replicates.} 
  \item{ndups}{a positive integer giving the number of times each gene is printed on an array. \code{nrow(M)} must be divisible by \code{ndups}.}
  \item{spacing}{the spacing between the rows of \code{M} corresponding to duplicate spots, \code{spacing=1} for consecutive spots.}
  \item{weights}{numeric matrix of the same dimension as \code{M} containing weights. If it is of different dimension to \code{M}, it will be filled out to the same size. \code{NULL} is equivalent to equal weights.}
  \item{...}{any other arguments are passed to \code{rlm.default}.}
}
\value{
  A list with components
  \item{coefficients}{numeric matrix containing the estimated coefficients for each linear model. Same number of rows as \code{M}, same number of columns as \code{design}.}
  \item{stdev.unscaled}{numeric matrix conformal with \code{coef} containing the unscaled standard deviations for the coefficient estimators. The standard errors are given by \code{stdev.unscaled * sigma}.}
  \item{sigma}{numeric vector containing the residual standard deviation for each gene.}
  \item{df.residual}{numeric vector giving the degrees of freedom corresponding to \code{sigma}.}
  \item{qr}{QR decomposition of \code{design}.}
}
\details{
The linear model is fit for each gene by calling the function \code{rlm} from the MASS library.

Warning: don't use weights with this function unless you understand how \code{rlm} treats weights.
The treatment of weights is somewhat different from that of \code{lm.series} and \code{gls.series}.

The function \code{rlm.series} is equivalent to \code{mrlm} but is deprecated and will be removed at some time in the future.
}
\seealso{
\code{\link[MASS]{rlm}}.

An overview of linear model functions in limma is given by \link{06.LinearModels}.
}
\author{Gordon Smyth}
\examples{
#  Simulate gene expression data,
#  6 microarrays and 100 genes with one gene differentially expressed in first 3 arrays
M <- matrix(rnorm(100*6,sd=0.3),100,6)
M[1,1:3] <- M[1,1:3] + 2
#  Design matrix includes two treatments, one for first 3 and one for last 3 arrays
design <- cbind(First3Arrays=c(1,1,1,0,0,0),Last3Arrays=c(0,0,0,1,1,1))
fit <- mrlm(M,design=design)
eb <- ebayes(fit)
#  Large values of eb$t indicate differential expression
qqt(eb$t[,1],df=fit$df+eb$df.prior)
abline(0,1)
}
\keyword{models}
\keyword{regression}
