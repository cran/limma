\name{rlm}
\alias{rlm.series}
\alias{mrlm}
\title{Fit Linear Model to Microrray Data by Robust Regression}
\description{Fit a linear model genewise to expression data from a series of arrays.
The fit is by robust M-estimation allowing for a small proportion of outliers.
This is a utility function for \code{lmFit}.}
\usage{
mrlm(M,design=NULL,ndups=1,spacing=1,weights=NULL,...)
rlm.series(x,...)
}
\arguments{
  \item{M}{numeric matrix containing log-ratio or log-expression values for a series of microarrays, rows correspond to genes and columns to arrays.}
  \item{x}{same as \code{M}}
  \item{design}{numeric design matrix defining the linear model, with rows corresponding to arrays and columns to comparisons to be estimated. The number of rows must match the number of columns of \code{M}. Defaults to the unit vector meaning that the arrays are treated as replicates.} 
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
This is a utility function used by the higher level function \code{\link{lmFit}}.
Most users should not use this function directly but should use \code{\link{lmFit}} instead.

This function fits linear models each gene by calling the function \code{rlm} from the MASS library.

Warning: don't use weights with this function unless you understand how \code{rlm} treats weights.
The treatment of weights is somewhat different from that of \code{lm.series} and \code{gls.series}.

The function \code{rlm.series} is equivalent to \code{mrlm} but is deprecated and will be removed at some time in the future.
}
\seealso{
\code{\link[MASS]{rlm}}.

An overview of linear model functions in limma is given by \link{06.LinearModels}.
}
\author{Gordon Smyth}
\keyword{models}
\keyword{regression}
