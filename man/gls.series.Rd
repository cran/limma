\name{gls.series}
\alias{gls.series}
\title{Fit Linear Model to Microarray Data by Generalized Least Squares}
\description{
Fit a linear model genewise to expression data from a series of microarrays.
The fit is by generalized least squares allowing for correlation between duplicate spots or related arrays.
This is a utility function for \code{lmFit}.
}
\usage{gls.series(M,design=NULL,ndups=2,spacing=1,block=NULL,correlation=NULL,weights=NULL,...)}
\arguments{
  \item{M}{numeric matrix containing log-ratio or log-expression values for a series of microarrays, rows correspond to genes and columns to arrays.}
  \item{design}{numeric design matrix defining the linear model, with rows corresponding to arrays and columns to comparisons to be estimated. The number of rows must match the number of columns of \code{M}. Defaults to the unit vector meaning that the arrays are treated as replicates.} 
  \item{ndups}{positive integer giving the number of times each gene is printed on an array. \code{nrow(M)} must be divisible by \code{ndups}.}
  \item{spacing}{the spacing between the rows of \code{M} corresponding to duplicate spots, \code{spacing=1} for consecutive spots}
  \item{block}{vector or factor specifying a blocking variable on the arrays.
  Same length as \code{ncol(M)}.}
  \item{correlation}{numeric value specifying the inter-duplicate or inter-block correlation.}
  \item{weights}{an optional numeric matrix of the same dimension as \code{M} containing weights for each spot. If it is of different dimension to \code{M}, it will be filled out to the same size.}
  \item{...}{other optional arguments to be passed to \code{dupcor.series}.}
}
\value{
  A list with components
  \item{coefficients}{numeric matrix containing the estimated coefficients for each linear model. Same number of rows as \code{M}, same number of columns as \code{design}.}
  \item{stdev.unscaled}{numeric matrix conformal with \code{coef} containing the unscaled standard deviations for the coefficient estimators. The standard errors are given by \code{stdev.unscaled * sigma}.}
  \item{sigma}{numeric vector containing the residual standard deviation for each gene.}
  \item{df.residual}{numeric vector giving the degrees of freedom corresponding to \code{sigma}}
  \item{correlation}{inter-duplicate or inter-block correlation}
  \item{qr}{QR decomposition of the generalized linear squares problem, i.e., the decomposition of \code{design} standardized by the Choleski-root of the correlation matrix defined by \code{correlation}}
}
\details{
This is a utility function used by the higher level function \code{\link{lmFit}}.
Most users should not use this function directly but should use \code{\link{lmFit}} instead.

This function is for fitting gene-wise linear models when some of the expression values are correlated.
The correlated groups may arise from replicate spots on the same array (duplicate spots) or from a biological or technical replicate grouping of the arrays.
This function is normally called by \code{lmFit} and is not normally called directly by users.

Note that the correlation is assumed to be constant across genes.
If \code{correlation=NULL} then a call is made to \code{duplicateCorrelation} to estimated the correlation.
}
\seealso{
\code{\link{duplicateCorrelation}}.

An overview of linear model functions in limma is given by \link{06.LinearModels}.
}
\author{Gordon Smyth}
\keyword{models}
\keyword{regression}

