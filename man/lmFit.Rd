\name{lmFit}
\alias{lmFit}
\title{Linear Model for Series of Arrays}
\description{Fit linear model for each gene given a series of arrays}
\usage{
lmFit(object,design=NULL,ndups=1,spacing=1,block=NULL,correlation=0.75,weights=NULL,method="ls",...) 
}
\arguments{
  \item{object}{object of class \code{numeric}, \code{matrix}, \code{MAList}, \code{marrayNorm}, \code{exprSet} or \code{PLMset} containing log-ratios or log-values of expression for a series of microarrays}
  \item{design}{the design matrix of the microarray experiment, with rows corresponding to arrays and columns to coefficients to be estimated.  Defaults to the unit vector meaning that the arrays are treated as replicates.} 
  \item{ndups}{positive integer giving the number of times each gene is printed on an array}
  \item{spacing}{positive integer giving the spacing between duplicate spots, \code{spacing=1} for consecutive spots}
  \item{block}{vector or factor specifying a blocking variable on the arrays. Has length equal to the number of arrays.}
  \item{correlation}{the inter-duplicate or inter-technical replicate correlation}
  \item{weights}{optional numeric matrix containing weights for each spot}
  \item{method}{character string, \code{"ls"} for least squares or \code{"robust"} for robust regression}
  \item{...}{other optional arguments to be passed to \code{lm.series}, \code{gls.series} or \code{rlm.series}}
}

\value{
Object of class \code{\link[limma:marraylm]{MArrayLM}}
}

\details{
This function accepts data from a experiment involving a series of microarrays with the same set of probes.
A linear model is fitted to the expression data for each probe.
The expression data should be log-ratios for two-color array platforms or log-expression values for one-channel platforms.
(To fit linear models to the individual channels of two-color array data, see \code{\link{lmscFit}}.)
The coefficients of the fitted models describe the differences between the RNA sources hybridized to the arrays.
The genewise fitted model results are stored in a compact form suitable for further processing by other functions in the limma package.

If \code{object} is a matrix then it should contain log-ratios or log-expression data with rows corresponding to probes and columns to arrays.
(A numeric vector is treated the same as a matrix with one column.)
For objects of other classes, a matrix of expression values is taken from the appropriate component or slot of the object.
If \code{object} is of class \code{MAList} or \code{marrayNorm}, then the matrix of log-ratios (M-values) is extracted.
If \code{object} is of class \code{exprSet} then the \code{exprs} slot is extracted.
(This may contain log-expression or log-ratio values, depending on the platform.)
If \code{object} is of class \code{PLMset} then the matrix of chip coefficients \code{chip.coefs} is extracted.

The arguments \code{design}, \code{ndups}, \code{spacing} and \code{weights} will be extracted from the data \code{object} if available and do not normally need to set explicitly in the call.
On the other hand, if any of these are set in the function call then they will over-ride the slots or components in the data \code{object}.
If \code{object} is an \code{PLMset}, then weights are computed as \code{1/pmax(object@se.chip.coefs, 1e-05)^2}.
If \code{object} is an \code{exprSet} object, then weights are not computed.

The actual linear model computations are done by passing the data to one the lower-level functions \code{lm.series}, \code{gls.series} or \code{rlm.series}.
The function \code{rlm.series} is used if \code{method="robust"}.
If \code{method="ls"}, then \code{gls.series} is used if a correlation structure has been specified, i.e., if \code{ndups>1} or \code{block} is non-null and \code{correlation} is different from zero.
If \code{method="ls"} and there is no correlation structure, \code{lm.series} is used.
}

\seealso{
An overview of linear model functions in limma is given by \link{5.LinearModels}.
}

\author{Gordon Smyth}
\keyword{models}
\keyword{regression}
