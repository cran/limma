\name{arrayWeights}
\alias{arrayWeights}
\title{Array Quality Weights}
\description{
  Estimates relative quality weights for each array in a multi-array
  experiment with replication.
}
\usage{
arrayWeights(object, design = NULL, weights = NULL, method = "genebygene", maxiter = 50, tol = 1e-15, trace=FALSE)
}
\arguments{
 \item{object}{object of class \code{numeric}, \code{matrix}, \code{MAList}, \code{marrayNorm},
           \code{exprSet} or \code{PLMset} containing log-ratios or log-values of
            expression for a series of microarrays.}
 \item{design}{the design matrix of the microarray experiment, with rows
           corresponding to arrays and columns to coefficients to be
           estimated.  Defaults to the unit vector meaning that the
           arrays are treated as replicates.}
 \item{weights}{optional numeric matrix containing prior weights for each spot.}
 \item{method}{character string specifying the estimating algorithm to be used. Choices
          are \code{"genebygene"} and \code{"reml"}.}
 \item{maxiter}{maximum number of iterations allowed.}
 \item{tol}{convergence tolerance.}
 \item{trace}{logical variable. If true then output diagnostic information
          at each iteration of '"reml"' algorithm.}
  }
\details{
     The relative reliability of each array is estimated by measuring how
     well the expression values for that array follow the linear model.

     A heteroscedastic model is fitted to the expression values for 
     each gene by calling the function \code{lm.wfit}.  The dispersion model 
     is fitted to the squared residuals from the mean fit, and is set up to 
     have array specific coefficients, which are updated in either full REML 
     scoring iterations, or using an efficient gene-by-gene update algorithm.  
     The final estimates of these array variances are converted to weights.

     The arguments \code{design} and \code{weights} will be extracted from the data 
     \code{object} if available and do not normally need to be set explicitly in 
     the call; if any of these are set in the call then they will over-ride 
     the slots or components in the data \code{object}.

     If \code{object} is a \code{PLMset}, then expression values will be taken
     from the slot \code{chip.coefs} and weights will be computed from
     \code{se.chip.coefs}. If \code{object} is an \code{exprSet}, then expression
     values will be taken from the \code{exprs} slot, but weights will not
     be computed.
}
\value{
	A matrix of array weights, suitable for use in the \code{weights} argument of \code{lmFit}.
 }
\seealso{
An overview of linear model functions in limma is given by \link{06.LinearModels}.
}
\examples{
\dontrun{
array.wts <- arrayWeights(MA, design)
fit.wts <- lmFit(MA, design, weights=array.wts)}
}

\author{Matthew Ritchie}
\keyword{regression}
\keyword{models}
