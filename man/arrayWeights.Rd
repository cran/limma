\name{arrayWeights}
\alias{arrayWeights}
\alias{arrayWeightsSimple}
\title{Array Quality Weights}
\description{
  Estimates relative quality weights for each array in a multi-array
  experiment with replication.
}
\usage{
arrayWeights(object, design = NULL, weights = NULL, method = "genebygene", maxiter = 50, tol = 1e-10, trace=FALSE)
arrayWeightsSimple(object, design = NULL, maxiter = 100, tol = 1e-6, maxratio = 100, trace=FALSE)
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
 \item{maxratio}{maximum ratio between largest and smallest weights before iteration stops}
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

The data object \code{object} is interpreted as for \code{lmFit}.
In particular, the arguments \code{design} and \code{weights} will be extracted from the data 
\code{object} if available and do not normally need to be set explicitly in 
the call; if any of these are set in the call then they will over-ride 
the slots or components in the data \code{object}.

\code{arrayWeightsSimple} is a fast version of \code{arrayWeights} with \code{method="reml"}, no prior weights and no missing values.
}
\value{
	A vector of array weights.
 }
\references{
Ritchie, M. E., Diyagama, D., Neilson, van Laar, R., J., Dobrovic, A., Holloway, A., and Smyth, G. K. (2006). Empirical array quality weights for microarray data.
}
\seealso{
An overview of linear model functions in limma is given by \link{06.LinearModels}.
}
\examples{
library(sma)
# Subset of data from ApoAI case study in Limma User's Guide
data(MouseArray)
# Avoid non-positive intensities
RG <- backgroundCorrect(mouse.data, method="half")
MA <- normalizeWithinArrays(RG, mouse.setup)
MA <- normalizeBetweenArrays(MA, method="Aq")
targets <- data.frame(Cy3=I(rep("Pool",6)),Cy5=I(c("WT","WT","WT","KO","KO","KO")))
design <- modelMatrix(targets, ref="Pool")
arrayw <- arrayWeightsSimple(MA, design)
fit <- lmFit(MA, design, weights=arrayw)
fit2 <- contrasts.fit(fit, contrasts=c(-1,1))
fit2 <- eBayes(fit2)
# Use of array weights increases the significance of the top genes
topTable(fit2)
}

\author{Matthew Ritchie and Gordon Smyth}
\keyword{regression}
\keyword{models}
