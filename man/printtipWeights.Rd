\name{printtipWeights}
\alias{printtipWeights}
\alias{printtipWeightsSimple}
\title{Sub-array Quality Weights}
\description{
Estimates relative quality weights for each sub-array in a multi-array
experiment.
}
\usage{
printtipWeights(object, design = NULL, weights = NULL, method = "genebygene", layout, maxiter = 50, tol = 1e-10, trace=FALSE)

}
\arguments{
 \item{object}{object of class \code{numeric}, \code{matrix}, \code{MAList}, \code{marrayNorm},
           or \code{exprSet} containing log-ratios or log-values of
           expression for a series of spotted microarrays.}
 \item{design}{the design matrix of the microarray experiment, with rows
           corresponding to arrays and columns to coefficients to be
           estimated.  Defaults to the unit vector meaning that the
           arrays are treated as replicates.}
 \item{weights}{optional numeric matrix containing prior weights for each spot.}
 \item{method}{character string specifying the estimating algorithm to be used. Choices
          are \code{"genebygene"} and \code{"reml"}.}
 \item{layout}{list specifying the dimensions of the spot matrix and the grid matrix. For details see \code{\link[limma:PrintLayout]{PrintLayout-class}}.}
 \item{maxiter}{maximum number of iterations allowed.}
 \item{tol}{convergence tolerance.}
 \item{trace}{logical variable. If true then output diagnostic information
          at each iteration of '"reml"' algorithm.}
  }
\details{
The relative reliability of each sub-array (print-tip group) is estimated by measuring how
well the expression values for that sub-array follow the linear model.

The method described in Ritchie et al (2006) and implemented in 
the \code{arrayWeights} function is adapted for this purpose.A heteroscedastic model is fitted to the expression values for 
each gene by calling the function \code{lm.wfit}.  The dispersion model 
is fitted to the squared residuals from the mean fit, and is set up to 
have sub-array specific coefficients, which are updated in either full REML 
scoring iterations, or using an efficient gene-by-gene update algorithm.  
The final estimates of the sub-array variances are converted to weights.

The data object \code{object} is interpreted as for \code{lmFit}.
In particular, the arguments \code{design}, \code{weights} and \code{layout} will 
be extracted from the data \code{object} if available and do not normally need to 
be set explicitly in the call; if any of these are set in the call then they 
will over-ride the slots or components in the data \code{object}.

}
\value{
	A matrix of sub-array weights which can be passed to \code{lmFit}.
 }
\references{
Ritchie, M. E., Diyagama, D., Neilson, van Laar, R., J., Dobrovic, A., Holloway, A., and Smyth, G. K. (2006). Empirical array quality weights in the analysis of microarray data. BMC Bioinformatics 7, 261. \url{http://www.biomedcentral.com/1471-2105/7/261/abstract}
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
subarrayw <- printtipWeights(MA, design, layout=mouse.setup)
fit <- lmFit(MA, design, weights=subarrayw)
fit2 <- contrasts.fit(fit, contrasts=c(-1,1))
fit2 <- eBayes(fit2)
# Use of sub-array weights increases the significance of the top genes
topTable(fit2)
# Create an image plot of sub-array weights from each array
zlim <- c(min(subarrayw), max(subarrayw))
par(mfrow=c(3,2), mai=c(0.1,0.1,0.3,0.1))
for(i in 1:6) 
	imageplot(subarrayw[,i], layout=mouse.setup, zlim=zlim, main=paste("Array", i))
}

\author{Matthew Ritchie and Gordon Smyth}
\keyword{regression}
\keyword{models}
