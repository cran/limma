\name{dupcor}
\alias{duplicateCorrelation}
\title{Correlation Between Duplicates}
\description{
Estimate the correlation between duplicate spots (regularly spaced replicate spots on the same array) or between technical replicates from a series of arrays.
}
\usage{
duplicateCorrelation(object,design=rep(1,ncol(as.matrix(object))),ndups=2,spacing=1,block=NULL,trim=0.15,weights=NULL)
}
\arguments{
  \item{object}{a numeric matrix of expression values, or any data object from which \code{as.matrix} will extract a suitable matrix such as an \code{\link[limma:malist]{MAList}}, \code{marrayNorm} or \code{exprSet} object.
  If \code{object} is an \code{MAList} object then the arguments \code{design}, \code{ndups}, \code{spacing} and \code{weights} will be extracted from it if available and do not have to be specified as arguments.
  Specifying these arguments explicitly will over-rule any components found in the data object.}
  \item{design}{the design matrix of the microarray experiment, with rows corresponding to arrays and columns to comparisons to be estimated. The number of rows must match the number of columns of \code{object}. Defaults to the unit vector meaning that the arrays are treated as replicates.} 
  \item{ndups}{a positive integer giving the number of times each gene is printed on an array. \code{nrow(object)} must be divisible by \code{ndups}.
  Will be ignored if \code{block} is specified.}
  \item{spacing}{the spacing between the rows of \code{object} corresponding to duplicate spots, \code{spacing=1} for consecutive spots}
  \item{block}{vector or factor specifying a blocking variable}
  \item{trim}{the fraction of observations to be trimmed from each end of \code{tanh(all.correlations)} when computing the trimmed mean.}
  \item{weights}{an optional numeric matrix of the same dimension as \code{object} containing weights for each spot. If smaller than \code{object} then it will be filled out the same size.}
}
\value{
  A list with components
  \item{consensus.correlation}{the average estimated inter-duplicate correlation. The average is the trimmed mean of the individual correlations on the atanh-transformed scale.}
  \item{cor}{same as \code{consensus.correlation}, for compatibility with earlier versions of the software}
  \item{atanh.correlations}{numeric vector of length \code{nrow(object)/ndups} giving the individual genewise atanh-transformed correlations.}
}
\details{
When \code{block=NULL}, this function estimates the correlation between duplicate spots (regularly spaced within-array replicate spots).
If \code{block} is not null, this function estimates the correlation between repeated observations on the blocking variable.
Typically the blocks are biological replicates and the repeated observations are technical replicates.
In either case, the correlation is estimated by fitting a mixed linear model by REML individually for each gene.
The function also returns a consensus correlation, which is a robust average of the individual correlations, which can be used as input for 
functions \code{lmFit} or \code{gls.series}.

At this time it is not possible to estimate correlations between duplicate spots and between technical replicates simultaneously.
If \code{block} is not null, then the function will set \code{ndups=1}.

For this function to return statistically useful results, there must be at least two more arrays than the number of coefficients to be estimated, i.e., two more than the column rank of \code{design}.

The function may take long time to execute as it fits a mixed linear model for each gene for an iterative algorithm.
It is not uncommon for the function to return a small number of warning messages that correlation estimates cannot be computed for some individual genes.
This is not a serious concern providing that there are only a few such warnings and the total number of genes is large.
The consensus estimator computed by this function will not be materially affected by a small number of genes.
}
\seealso{
These functions use \code{\link[statmod]{randomizedBlockFit}} from the statmod package.

An overview of linear model functions in limma is given by \link{06.LinearModels}.
}
\author{Gordon Smyth}
\references{
Smyth, G. K., Michaud, J., and Scott, H. (2005). The use of within-array replicate spots for assessing differential expression in microarray experiments. \emph{Bioinformatics} 21(9), 2067-2075.
\url{http://www.statsci.org/smyth/pubs/dupcor.pdf}
}
\examples{
#  Also see lmFit examples

\dontrun{
corfit <- duplicateCorrelation(MA, ndups=2, design)
all.correlations <- tanh(corfit$atanh.correlations)
boxplot(all.correlations)
fit <- lmFit(MA, design, ndups=2, correlation=corfit$consensus)
}
}
\keyword{multivariate}
