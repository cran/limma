\name{dupcor.series}
\alias{duplicateCorrelation}
\alias{dupcor.series}
\title{Correlation Between Duplicates}
\description{
Estimate the correlation between duplicate spots (regularly spaced replicate spots on the same array) or between technical replicates from a series of arrays.
}
\usage{
duplicateCorrelation(object,design=rep(1,ncol(M)),ndups=2,spacing=1,block=NULL,trim=0.15,weights=NULL)
dupcor.series(M,design=rep(1,ncol(M)),ndups=2,spacing=1,initial=0.8,trim=0.15,weights=NULL)
}
\arguments{
  \item{object}{a numeric matrix of log-ratios or an \code{\link[limma:malist]{MAList}} object from which the log-ratios can be extracted.
  If \code{object} is an \code{MAList} then the arguments \code{design}, \code{ndups}, \code{spacing} and \code{weights} will be extracted from it if available and do not have to be specified as arguments.}
  \item{M}{a numeric matrix. Usually the log-ratios of expression for a series of cDNA microarrrays with rows corresponding to genes and columns to arrays.}
  \item{design}{the design matrix of the microarray experiment, with rows corresponding to arrays and columns to comparisons to be estimated. The number of rows must match the number of columns of \code{M}. Defaults to the unit vector meaning that the arrays are treated as replicates.} 
  \item{ndups}{a positive integer giving the number of times each gene is printed on an array. \code{nrow(M)} must be divisible by \code{ndups}.
  Will be ignored if \code{block} is specified.}
  \item{spacing}{the spacing between the rows of \code{M} corresponding to duplicate spots, \code{spacing=1} for consecutive spots}
  \item{block}{vector or factor specifying a blocking variable}
  \item{initial}{a numeric value between -1 and 1 giving an initial estimate for the correlation.
  Not currently used.}
  \item{trim}{the fraction of observations to be trimmed from each end of \code{tanh(all.correlations)} when computing the trimmed mean.}
  \item{weights}{an optional numeric matrix of the same dimension as \code{M} containing weights for each spot. If smaller than \code{M} then it will be filled out the same size.}
}
\value{
  A list with components
  \item{consensus.correlation}{the average estimated inter-duplicate correlation. The average is the 0.1 trimmed mean of the correlations for individual genes on the tanh-transformed scale.}
  \item{cor}{same as \code{consensus.correlation}, for compatibility with earlier versions of the software}
  \item{all.correlations}{a numeric vector of length \code{nrow(M)/ndups} giving the individual genewise correlations.}
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

The function may take long time to execute as it fits a mixed linear model for each gene.

\code{dupcor.series} produces the same values as \code{duplicateCorrelation} and it retained for compatibility with earlier releases of the software.
}
\seealso{
These functions use \code{\link[statmod]{randomizedBlockFit}} from the statmod package.

An overview of linear model functions in limma is given by \link{5.LinearModels}.
}
\author{Gordon Smyth}
\references{
Smyth, G. K., Michaud, J., and Scott, H. (2003). The use of within-array duplicate spots for assessing differential expression in microarray experiments.
\url{http://www.statsci.org/smyth/pubs/dupcor.pdf}
}
\examples{
#  See gls.series for an example
}
\keyword{multivariate}
