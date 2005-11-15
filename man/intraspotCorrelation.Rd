\name{intraspotCorrelation}
\alias{intraspotCorrelation}
\title{Intra-Spot Correlation for Two Color Data}
\description{
Estimate the within-block correlation associated with spots for spotted two color microarray data.
}
\usage{
intraspotCorrelation(object, design, trim=0.15)
}
\arguments{
  \item{object}{an \code{\link[limma:malist]{MAList}} object or a list from which \code{M} and \code{A} values may be extracted}
  \item{design}{a numeric matrix containing the design matrix for linear model in terms of the individual channels.
  The number of rows should be twice the number of arrays.
  The number of columns will determine the number of coefficients estimated for each gene.}
  \item{trim}{the fraction of observations to be trimmed from each end of the atanh-correlations when computing the consensus correlation. See \code{\link{mean}}.}
}
\value{
  A list with components
  \item{consensus.correlation}{robust average of the estimated inter-duplicate correlations. The average is the trimmed mean of the correlations for individual genes on the atanh-transformed scale.}
  \item{atanh.correlations}{a numeric vector giving the individual genewise correlations on the atanh scale}
  \item{df}{numeric matrix of degrees of freedom associated with the correlations. The first column gives the degrees of freedom for estimating the within-spot or M-value mean square while the second gives the degrees of freedom for estimating the between spot or A-value mean square.}
}
\details{
This function estimates the correlation between two channels observed on each spot.
The correlation is estimated by fitting a heteroscedastic regression model to the M and A-values of each gene.
The function also returns a consensus correlation, which is a robust average of the individual correlations, which can be used as input for 
functions \code{lmscFit}.

The function may take long time to execute.
}
\seealso{
This function uses \code{\link[statmod]{remlscore}} from the statmod package.

An overview of methods for single channel analysis in limma is given by \link{07.SingleChannel}.
}
\examples{
#  See lmscFit
\dontrun{
corfit <- dupclicateCorrelation(MA, design)
all.correlations <- tanh(corfit$atanh.correlations)
boxplot(all.correlations)
}
}
\author{Gordon Smyth}
\keyword{multivariate}
