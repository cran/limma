\name{normalizeRobustSpline}
\alias{normalizeRobustSpline}
\title{Normalize Single Microarray Using Shrunk Robust Splines}
\description{
Normalize the M-values for a single microarray using robustly fitted regression splines and empirical Bayes shrinkage.
}
\usage{
normalizeRobustSpline(M,A,layout,df=5,method="M")
}
\arguments{
  \item{M}{numeric vector of M-values}
  \item{A}{numeric vector of A-values}
  \item{layout}{list specifying the dimensions of the spot matrix and the grid matrix}
  \item{df}{degrees of freedom for regression spline, i.e., the number of regression coefficients and the number of knots}
  \item{method}{choices are \code{"M"} for M-estimation or \code{"MM"} for high breakdown point regression}
}
\details{
This function implements an idea similar to print-tip loess normalization but uses regression splines in place of the loess curves and uses empirical Bayes ideas to shrink the individual prtin-tip curves towards a common value.
This allows the technique to introduce less noise into good quality arrays with little spatial variation while still giving good results on arrays with strong spatial variation.
}
\value{
Numeric vector containing normalized M-values.
}
\author{Gordon Smyth}
\references{
The function is based on unpublished work by the author.
}
\seealso{
  An overview of LIMMA functions for normalization is given in \link{05.Normalization}.
}
\examples{
library(sma)
data(MouseArray)
MA <- MA.RG(mouse.data)
normM <- normalizeRobustSpline(MA$M[,1],MA$A[,1],mouse.setup)
}
\keyword{models}
