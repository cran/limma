\name{normalizeWithinArrays}
\alias{normalizeWithinArrays}
\alias{MA.RG}
\alias{RG.MA}
\title{Normalize Within Arrays}
\description{
Normalize the expression log-ratios for one or more two-colour spotted microarray experiments so that the log-ratios average to zero within each array or sub-array.
}
\usage{
normalizeWithinArrays(object, layout, method="printtiploess", weights=object$weights, span=0.3, iterations=4, controlspots=NULL, df=5, robust="M", bc.method="subtract", offset=0)
MA.RG(object, bc.method="subtract", offset=0)
RG.MA(object)
}
\arguments{
  \item{object}{object of class \code{list}, \code{RGList} or \code{MAList} containing red and green intensities constituting two-color microarray data.}
  \item{layout}{list specifying the dimensions of the spot matrix and the grid matrix. For details see \code{\link[limma:PrintLayout]{PrintLayout-class}}.}
  \item{method}{character string specifying the normalization method.
  Choices are \code{"none"}, \code{"median"}, \code{"loess"}, \code{"printtiploess"}, \code{"composite"}, \code{"control"} and \code{"robustspline"}.
  A partial string sufficient to uniquely identify the choice is permitted.}
  \item{weights}{numeric matrix or vector of the same size and shape as the components of \code{object} containing spot quality weights.}
  \item{span}{numeric scalar giving the smoothing parameter for the \code{loess} fit}
  \item{iterations}{number of iterations used in loess fitting.  More iterations give a more robust fit.}
  \item{controlspots}{numeric or logical vector specifying the subset of spots which are non-differentially-expressed control spots, for use with \code{method="composite"} or \code{method="control"}.}
  \item{df}{degrees of freedom for spline if \code{method="robustspline"}.}
  \item{robust}{robust regression method if \code{method="robustspline"}.  Choices are \code{"M"} or \code{"MM"}.}
  \item{bc.method}{character string specifying background correct method, see \code{\link{backgroundCorrect}} for options.}
  \item{offset}{numeric value, intensity offset used when computing log-ratios, see \code{\link{backgroundCorrect}}.}
}

\details{
Normalization is intended to remove from the expression measures any systematic trends which arise from the microarray technology rather than from differences between the probes or between the target RNA samples hybridized to the arrays.

This function normalizes M-values (log-ratios) for dye-bias within each array.
Apart from \code{method="none"} and \code{method="median"}, all the normalization methods make use of the relationship between dye-bias and intensity.
Method \code{"none"} computes M-values and A-values but does no normalization.
Method \code{"median"} subtracts the weighted median from the M-values for each array.

The loess normalization methods (\code{"loess"}, \code{"printtiploess"} and \code{"composite"}) were proposed by Yang et al (2001, 2002).
Smyth and Speed (2003) review these methods and describe how the methods are implemented in the limma package, including choices of tuning parameters.
More information on the loess control parameters \code{span} and \code{iterations} can be found under \code{\link{loessFit}}.
The default values used here are equivalent to those for the older function \code{stat.ma} in the sma package.

The \code{"control"} method fits a global loess curve through a set of control spots, such as a whole-library titration series, and applies that curve to all the other spots.

The \code{"robustspline"} method calls \code{\link{normalizeRobustSpline}}.
See that function for more documentation.

\code{MA.RG} converts an unlogged \code{RGList} object into an \code{MAList} object.
\code{MA.RG(object)} is equivalent to \code{normalizeWithinArrays(object,method="none")}.

\code{RG.MA(object)} converts back from an \code{MAList} object to a \code{RGList} object with unlogged intensities.

\code{weights} is normally a matrix giving a quality weight for every spot on every array.
If \code{weights} is instead a vector or a matrix with only one column, then the weights will be assumed to be the same for every array, i.e., the weights will be probe-specific rather than spot-specific.
}
\value{
An object of class \code{\link[limma:malist]{MAList}}.
Any components found in \code{object} will preserved except for \code{R}, \code{G}, \code{Rb}, \code{Gb} and \code{other}.
}

\author{Gordon Smyth}

\references{
Yang, Y. H., Dudoit, S., Luu, P., and Speed, T. P. (2001). Normalization for cDNA microarray data. In \emph{Microarrays: Optical Technologies and Informatics}, M. L. Bittner, Y. Chen, A. N. Dorsel, and E. R. Dougherty (eds), Proceedings of SPIE, Vol. 4266, pp. 141-152. 

Yang, Y. H., Dudoit, S., Luu, P., Lin, D. M., Peng, V., Ngai, J., and Speed, T. P. (2002). Normalization for cDNA microarray data: a robust composite method addressing single and multiple slide systematic variation. \emph{Nucleic Acids Research} \bold{30}(4):e15.

Smyth, G. K., and Speed, T. P. (2003). Normalization of cDNA microarray data. \emph{Methods} \bold{31}, 265-273.
}

\seealso{
An overview of limma functions for normalization is given in \link{05.Normalization}.
In particular, see \code{\link{normalizeBetweenArrays}} for between-array normalization.

The original loess normalization function was the \code{\link[sma]{statma}} funtion in the sma package.
\code{normalizeWithinArrays} is a direct generalization of that function, with more options and with support for quantitative spot quality weights.

A different implementation of loess normalization methods is provided by the \code{\link[marray]{maNorm}} in the marray package.

An alternative wrapper for loess normalization, using different data classes but calling the limma functions, is provided by the \code{\link[arrayMagic]{normalise}} function in the arrayMagic package.
}

\examples{
#  See normalizeBetweenArrays
}

\keyword{models}
