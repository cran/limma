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
  \item{object}{object of class \code{list}, \code{RGList} or \code{MAList} containing two-color microarray data}
  \item{layout}{list specifying the dimensions of the spot matrix and the grid matrix. For details see \code{\link[limma:PrintLayout]{PrintLayout-class}}.}
  \item{method}{character string specifying the normalization method.
  Choices are \code{"none"}, \code{"median"}, \code{"loess"}, \code{"printtiploess"}, \code{"composite"} and \code{"robustspline"}.
  A partial string sufficient to uniquely identify the choice is permitted.}
  \item{weights}{numeric matrix or vector of the same size and shape as the components of \code{object}. Will use by default weights found in \code{object} if they exist.}
  \item{span}{numeric scalar giving the smoothing parameter for the \code{loess} fit}
  \item{iterations}{number of iterations used in loess fitting.  More iterations give a more robust fit.}
  \item{controlspots}{numeric or logical vector specifying the subset of spots which are non-differentially expressed control spots, for use with \code{method="composite"}}
  \item{df}{degrees of freedom for spline if \code{method="robustspline"}}
  \item{robust}{robust regression method if \code{method="robustspline"}.  Choices are \code{"M"} or \code{"MM"}.}
  \item{bc.method}{character string specifying background correct method, see \code{\link{backgroundCorrect}} for options}
  \item{offset}{numeric value, intensity offset used when computing log-ratios, see \code{\link{backgroundCorrect}}}
}

\details{
Normalization is intended to remove from the expression measures any systematic trends which arise from the microarray technology rather than from differences between the probes or between the target RNA samples hybridized to the arrays.

This function normalizes M-values (log-ratios) for dye-bias within each array.
Apart from \code{method="none"} and \code{method="median"}, all the normalization methods make use of the relationship between dye-bias and intensity.
The loess normalization methods were proposed by Yang et al (2001, 2002).
Smyth and Speed (2003) give a detailed statement of the methods.

More information on the loess control parameters \code{span} and \code{iterations} can be found under \code{\link{loessFit}}.
The default values given here are equivalent to those for the older function \code{stat.ma} in the SMA package.

The \code{"robustspline"} method calls \code{\link{normalizeRobustSpline}}.

\code{MA.RG} converts an unlogged \code{RGList} object into an \code{MAList} object.
\code{MA.RG(object)} is equivalent to \code{normalizeWithinArrays(object,method="none")}.

\code{RG.MA(object)} converts back from an \code{MAList} object to a \code{RGList} object with intensities on the log2 scale.
}
\value{
An object of class \code{\link[limma:malist]{MAList}}.
Any components found in \code{object} will preserved except for \code{R}, \code{G}, \code{Rb}, \code{Gb} and \code{other}.
}

\author{Gordon Smyth}

\references{
Yang, Y. H., Dudoit, S., Luu, P., and Speed, T. P. (2001). Normalization for cDNA microarray data. In \emph{Microarrays: Optical Technologies and Informatics}, M. L. Bittner, Y. Chen, A. N. Dorsel, and E. R. Dougherty (eds), Proceedings of SPIE, Vol. 4266, pp. 141-152. 

Yang, Y. H., Dudoit, S., Luu, P., Lin, D. M., Peng, V., Ngai, J., and Speed, T. P. (2002). Normalization for cDNA microarray data: a robust composite method addressing single and multiple slide systematic variation. \emph{Nucleic Acids Research} \bold{30}(4):e15.

Smyth, G. K., and Speed, T. P. (2003). Normalization of cDNA microarray data. In: \emph{METHODS: Selecting Candidate Genes from DNA Array Screens: Application to Neuroscience}, D. Carter (ed.). Methods Volume 31, Issue 4, December 2003, pages 265-273.
}

\seealso{
An overview of LIMMA functions for normalization is given in \link{4.Normalization}.

See also \code{\link{normalizeBetweenArrays}} and \code{\link[marrayNorm]{maNorm}} in the marrayNorm package.
}

\examples{
#  See normalizeBetweenArrays
}

\keyword{models}
