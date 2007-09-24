\name{backgroundCorrect}
\alias{backgroundCorrect}
\title{Correct Intensities for Background}
\description{
Background correct microarray expression intensities.
}
\usage{
backgroundCorrect(RG, method="subtract", offset=0, printer=RG$printer, verbose=TRUE)
}
\arguments{
  \item{RG}{an \code{\link[limma:rglist]{RGList}} object, an unclassed list containing the same components as an \code{RGList}, or a numeric matrix.}
  \item{method}{character string specifying correction method.  Possible values are \code{"none"}, \code{"subtract"}, \code{"half"}, \code{"minimum"}, \code{"movingmin"}, \code{"edwards"}, \code{"normexp"} or \code{"rma"}.
  If \code{RG} is a matrix, possible values are restricted to \code{"none"}, \code{"normexp"} or \code{"rma"}.}
  \item{offset}{numeric value to add to intensities}
  \item{printer}{a list containing printer layout information, see \code{\link{PrintLayout-class}}}
  \item{verbose}{logical, should progress be reported to standard output}
}
\details{
If \code{method="none"} then the corrected intensities are equal to the foreground intensities, i.e., the background intensities are treated as zero.
If \code{method="subtract"} then this function simply subtracts the background intensities from the foreground intensities which is the usual background correction method.
If \code{method="movingmin"} then the background estimates are replaced with the minimums of the backgrounds of the spot and its eight neighbors, i.e., the background is replaced by a moving minimum of 3x3 grids of spots.

The remaining methods are all designed to produce positive corrected intensities.
If \code{method="half"} then any intensity which is less than 0.5 after background subtraction is reset to be equal to 0.5.
If \code{method="minimum"} then any intensity which is zero or negative after background subtraction is set equal to half the minimum of the positive corrected intensities for that array.
If \code{method="edwards"} a log-linear interpolation method is used to adjust lower intensities as in Edwards (2003).
If \code{method="normexp"} a convolution of normal and exponential distributions is fitted to the foreground intensities using the background intensities as a covariate, and the expected signal given the observed foreground becomes the corrected intensity.
This results in a smooth monotonic transformation of the background subtracted intensities such that all the corrected intensities are positive.
See Ritchie et al (2007) and \code{\link{normexp.fit}} for more details.

The \code{offset} can be used to add a constant to the intensities before log-transforming, so that the log-ratios are shrunk towards zero at the lower intensities.
This may eliminate or reverse the usual 'fanning' of log-ratios at low intensities associated with local background subtraction.

Background correction (background subtraction) is also performed by the \code{\link{normalizeWithinArrays}} method for \code{RGList} objects, so it is not necessary to call \code{backgroundCorrect} directly unless one wants to use a method other than simple subtraction.
Calling \code{backgroundCorrect} before \code{normalizeWithinArrays} will over-ride the default background correction.
}
\value{
An \code{RGList} object in which components \code{R} and \code{G} are background corrected
and components \code{Rb} and \code{Gb} are removed.
}
\references{
Edwards, D. E. (2003). Non-linear normalization and background correction in one-channel cDNA microarray studies
\emph{Bioinformatics} 19, 825-833. 

Ritchie, M. E., Silver, J., Oshlack, A., Silver, J., Holmes, M., Diyagama, D., Holloway, A., and Smyth, G. K. (2007).
A comparison of background correction methods for two-colour microarrays.
\emph{Bioinformatics}
\url{http://bioinformatics.oxfordjournals.org/cgi/content/abstract/btm412}

Yang, Y. H., Buckley, M. J., Dudoit, S., and Speed, T. P. (2002). Comparison of methods for image analysis on cDNA microarray data. \emph{Journal of Computational and Graphical Statistics} 11, 108-136.

Yang, Y. H., Buckley, M. J., and Speed, T. P. (2001). Analysis of microarray images. \emph{Briefings in Bioinformatics} 2, 341-349.
}
\author{Gordon Smyth}
\examples{
RG <- new("RGList", list(R=c(1,2,3,4),G=c(1,2,3,4),Rb=c(2,2,2,2),Gb=c(2,2,2,2)))
backgroundCorrect(RG)
backgroundCorrect(RG, method="half")
backgroundCorrect(RG, method="minimum")
backgroundCorrect(RG, offset=5)
}
\seealso{
An overview of background correction functions is given in \code{\link{04.Background}}.
}
\keyword{models}
