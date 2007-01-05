\name{normalizeForPrintorder}
\alias{normalizeForPrintorder}
\alias{normalizeForPrintorder.rg}
\alias{plotPrintorder}
\title{Print-Order Normalization}

\description{
Normalize intensity values on one or more spotted microarrays to adjust for print-order effects.
}
\usage{
normalizeForPrintorder(object, layout, start="topleft", method = "loess", separate.channels = FALSE, span = 0.1, plate.size = 32)
normalizeForPrintorder.rg(R, G, printorder, method = "loess", separate.channels = FALSE, span = 0.1, plate.size = 32, plot = FALSE)
plotPrintorder(object, layout, start="topleft", slide = 1, method = "loess", separate.channels = FALSE, span = 0.1, plate.size = 32)
}
\arguments{
  \item{object}{an \code{RGList} or \code{list} object containing components \code{R} and \code{G} which are matrices containing the red and green channel intensities for a series of arrays}
  \item{R}{numeric vector containing red channel intensities for a single microarray}
  \item{G}{numeric vector containing the green channel intensities for a single microarray}
  \item{layout}{list specifying the printer layout, see \code{\link{PrintLayout-class}}}
  \item{start}{character string specifying where printing starts in each pin group.  Choices are \code{"topleft"} or \code{"topright"}.}
  \item{printorder}{numeric vector specifying order in which spots are printed.
  Can be computed from \code{printorder(layout,start=start)}.}
  \item{slide}{positive integer giving the column number of the array for which a plot is required}
  \item{method }{character string, "loess" if a smooth loess curve should be fitted through the print-order trend or "plate" if plate effects are to be estimated}
  \item{separate.channels}{logical, \code{TRUE} if normalization should be done separately for the red and green channel and \code{FALSE} if the normalization should be proportional for the two channels}
  \item{span}{numerical constant between 0 and 1 giving the smoothing span for the loess the curve.  Ignored if \code{method="plate"}.}
  \item{plate.size}{positive integer giving the number of consecutive spots corresponding to one plate or plate pack.  Ignored if \code{method="loess"}.}
  \item{plot}{logical. If \code{TRUE} then a scatter plot of the print order effect is sent to the current graphics device.}
}
\details{
Print-order is associated with the 384-well plates used in the printing of spotted microarrays.
There may be variations in DNA concentration or quality between the different plates.
The may be variations in ambient conditions during the time the array is printed.

This function is intended to pre-process the intensities before other normalization methods are applied to adjust for variations in DNA quality or concentration and other print-order effects.

Printorder means the order in which spots are printed on a microarray.
Spotted arrays are printed using a print head with an array of print-tips.
Spots in the various tip-groups are printed in parallel.
Printing is assumed to start in the top left hand corner of each tip-groups and to proceed right and down by rows, or else to start in the top right hand and to proceed left and down by rows.
See \code{\link{printorder}} for more details.
(WARNING: this is not always the case.)
This is true for microarrays printed at the Australian Genome Research Facility but might not be true for arrays from other sources.

If \code{object} is an \code{RGList} then printorder is performed for each intensity in each array.

\code{plotPrintorder} is a non-generic function which calls \code{normalizeForPrintorder} with \code{plot=TRUE}.
}
\value{
\code{normalizeForPrintorder} produces an \code{RGList} containing normalized intensities.

The function \code{plotPrintorder} or \code{normalizeForPrintorder.rg} with \code{plot=TRUE} returns no value but produces a plot as a side-effect.

\code{normalizeForPrintorder.rg} with \code{plot=FALSE} returns a list with the following components: 
  \item{R}{numeric vector containing the normalized red channel intensities}
  \item{G}{numeric vector containing the normalized red channel intensites}
  \item{R.trend}{numeric vector containing the fitted printorder trend for the red channel}
  \item{G.trend}{numeric vector containing the fitted printorder trend for the green channe}
}
\references{
Smyth, G. K. Print-order normalization of cDNA microarrays. March 2002.  \url{http://www.statsci.org/smyth/pubs/porder/porder.html}
}
\author{Gordon Smyth}
\seealso{
\code{\link{printorder}}.

An overview of LIMMA functions for normalization is given in \link{05.Normalization}.
}
\examples{
library(sma)
data(MouseArray)
plotPrintorder(mouse.data,mouse.setup,slide=1,separate=TRUE)
RG <- normalizeForPrintorder(mouse.data,mouse.setup)
}
\keyword{models}
