\title{Kooperberg Model-Based Background Correction}
\name{kooperberg}
\alias{kooperberg}
\description{
This function uses a Bayesian model to background correct
data from a series of microarray experiments.
It currently works only with GenePix data.
}

\usage{
kooperberg(names, fg="mean", bg="median", a=TRUE, layout)
}
\arguments{
\item{names}{character vector giving the names of data.frames containing GenePix data}
\item{fg}{character string giving foreground estimator.
Choices are \code{"mean"} or \code{"median"}.}
\item{bg}{character string giving foreground estimator.
Choices are \code{"mean"} or \code{"median"}.}
\item{a}{logical.  If \code{TRUE}, the 'a' parameters in the model (equation 3 and 4) are estimated for each slide.  If \code{FALSE} the 'a' parameters are set to unity.}
\item{layout}{list containing print layout with components \code{ngrid.r}, \code{ngrid.c}, \code{nspot.r} and \code{nspot.c}}
}

\details{
This function is for use with GenePix data and is designed to cope with the problem of large numbers of negative intensities and hence missing values on the log-intensity scale.
It avoids missing values in most cases and at the same time dampens down the variability of log-ratios for low intensity spots.
See Kooperberg et al (2003) for more details.

\code{kooperberg} serially extracts the foreground and background intensities, standard
deviations and number of pixels from GenePix data frames.
This information is used to compute empirical estimates of the model parameters
as described in equation 2 of Kooperberg et al (2003).

The foreground and background estimates extracted from the GenePix files may be based on means or medians of pixel values.
Setting \code{fg="mean"} uses the GenePix column \code{F635.Mean} for red foreground and
the GenePix column \code{F532.Mean} for green foreground.
Setting \code{fg="median"} uses columns \code{F635.Median} and \code{F532.Median}.
Similarly for the background, \code{bg="mean"} uses columns \code{B635.Mean} and \code{B532.Mean} while \code{bg="median"} uses columns \code{B635.Median} and \code{B532.Median}.
}

\value{
An \code{RGList} containing the components
\item{R}{matrix containing the background adjusted intensities for
the red channel for each spot for each array}
\item{G}{matrix containing the background adjusted intensities for the green channel for each spot for each array}
\item{printer}{list containing printer layout, as input in \code{layout}}
}

\author{Matthew Ritchie}

\references{
Kooperberg, C., Fazzio, T. G., Delrow, J. J., and Tsukiyama, T. (2002)
Improved background correction for spotted DNA microarrays.
\emph{Journal of Computational Biology} \bold{9}, 55-66.
}
	
\seealso{
\link{4.Normalization} gives an overview of normalization and background correction functions defined in the LIMMA package.
}

\examples{
#  This is example code for reading and background correcting GenePix data
#  given GenePix Results (gpr) files in the working directory (data not
#  provided).
\dontrun{
genepixFiles <- dir(pattern="*\\\\.gpr$") # get the names of the GenePix image analysis output files in the current directory
read.series(genepixFiles, suffix=NULL, skip=26, sep="\t") # read the files into data.frames, assuming headers of 26 lines
layout <- getLayout(get(genepixFiles[1]))
RGmodel <- kooperberg(genepixFiles, layout=layout) # model-based background correction
MA <- normalizeWithinArrays(RGModel) # normalize the data
}
}

\keyword{models}
