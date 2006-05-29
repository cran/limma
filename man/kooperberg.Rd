\title{Kooperberg Model-Based Background Correction}
\name{kooperberg}
\alias{kooperberg}
\description{
This function uses a Bayesian model to 
background correct GenePix microarray data.
}

\usage{
kooperberg(RG, a=TRUE, layout=RG$printer, verbose=TRUE)
}
\arguments{
\item{RG}{an RGList of GenePix data, read in using \code{read.maimages}, with \code{other.columns=c("F635 SD","B635 SD","F532 SD","B532 SD","B532 Mean","B635 Mean","F Pixels","B Pixels")}.}
\item{a}{logical.  If \code{TRUE}, the 'a' parameters in the model (equation 3 and 4) are estimated for each slide.  If \code{FALSE} the 'a' parameters are set to unity.}
\item{layout}{list containing print layout with components \code{ngrid.r}, \code{ngrid.c}, \code{nspot.r} and \code{nspot.c}.  Defaults to \code{RG$printer}.}
\item{verbose}{logical.  If \code{TRUE}, progress is reported to standard output.}
}

\details{
This function is for use with GenePix data and is designed to cope with the problem of large numbers of negative intensities and hence missing values on the log-intensity scale.
It avoids missing values in most cases and at the same time dampens down the variability of log-ratios for low intensity spots.
See Kooperberg et al (2002) for more details.

\code{kooperberg} uses the foreground and background intensities, standard
deviations and number of pixels to compute empirical estimates of the model 
parameters as described in equation 2 of Kooperberg et al (2002).
}

\value{
An \code{RGList} containing the components
\item{R}{matrix containing the background adjusted intensities for
the red channel for each spot for each array}
\item{G}{matrix containing the background adjusted intensities for the green channel for each spot for each array}
\item{printer}{list containing print layout}
}

\author{Matthew Ritchie}

\references{
Kooperberg, C., Fazzio, T. G., Delrow, J. J., and Tsukiyama, T. (2002)
Improved background correction for spotted DNA microarrays.
\emph{Journal of Computational Biology} \bold{9}, 55-66.
}
	
\seealso{
\link{04.Background} gives an overview of background correction functions defined in the LIMMA package.
}

\examples{
#  This is example code for reading and background correcting GenePix data
#  given GenePix Results (gpr) files in the working directory (data not
#  provided).
\dontrun{
genepixFiles <- dir(pattern="*\\\\.gpr$") # get the names of the GenePix image analysis output files in the current directory
RG <- read.maimages(genepixFiles, source="genepix", other.columns=c("F635 SD","B635 SD","F532 SD","B532 SD","B532 Mean","B635 Mean","F Pixels","B Pixels"))
RGmodel <- kooperberg(RG)
MA <- normalizeWithinArrays(RGmodel)
}
}

\keyword{models}
