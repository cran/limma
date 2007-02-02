\name{mergeScans}
\alias{mergeScans}
\alias{mergeScansRG}
\title{Merge two scans of two-color arrays}
\description{
Merge two sets of intensities of two-color arrays that are scanned twice at two different scanner settings, one at a lower gain setting with no saturated spot intensities and the other at a higher gain setting with a higher signal-to-noise ratio and some saturated spot intensities. 
}

\usage{
mergeScansRG(RGlow, RGhigh, AboveNoiseLowG=NULL, AboveNoiseLowR=NULL, outlierp=0.01)
}
\arguments{
  \item{RGlow}{object of class \code{RGList} containing red and green intensities constituting two-color microarray data scanned at a lower gain setting.}
  \item{RGhigh}{object of class \code{RGList} containing red and green intensities constituting two-color microarray data scanned at a higher gain setting.}
  \item{AboveNoiseLowG}{matrix of 1 or 0 for low scan intensities of green color, 1 for spots above noise level or 0 otherwise. One column per array.}
  \item{AboveNoiseLowR}{matrix of 1 or 0 for low scan intensities of red color, 1 for spots above noise level or 0 otherwise. One column per array. }
  \item{outlierp}{p-value for outliers. 0 for no outlier detection or any value between 0 and 1. Default p-value is 0.01.}
} 

\details{
This function merges two separate scans of each fluorescent label on a two-color array scanned at two different scanner settings by using a nonlinear regression model consisting of two linear regression lines and a quadratic function connecting the two, which looks like a hockey stick.
The changing point, i.e. the saturation point, in high scan is also estimated as part of model.
Signals produced for certain spots can sometimes be very low (below noise) or too high (saturated) to be accurately read by the scanner.
The proportions of spots that are below noise or above saturation are affected by the settings of the laser scanner used to read the arrays, with low scans minimizing saturation effects and high scans maximizing signal-to-noise ratios. Saturated spots can cause bias in intensity ratios that cannot be corrected for using conventional normalization methods.

Each fluorescent label on a two-color array can be scanned twice: for example, a high scan targeted at reaching saturation level for the brightest 1 percent of the spots on the array, and a low scan targeted at the lowest level of intensity which still allowed accurate grid placement on the arrays.
By merging data from two separate laser scans of each fluorescent label on an array, we can avoid the potential bias in signal intensities due to below noise or above saturation and, thus provide better estimates of true differential expression as well as increase usable spots.

The merging process is designed to retain signal intensities from the high scan except when scanner saturation causes the high scan signal to be under-measured. The saturated spots are predicted from the corresponding low scans by the fitted regression model.
It also checks any inconsistency between low and high scans.
}

\value{
An object of class \code{\link{RGList-class}} with the following components:

\item{G}{numeric matrix containing the merged green (cy3) foreground intensities. Rows correspond to spots and columns to arrays.}

\item{R}{numeric matrix containing the merged red (cy5) foreground intensities. Rows correspond to spots and columns to arrays.}

\item{Gb}{numeric matrix containing the green (cy3) background intensities from high scan.}

\item{Rb}{numeric matrix containing the red (cy5) background intensities from high scan.}

\item{other}{list numeric matrices \code{Gsaturated}, \code{Rsatured}, \code{Goutlier} and \code{Routlier}.
The first two contain saturation flags (1=saturated, 0=otherwise) for the green (cy3) and red (Cy5) channels of the high scan.
The second two contain outlier flags (1=outlier, 0=otherwise) for the green (cy3) and red (Cy5) channels.}
}

\author{Dongseok Choi \email{choid@ohsu.edu}.}

\references{
Choi D, O'Malley JP, Lasarev MR, Lapidus J, Lu X, Pattee P,  Nagalla SR (2006). Extending the Dynamic Range of Signal Intensities in DNA Microarrays. \emph{Online Journal of Bioinformatics}, \bold{7}, 46-56.
}

\examples{
\dontrun{
#RG1: An RGList from low scan
#RG2: An RGList from high scan
RGmerged <- mergeScansRG(RG1,RG2,AboveNoiseLowG=ANc3,AboveNoiseLowR=ANc5)

#merge two scans when all spots are above noise in low scan and no outlier detection.
RGmerged <- mergeScansRG(RG1,RG2,outlierp=0)
}
}

\keyword{models}
