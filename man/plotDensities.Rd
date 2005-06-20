\title{Individual-channel Densities Plot}
\name{plotDensities}
\alias{plotDensities}
\description{
Plots the densities of individual-channel intensities for two-color microarray data.
}
\usage{
plotDensities(object, log=TRUE, arrays=NULL, singlechannels=NULL, groups=NULL, col=NULL) 
}
\arguments{
  \item{object}{an \code{RGList} or \code{MAList} object. \code{RGList} objects containing logged or unlogged intensities can be accommodated using the \code{log.transform} argument.}

  \item{log}{logical, should densities be formed and plotted for the log-intensities (\code{TRUE}) or raw intensities (\code{FALSE})?}
  
  \item{arrays}{vector of integers giving the arrays from which the
    individual-channels will be selected to be plotted.
    Corresponds to columns of \code{M}
    and \code{A} (or \code{R} and \code{G}).  Defaults to all arrays.}
  
  \item{singlechannels}{vector of integers indicating which
    individual-channels will be selected to be plotted.  Values correspond
    to the columns of the matrix of \code{cbind(R,G)} and range
    between \code{1:ncol(R)} for red channels and
    \code{( (ncol(R)+1):(ncol(R)+ncol(G)) )} for the green
    channels in \code{object}.
    Defaults to all channels.}
  
  \item{groups}{vector of consecutive integers beginning at 1 indicating
    the groups of arrays or individual-channels (depending on which of
    \code{arrays} or \code{singlechannels} are non \code{NULL}).  This is used
    to color any groups of the individual-channel densities.
    If \code{NULL} (default), \code{groups} correspond to the
    red and green channels.  If both \code{arrays} and
    \code{singlechannels} are \code{NULL} all arrays are selected and
    groups (if specified) must correspond to the arrays.}

  \item{col}{vector of colors of the same length as the number of
    different groups. If \code{NULL} (default) the \code{col} equals
    \code{c("red","green")}.  See details for more specifications.}
}

\details{
This function is used as a data display technique associated with between-array normalization, especially individual-channel normalization methods such as quantile-normalization.
See the section on between-array normalization in the LIMMA User's Guide.

If no \code{col} is specified, the default is to color individual channels
according to red and green. If both \code{arrays} and \code{groups} are
non-\code{NULL}, then the length of \code{groups} must equal the length
of \code{arrays} and the maximum of \code{groups} (i.e. the number of
groups) must equal the length of \code{col} otherwise the default color
of black will be used for all individual-channels.
If \code{arrays} is \code{NULL} and both \code{singlechannels} and
\code{groups} are non-\code{NULL}, then the length of \code{groups} must
equal the length of \code{singlechannels} and the maximum of \code{groups}
(i.e. the number of groups) must equal the length of \code{col}
otherwise the default color of black will be used for all individual-channels.
}
\value{A plot is created on the current graphics device.}
\author{Natalie Thorne}
\seealso{
An overview of diagnostic plots in LIMMA is given in \link{09.Diagnostics}.
There is a section using \code{plotDensities} in conjunction with between-array normalization
in the \emph{\link[limma:../doc/usersguide]{LIMMA User's Guide}}.
}
\examples{
library(sma)
data(MouseArray)

#  no normalization but background correction is done
MA.n <- MA.RG(mouse.data)

#  Default settings for plotDensities.
plotDensities(MA.n)

#  One can reproduce the default settings.
plotDensities(MA.n,arrays=c(1:6),groups=c(rep(1,6),rep(2,6)),
col=c("red","green"))

#  Color R and G individual-channels by blue and purple.
plotDensities(MA.n,arrays=NULL,groups=NULL,col=c("blue","purple"))

#  Indexing individual-channels using singlechannels (arrays=NULL).
plotDensities(MA.n,singlechannels=c(1,2,7))

#  Change the default colors from c("red","green") to c("pink","purple")
plotDensities(MA.n,singlechannels=c(1,2,7),col=c("pink","purple"))

#  Specified too many colors since groups=NULL defaults to two groups.
plotDensities(MA.n,singlechannels=c(1,2,7),col=c("pink","purple","blue"))

#  Three individual-channels, three groups, three colors.
plotDensities(MA.n,singlechannels=c(1,2,7),groups=c(1,2,3),
col=c("pink","purple","blue"))

#  Three individual-channels, one group, one color.
plotDensities(MA.n,singlechannels=c(1,2,7),groups=c(1,1,1),
col=c("purple"))

#  All individual-channels, three groups (ctl,tmt,reference), three colors.
plotDensities(MA.n,singlechannels=c(1:12),
groups=c(rep(1,3),rep(2,3),rep(3,6)),col=c("darkred","red","green"))

}
\keyword{hplot}
