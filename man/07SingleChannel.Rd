\name{07.SingleChannel}
\alias{07.SingleChannel}
\title{Individual Channel Analysis of Two-Color Microarrays}

\description{
This page gives an overview of the LIMMA functions fit linear models to two-color microarray data in terms of the log-intensities rather than log-ratios.

The function \code{\link{intraspotCorrelation}} estimates the intra-spot correlation between the two channels.
The regression function \code{\link{lmscFit}} takes the correlation as an argument and fits linear models to the two-color data in terms of the individual log-intensities.
The output of \code{lmscFit} is an \code{MArrayLM} object just the same as from \code{lmFit}, so inference proceeds in the same way as for log-ratios once the linear model is fitted.
See \link{06.LinearModels}.

The function \code{\link{targetsA2C}} converts two-color format target data frames to single channel format, i.e, converts from array-per-line to channel-per-line, to facilitate the formulation of the design matrix.
}

\author{Gordon Smyth}
\keyword{documentation}
