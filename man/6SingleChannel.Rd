\name{6.SingleChannel}
\alias{6.SingleChannel}
\title{Individual Channel Analysis of Two-Color Microarrays}

\description{
This page gives an overview of the LIMMA functions fit linear models to two-color microarray data in terms of the log-intensities rather than log-ratios.

The function \code{\link{intraspotCorrelation}} estimates the intra-spot correlation between the two channels.
The regression function \code{\link{lmscFit}} takes the correlation as an argument and fits linear models to the two-color data in terms of the individual log-intensities.
The output of \code{lmscFit} is an \code{MArrayLM} object just the same as from \code{lmFit}, so inference proceeds in the same way as for log-ratios once the linear model is fitted.
See \link{5.LinearModels}.

The function \code{\link{array2channel}} converts two-color format target data frames to single channel format, to facilitate the formulation of the design matrix.
}

\author{Gordon Smyth}
\keyword{documentation}
