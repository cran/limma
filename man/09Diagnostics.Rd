\name{09.Diagnostics}
\alias{09.Diagnostics}
\title{Diagnostics and Quality Assessment}

\description{
This page gives an overview of the LIMMA functions available for microarray quality assessment and diagnostic plots.

This package provides an \code{\link[limma:anova-method]{anova}} method which is designed for assessing the quality of an array series or of a normalization method.
It is not designed to assess differential expression of individual genes.
\code{\link[limma:anova-method]{anova}} uses utility functions \code{\link{bwss}} and \code{\link{bwss.matrix}}.

The function \code{\link{arrayWeights}} estimates the empirical reliability of each array following a linear model fit.

Diagnostic plots can be produced by

\describe{

\item{ \code{\link{imageplot}} }{
Produces a spatial picture of any spot-specific measure from an array image. If the log-ratios are plotted, then this produces an in-silico representation of the well known false-color TIFF image of an array.
\code{\link{imageplot3by2}} will write imageplots to files, six plots to a page.}

\item{ \code{\link{plotFB}} }{
Plots foreground versus background log-intensies for an array.}

\item{ \code{\link{plotMA}} }{
MA-plots.
One of the most useful plots of a two-color array.
\code{\link{plotMA3by2}} will write MA-plots to files, six plots to a page.
\code{\link{mdplot}} can also be useful for comparing two one-channel microarrays.}

\item{ \code{\link{plotPrintTipLoess}} }{
Produces a grid of MA-plots, one for each print-tip group on an array, together with the corresponding lowess curve.
Intended to help visualize print-tip loess normalization.}

\item{ \code{\link{plotPrintorder}} }{
For an array, produces a scatter plot of log-ratios or log-intensities by print order.}

\item{ \code{\link{plotDensities}} }{
Individual channel densities for one or more arrays.
An essential plot to accompany between array normalization, especially quantile normalization.}
}

\code{plotPrintTipLoess} uses utility functions \code{\link{gridr}} and \code{\link{gridc}}.
\code{\link{plotDensities}} uses utility function \code{\link{RG.MA}}.
}

\author{Gordon Smyth}
\keyword{documentation}
