\name{4.Normalization}
\alias{4.Normalization}
\title{Normalization of Microarray Data}

\description{
This page gives an overview of the LIMMA functions available to normalize data from spotted two-colour microarrays.
Smyth and Speed (2003) give an overview of the normalization techniques implemented in the functions.

Usually data from spotted microarrays will be normalized using \code{\link{normalizeWithinArrays}}.
A minority of data will also be normalized using \code{\link{normalizeBetweenArrays}} if diagnostic plots suggest a difference in scale between the arrays.

In rare circumstances, data might be normalized using \code{\link{normalizeForPrintorder}} before using \code{\link{normalizeWithinArrays}}.

All the normalization routines take account of spot quality weights which might be set in the data objects.
The weights can be temporarily modified using \code{\link{modifyWeights}} to, for example, remove ratio control spots from the normalization process.

If one is planning analysis of single-channel information from the microarrays rather than analysis of differential expression based on log-ratios, then the data should be normalized using a single channel-normalization technique.
Single channel normalization uses further options of the \code{\link{normalizeBetweenArrays}} function.
For more details see the \emph{\link[limma:../doc/usersguide]{LIMMA User's Guide}} which includes a section on single-channel normalization.

\code{normalizeWithinArrays} uses utility functions \code{\link{MA.RG}}, \code{\link{loessFit}} and \code{\link{normalizeRobustSpline}}.
\code{normalizeBetweenArrays} uses utility functions \code{\link{normalizeMedians}}, \code{\link{normalizeMedianDeviations}} and \code{\link{normalizeQuantiles}}, none of which need to be called directly by users.
}

\section{Backgound Correction}{
Usually one doesn't need to explicitly ask for background correction of the intensities because this is done by default by \code{\link{normalizeWithinArrays}},
which subtracts the background from the foreground intensities before applying the normalization method.
This default background correction method can be over-ridden by using \code{\link{backgroundCorrect}} which offers a number of alternative
background correct methods to simple subtraction.
Simply use \code{backgroundCorrect} to correct the \code{RGList} before applying \code{normalizeWithinArrays}.

\code{backgroundCorrect} uses utility functions \code{\link{ma3x3.matrix}}, \code{\link{ma3x3.spottedarray}}, \code{\link{fit.normexp}}, \code{\link{signal.normexp}} and \code{\link{m2loglik.normexp}}.

\code{\link{kooperberg}} is a Bayesian background correction tool designed specifically for GenePix data.
\code{kooperberg} is not currently used as the default method for GenePix data because it is computationally intensive.
It requires several columns of the GenePix data files which are not read in by read.maimages, so you will need to use \code{read.series} instead of \code{read.maimages} if you wish to use \code{kooperberg}.
}

\author{Gordon Smyth}
\references{
Smyth, G. K., and Speed, T. P. (2003). Normalization of cDNA microarray data. In: \emph{METHODS: Selecting Candidate Genes from DNA Array Screens: Application to Neuroscience}, D. Carter (ed.). Methods Volume 31, Issue 4, December 2003, pages 265-273. \url{http://www.statsci.org/smyth/pubs/normalize.pdf}
}
\keyword{documentation}
