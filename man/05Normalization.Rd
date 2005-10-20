\name{05.Normalization}
\alias{05.Normalization}
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
\code{normalizeBetweenArrays} uses utility functions \code{\link{normalizeMedianAbsValues}} and \code{\link{normalizeQuantiles}}, none of which need to be called directly by users.
}

\author{Gordon Smyth}
\references{
Smyth, G. K., and Speed, T. P. (2003). Normalization of cDNA microarray data. In: \emph{METHODS: Selecting Candidate Genes from DNA Array Screens: Application to Neuroscience}, D. Carter (ed.). Methods Volume 31, Issue 4, December 2003, pages 265-273. \url{http://www.statsci.org/smyth/pubs/normalize.pdf}
}
\keyword{documentation}
