\name{04.Background}
\alias{04.Background}
\title{Background Correction}

\description{
This page deals with background correction methods for two-color microarray data.

Usually one doesn't need to explicitly ask for background correction of the intensities because this is done by default by \code{\link{normalizeWithinArrays}},
which subtracts the background from the foreground intensities before applying the normalization method.
This default background correction method can be over-ridden by using \code{\link{backgroundCorrect}} which offers a number of alternative
background correct methods to simple subtraction.
The function \code{backgroundCorrect} is used to correct the \code{RGList} before applying \code{normalizeWithinArrays}.

The \code{movingmin} method of \code{backgroundCorrect} uses utility functions \code{\link{ma3x3.matrix}} and \code{\link{ma3x3.spottedarray}}.

The \code{normexp} method of \code{backgroundCorrect} uses ulility functions \code{\link{normexp.fit}}, \code{\link{normexp.signal}}, \code{\link{normexp.m2loglik}}, \code{\link{normexp.grad}}, \code{\link{dnormexp}} and \code{\link{dnormexp.saddle}}.

\code{\link{kooperberg}} is a Bayesian background correction tool designed specifically for GenePix data.
\code{kooperberg} is not currently used as the default method for GenePix data because it is computationally intensive.
It requires several additional columns from the GenePix data files which can be read in using \code{read.maimages} and specifying the \code{other.columns} argument.
}

\author{Gordon Smyth}
\keyword{documentation}
