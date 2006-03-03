\name{anova.MAList-method}
\docType{methods}
\alias{anova.MAList}
\title{ANOVA Table - method}
\description{
Analysis of variance method for objects of class \code{MAList}.
Produces an ANOVA table useful for quality assessment by decomposing between and within gene sums of squares for a series of replicate arrays.
This method produces a single ANOVA Table rather than one for each gene and is not used to identify differentially expressed genes.
}
\section{Usage}{
\code{anova(object,design=NULL,ndups=2,...)}
}
\section{Arguments}{
\describe{
  \item{\code{object}}{object of class \code{MAList}. Missing values in the M-values are not allowed.}
  \item{\code{design}}{numeric vector or single-column matrix containing the design matrix for linear model. The length of the vector or the number of rows of the matrix should agree with the number of columns of M.}
  \item{\code{ndups}}{number of duplicate spots. Each gene is printed ndups times in adjacent spots on each array.}
  \item{\code{...}}{other arguments are not used}
}
}
\section{Details}{
This function aids in quality assessment of microarray data and in the comparison of normalization methodologies.
It applies only to replicated two-color experiments in which all the arrays are hybridized with the same RNA targets, possibly with dye-swaps, so the design matrix should have only one column.
The function has not been heavily used and is somewhat experimental.
}
\section{Value}{
  An object of class \code{anova} containing rows for between genes, between arrays, gene x array interaction, and between duplicate with array sums of squares.
  Variance components are estimated for each source of variation.
}
\note{
This function does not give valid results in the presence of missing M-values.
}
\seealso{
\code{\link{MAList-class}}, \code{\link{bwss.matrix}}, \code{\link[base:anova]{anova}}.

An overview of quality assessment and diagnostic functions in LIMMA is given by \link{09.Diagnostics}.
}
\author{Gordon Smyth}
\keyword{models}
