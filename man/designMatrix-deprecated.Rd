\name{designMatrix-deprecated}
\alias{designMatrix}
\title{Construct Design Matrix}
\description{
Construct design matrix from RNA target information.
}
\usage{
designMatrix(targets, parameters, ref, verbose=TRUE)
}
\arguments{
  \item{targets}{matrix or data.frame with columns \code{Cy3} and \code{Cy5} specifying which RNA was hybridized to each array}
  \item{parameters}{matrix specifying contrasts between RNA samples which should correspond to regression coefficients.
  Row names should correspond to unique RNA sample names found in \code{targets}.}
  \item{ref}{character string giving name of common reference RNA if such exists.
  Exactly one argument of \code{parameters} or \code{ref} should be specified.}
  \item{verbose}{logical, if \code{TRUE} then unique names found in \code{targets} will be printed to standard output}
}

\details{
This function is intended to produce a design matrix for use in functions \code{lmFit} etc for two-color microarray experiments.
}

\value{
\code{designMatrix} produces a design matrix with row names as in \code{targets} and column names as in \code{parameters}.

\code{uniqueTargets} produces a character vector of unique target names from the columns \code{Cy3} and \code{Cy5} of \code{targets}.
}

\seealso{
\code{\link[stats]{model.matrix}} in the base package.

An overview of linear model functions in limma is given by \link{5.LinearModels}.
}

\author{Gordon Smyth}

\keyword{regression}
