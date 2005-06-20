\name{modelMatrix}
\alias{modelMatrix}
\alias{uniqueTargets}
\title{Construct Design Matrix}
\description{
Construct design matrix from RNA target information for a two colour microarray experiment.
}
\usage{
modelMatrix(targets, parameters, ref, verbose=TRUE)
uniqueTargets(targets)
}
\arguments{
  \item{targets}{matrix or data.frame with columns \code{Cy3} and \code{Cy5} specifying which RNA was hybridized to each array}
  \item{parameters}{matrix specifying contrasts between RNA samples which should correspond to regression coefficients.
  Row names should correspond to unique RNA sample names found in \code{targets}.}
  \item{ref}{character string giving name of one of the RNA sources to be treated as reference.
  Exactly one argument of \code{parameters} or \code{ref} should be specified.}
  \item{verbose}{logical, if \code{TRUE} then unique names found in \code{targets} will be printed to standard output}
}

\details{
This function computes a design matrix for input to \code{lmFit} when analysing two-color microarray experiments in terms of log-ratios.

If the argument \code{ref} is used, then the experiment is treated as a one-way layout and the coefficients measure expression changes relative to the RNA source specified by \code{ref}.
The RNA source \code{ref} is often a common reference which appears on every array or is a control sample to which all the others are compared.
There is no restriction however.
One can choose \code{ref} to be any of the RNA sources appearing the \code{Cy3} or \code{Cy5} columns of \code{targets}.

If the \code{parameters} argument is set, then the columns of this matrix specify the comparisons between the RNA sources which are of interest.
This matrix must be of size n by (n-1), where n is the number of unique RNA sources found in \code{Cy3} and \code{Cy5}, and must have row names which correspond to the RNA sources.
}

\value{
\code{modelMatrix} produces a numeric design matrix with row names as in \code{targets} and column names as in \code{parameters}.

\code{uniqueTargets} produces a character vector of unique target names from the columns \code{Cy3} and \code{Cy5} of \code{targets}.
}

\seealso{
\code{\link[stats]{model.matrix}} in the stats package.

An overview of linear model functions in limma is given by \link{06.LinearModels}.
}

\author{Gordon Smyth}

\examples{
targets <- cbind(Cy3=c("Ref","Control","Ref","Treatment"),Cy5=c("Control","Ref","Treatment","Ref"))
rownames(targets) <- paste("Array",1:4)

parameters <- cbind(C=c(-1,1,0),T=c(-1,0,1))
rownames(parameters) <- c("Ref","Control","Treatment")

modelMatrix(targets, parameters)
modelMatrix(targets, ref="Ref")
}

\keyword{regression}
