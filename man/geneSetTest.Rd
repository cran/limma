\name{geneSetTest}
\alias{geneSetTest}
\title{Gene Set Test}
\description{
Test whether a given statistic is larger or over-represented in a given subset of genes.
}
\usage{
geneSetTest(selected,statistics,nsim=10000)
}
\arguments{
  \item{selected}{vector specifying the elements of \code{statistic} in the test group.  This can be a vector of indices, or a logical vector of the same length as \code{statistics}, or any vector such as \code{statistic[selected]} contains the statistic values for the selected group.}
  \item{statistics}{numeric vector giving the values of the test statistic for every gene or probe in the reference set, usually every probe on the microarray.}
  \item{nsim}{number of random samples to take in computing the p-value.}
}
\value{
Numeric value giving the estimated p-value.
}
\details{
This function computes the proportion of randomly selected groups from the set of all statistics which have mean statistic greater or equal to that of the test group.
}
\author{Gordon Smyth}
\keyword{htest}
