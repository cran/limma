\name{geneSetTest}
\alias{geneSetTest}
\title{Gene Set Test}
\description{
Test whether a given statistic is larger or over-represented in a given subset of genes.
}
\usage{
geneSetTest(selected,statistics,alternative="two.sided",nsim=10000,ranks.only=FALSE)
}
\arguments{
  \item{selected}{vector specifying the elements of \code{statistic} in the test group.  This can be a vector of indices, or a logical vector of the same length as \code{statistics}, or any vector such as \code{statistic[selected]} contains the statistic values for the selected group.}
  \item{statistics}{numeric vector giving the values of the test statistic for every gene or probe in the reference set, usually every probe on the microarray.}
  \item{alternative}{character string specifying the alternative hypothesis,
          must be one of \code{"two.sided"} (default), \code{"greater"} or
          \code{"less"}.  You can specify just the initial letter.}
  \item{nsim}{number of random samples to take in computing the p-value.
  Not used if \code{ranks.only=TRUE}.}
  \item{ranks.only}{logical, should the values \code{statistics} be used only to rank the genes or does it make sense to average statistics for selected sets?} 
}
\value{
Numeric value giving the estimated p-value.
}
\details{
This function computes a p-value to test the hypothesis that the selected genes tend to be more highly ranked on the given statistic.
If it makes sense to average values of the statistic, which would be so for example if the statistic was a t-statistic, then a permutation test is conducted.
In that case the function returns the proportion of \code{nsim} randomly selected groups from the set of all statistics which have mean statistic equal or more extreme than that of the test group.

If it doesn't make sense to average the values of the statistic for any reason, then only the ranks of the statistics are used and a Wilcoxon two-sample test, also known as a Mann-Whitney test, is performed.

This is essentially a stream-lined approach to Gene Set Enrichment Analysis introduced by Mootha et al (2003).
}
\seealso{\code{\link[stats]{wilcox.test}}}
\author{Gordon Smyth}

\references{
Mootha, V. K.,
Lindgren, C. M.,
Eriksson, K. F.,
Subramanian, A.,
Sihag, S.,
Lehar, J.,
Puigserver, P.,
Carlsson, E.,
Ridderstrale, M.,
Laurila, E.,
Houstis, N.,
Daly, M. J.,
Patterson, N.,
Mesirov, J. P.,
Golub, T. R.,
Tamayo, P.,
Spiegelman, B.,
Lander, E. S.,
Hirschhorn, J. N.,
Altshuler, D.,
Groop, L. C.
(2003). 
PGC-1alpha-responsive genes involved in oxidative phosphorylation are coordinately downregulated in human diabetes.
\emph{Nature Genetics} 34, 267-273.
}

\examples{
geneSetTest(c(0.1,3.3,2.2),-9:9,nsim=100)
geneSetTest(c(0.1,3.3,2.2),-9:9,ranks.only=TRUE)
}
\keyword{htest}
