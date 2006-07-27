\name{geneSetTest}
\alias{geneSetTest}
\title{Gene Set Test}
\description{
Test whether a set of genes is enriched for differential expression.
}
\usage{
geneSetTest(selected,statistics,alternative="mixed",type="auto",ranks.only=TRUE,nsim=10000)
}
\arguments{
  \item{selected}{vector specifying the elements of \code{statistic} in the test group.  This can be a vector of indices, or a logical vector of the same length as \code{statistics}, or any vector such as \code{statistic[selected]} contains the statistic values for the selected group.}
  \item{statistics}{numeric vector giving the values of the test statistic for every gene or probe in the reference set, usually every probe on the microarray.}
  \item{alternative}{character string specifying the alternative hypothesis, must be one of \code{"mixed"} (default), \code{"either"}, \code{"up"} or \code{"down"}. \code{two.sided}, \code{"greater"} and \code{"less"} are also permitted as synonyms for \code{"either"}, \code{"up"} and \code{"down"} respectively.}
  \item{type}{character string specifying whether the statistics are t-like (\code{"t"}), F-like \code{"f"} or whether the function should make an educated guess (\code{"auto"})}
  \item{ranks.only}{logical, if \code{TRUE} only the ranks of the \code{statistics} are used.} 
  \item{nsim}{number of random samples to take in computing the p-value.
  Not used if \code{ranks.only=TRUE}.}
}
\value{
Numeric value giving the estimated p-value.
}
\details{
This function computes a p-value to test the hypothesis that the selected genes tend to be differentially expressed.
The aim here is the same as for Gene Set Enrichment Analysis introduced by Mootha et al (2003), but the statistical tests used are different.

The \code{statistics} are usually a set of probe-wise statistics arising for some comparison from a microarray experiment.
They may be t-statistics, meaning that the genewise null hypotheses would be rejected for large positive or negative values, or they may be F-statistics, meaning that only large values are significant.
Any set of signed statistics, such as log-ratios, M-values or moderated t-statistics, are treated as t-like.
Any set of unsigned statistics, such as F-statistics, posterior probabilities or chi-square tests are treated as F-like.
If \code{type="auto"} then the statistics will be taken to be t-like if they take both positive and negative values and otherwise will be taken to be F-like.

There are four possible alternatives to test for.
\code{alternative=="up"} means the genes in the set tend to be up-regulated, with positive t-statistics.
\code{alternative=="down"} means the genes in the set tend to be down-regulated, with negative t-statistics.
\code{alternative=="either"} means the set is either up or down-regulated as a whole.
\code{alternative=="mixed"} test whether the genes in the set tend to be differentially expressed, without regard for direction.
In this case, the test will be significant if the set contains mostly large test statistics, even if some are positive and some are negative.

The first three alternativea appropriate if you have a prior expection that all the genes in the set will react in the same direction.
The \code{"mixed"} alternative is appropriate if you know only that the genes are involved in the relevant pathways, without knowing the direction of effect for each gene.
The \code{"mixed"} alternative is the only one possible with F-like statistics.

The test statistic used for the gene-set-test is the mean of the statistics in the set.
If \code{ranks.only} is \code{TRUE} the only the ranks of the statistics are used.
In this case the p-value is obtained from a Wilcoxon test.
If \code{ranks.only} is \code{FALSE}, then the p-value is obtained by simulation using \code{nsim} random selected sets of genes.
}

\note{
This function can be used to detect differential expression for a group of genes, even when the effects are too small or there is too little data to detect the genes individually.
It is also provides a means to compare the results between different experiments.
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
stat <- -9:9
sel <- c(2,4,5)
geneSetTest(sel,stat,alternative="down")
geneSetTest(sel,stat,alternative="either")
geneSetTest(sel,stat,alternative="down",ranks=FALSE)
sel <- c(1,19)
geneSetTest(sel,stat,alternative="mixed")
geneSetTest(sel,stat,alternative="mixed",ranks=FALSE)
}
\keyword{htest}
