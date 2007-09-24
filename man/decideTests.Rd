\name{decideTests}
\alias{decideTests}
\title{Multiple Testing Across Genes and Contrasts}
\description{
Classify a series of related t-statistics as up, down or not significant.
A number of different multiple testing schemes are offered which adjust for multiple testing down the genes as well as across contrasts for each gene.
}
\usage{
decideTests(object,method="separate",adjust.method="BH",p.value=0.05,lfc=0)
}
\arguments{
  \item{object}{\code{MArrayLM} object output from \code{eBayes} from which the t-statistics may be extracted.}
  \item{method}{character string specify how probes and contrasts are to be combined in the multiple testing strategy.  Choices are \code{"separate"}, \code{"global"}, \code{"hierarchical"}, \code{"nestedF"} or any partial string.}
  \item{adjust.method}{character string specifying p-value adjustment method.  Possible values are \code{"none"}, \code{"BH"}, \code{"fdr"} (equivalent to \code{"BH"}), \code{"BY"} and \code{"holm"}. See \code{\link[stats]{p.adjust}} for details.}
  \item{p.value}{numeric value between 0 and 1 giving the desired size of the test}
  \item{lfc}{minimum log2-fold-change required}
}
\value{
An object of class \code{\link[limma:TestResults-class]{TestResults}}.
This is essentially a numeric matrix with elements \code{-1}, \code{0} or \code{1} depending on whether each t-statistic is classified as significantly negative, not significant or significantly positive respectively.

If \code{lfc>0} then contrasts are judged significant only when the log2-fold change is at least this large in absolute value.
For example, one might choose \code{lfc=log2(1.5)} to restrict to 50\% changes or \code{lfc=1} for 2-fold changes.
In this case, contrasts must satisfy both the p-value and the fold-change cutoff to be judged significant.
}
\details{
These functions implement multiple testing procedures for determining whether each statistic in a matrix of t-statistics should be considered significantly different from zero.
Rows of \code{tstat} correspond to genes and columns to coefficients or contrasts.

The setting \code{method="separate"} is equivalent to using \code{topTable} separately for each coefficient in the linear model fit, and will give the same lists of probes if \code{adjust.method} is the same.
\code{method="global"} will treat the entire matrix of t-statistics as a single vector of unrelated tests.
\code{method="hierarchical"} adjusts down genes and then across contrasts.
\code{method="nestedF"} adjusts down genes and then uses \code{classifyTestsF} to classify contrasts as significant or not for the selected genes.
}
\seealso{
An overview of multiple testing functions is given in \link{08.Tests}.
}
\author{Gordon Smyth}
\keyword{htest}
