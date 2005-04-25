\name{decideTests}
\alias{decideTests}
\title{Compute Matrix of Hypothesis Test Results}
\description{
Classify a series of related t-statistics as up, down or not significant.
}
\usage{
decideTests(object,method="separate",adjust.method="fdr",p.value=0.05)
}
\arguments{
  \item{object}{\code{MArrayLM} object output from \code{eBayes} from which the t-statistics may be extracted.}
  \item{method}{character string specify how probes and contrasts are to be combined in the multiple testing strategy.  Choices are \code{"separate"}, \code{"global"}, \code{"heirarchical"}, \code{"nestedF"} or any partial string.}
  \item{adjust.method}{character string specifying p-value adjustment method.  See \code{\link[stats]{p.adjust}} for possible values.}
  \item{p.value}{numeric value between 0 and 1 giving the desired size of the test}
}
\value{
An object of class \code{\link[limma:TestResults-class]{TestResults}}.
This is essentially a numeric matrix with elements \code{-1}, \code{0} or \code{1} depending on whether each t-statistic is classified as significantly negative, not significant or significantly positive respectively.
}
\details{
These functions implement multiple testing procedures for determining whether each statistic in a matrix of t-statistics should be considered significantly different from zero.
Rows of \code{tstat} correspond to genes and columns to coefficients or contrasts.

The setting \code{method="separate"} is equivalent to using \code{topTable} separately for each coefficient in the linear model fit, and will give the same lists of probes if \code{adjust.method} is the same.
Note that the defaults for \code{adjust.method} are different for \code{decideTests} and \code{topTable}.
\code{method="global"} will treat the entire matrix of t-statistics as a single vector of unrelated tests.
\code{method="heirarchical"} adjusts down genes and then across contrasts.
\code{method="nestedF"} adjusts down genes and then uses \code{classifyTestsF} to classify contrasts as significant or not for the selected genes.
}
\seealso{
An overview of linear model functions in limma is given by \link{5.LinearModels}.
}
\author{Gordon Smyth}
\keyword{htest}
