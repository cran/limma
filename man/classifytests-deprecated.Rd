\name{classifyTests-deprecated}
\alias{classifyTests}
\title{Treat Simultaneous T-Tests as Classification Problem}
\description{
Classify a series of related t-statistics as up, down or not significant.
}
\usage{
classifyTests(object, cor.matrix=NULL, design=NULL, contrasts=diag(ncol(design)), df=Inf, p.value=0.01, fstat.only=FALSE)
}
\arguments{
  \item{object}{numeric matrix of t-statistics or an \code{MArrayLM} object from which the t-statistics may be extracted.}
  \item{cor.matrix}{covariance matrix of each row of t-statistics.  Defaults to the identity matrix.}
  \item{design}{full rank numeric design matrix.  Not used if \code{cor.matrix} is specified.}
  \item{contrasts}{numeric matrix with columns specifying contrasts of the coefficients of interest.  Not used if \code{cor.matrix} is specified.}
  \item{df}{numeric vector giving the degrees of freedom for the t-statistics.
  May have length 1 or length equal to the number of rows of \code{tstat}.}
  \item{p.value}{numeric value between 0 and 1 giving the desired size of the test}
  \item{fstat.only}{logical, if \code{TRUE} then return the overall F-statistic as for \code{FStat} instead of classifying the test results}
}
\value{
An object of class \code{TestResults}.
This is essentially a numeric matrix with elements \code{-1}, \code{0} or \code{1} depending on whether each t-statistic is classified as significantly negative, not significant or significantly positive respectively.
}
\details{
\code{FStat} computes the F-statistic for testing all the contrasts equal to zero.
It is equivalent to \code{classifyTests} with \code{fstat.only=TRUE}.

\code{classifyTests} classifies using a nested F-test approach giving particular attention to correctly classifying genes which have two or more significant t-statistics, i.e., are differential expressed under two or more conditions.
The overall F-statistic used is computed by \code{FStat}.
At least one constrast will be classified as significant if and only if the overall F-statistic is significant.
\code{classifyTestsT} and \code{classifyTestsP} implement simpler classification schemes based on threshold or critical values for the individual t-statistics in the case of \code{classifyTestsT} or p-values obtained from the t-statistics in the case of \code{classifyTestsP}.

Rows of \code{tstat} correspond to genes and columns to coefficients or contrasts.
For each row of \code{tstat}, F-statistics are constructed from the t-statistics.
If the overall F-statistic is significant, then the function makes a best choice as to which t-statistics contributed to this result.
The methodology is based on the principle that any t-statistic should be called significant if the F-test is still significant for that row when all the larger t-statistics are set to the same absolute size as the t-statistic in question.

If \code{tstat} is an \code{MArrayLM} object, then all arguments except for \code{p.value} are extracted from it.

\code{cor.matrix} is the same as the correlation matrix of the coefficients from which the t-statistics are calculated.
If \code{cor.matrix} is not specified, then it is calculated from \code{design} and \code{contrasts} if at least \code{design} is specified or else defaults to the identity matrix.
In terms of \code{design} and \code{contrasts}, \code{cor.matrix} is obtained by standardizing the matrix
\code{ t(contrasts) \%*\% solve(t(design) \%*\% design) \%*\% contrasts }
to a correlation matrix.
}
\seealso{
An overview of linear model functions in limma is given by \link{5.LinearModels}.
}
\author{Gordon Smyth}
\keyword{htest}
