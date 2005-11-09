\name{classifyTests}
\alias{FStat}
\alias{classifyTestsF}
\alias{classifyTestsT}
\alias{classifyTestsP}
\title{Multiple Testing Genewise Across Contrasts}
\description{
For each gene, classify a series of related t-statistics as up, down or not significant.
}
\usage{
classifyTestsF(object, cor.matrix=NULL, df=Inf, p.value=0.01, fstat.only=FALSE)
classifyTestsT(object, t1=4, t2=3)
classifyTestsP(object, df=Inf, p.value=0.05, method="holm")
FStat(object, cor.matrix=NULL)
}
\arguments{
  \item{object}{numeric matrix of t-statistics or an \code{MArrayLM} object from which the t-statistics may be extracted.}
  \item{cor.matrix}{covariance matrix of each row of t-statistics.  Defaults to the identity matrix.}
  \item{df}{numeric vector giving the degrees of freedom for the t-statistics.
  May have length 1 or length equal to the number of rows of \code{tstat}.}
  \item{p.value}{numeric value between 0 and 1 giving the desired size of the test}
  \item{fstat.only}{logical, if \code{TRUE} then return the overall F-statistic as for \code{FStat} instead of classifying the test results}
  \item{t1}{first critical value for absolute t-statistics}
  \item{t2}{second critical value for absolute t-statistics}
  \item{method}{character string specifying p-value adjustment method.  See \code{\link[stats]{p.adjust}} for possible values.}
}
\value{
An object of class \code{TestResults}.
This is essentially a numeric matrix with elements \code{-1}, \code{0} or \code{1} depending on whether each t-statistic is classified as significantly negative, not significant or significantly positive respectively.

\code{FStat} produces a numeric vector of F-statistics with attributes \code{df1} and \code{df2} giving the corresponding degrees of freedom.
}
\details{
Note that these functions do not adjust for multiple testing across genes.
The adjustment for multiple testing is across the contrasts rather than the more usual control across genes.
The functions described here are called by \code{\link{decideTests}}.
Most users should use \code{decideTests} rather than using these functions directly.

These functions implement multiple testing procedures for determining whether each statistic in a matrix of t-statistics should be considered significantly different from zero.
Rows of \code{tstat} correspond to genes and columns to coefficients or contrasts.

\code{FStat} computes the gene-wise F-statistics for testing all the contrasts equal to zero.
It is equivalent to \code{classifyTestsF} with \code{fstat.only=TRUE}.

\code{classifyTestsF} uses a nested F-test approach giving particular attention to correctly classifying genes which have two or more significant t-statistics, i.e., are differential expressed under two or more conditions.
For each row of \code{tstat}, the overall F-statistics is constructed from the t-statistics as for \code{FStat}.
At least one constrast will be classified as significant if and only if the overall F-statistic is significant.
If the overall F-statistic is significant, then the function makes a best choice as to which t-statistics contributed to this result.
The methodology is based on the principle that any t-statistic should be called significant if the F-test is still significant for that row when all the larger t-statistics are set to the same absolute size as the t-statistic in question.

\code{classifyTestsT} and \code{classifyTestsP} implement simpler classification schemes based on threshold or critical values for the individual t-statistics in the case of \code{classifyTestsT} or p-values obtained from the t-statistics in the case of \code{classifyTestsP}.
For \code{classifyTestsT}, classifies any t-statistic with absolute greater than \code{t2} as significant provided that at least one t-statistic for that gene is at least \code{t1} in absolute value.
\code{classifyTestsP} applied p-value adjustment from \code{p.adjust} to the p-values for each gene.

If \code{tstat} is an \code{MArrayLM} object, then all arguments except for \code{p.value} are extracted from it.

\code{cor.matrix} is the same as the correlation matrix of the coefficients from which the t-statistics are calculated.
If \code{cor.matrix} is not specified, then it is calculated from \code{design} and \code{contrasts} if at least \code{design} is specified or else defaults to the identity matrix.
In terms of \code{design} and \code{contrasts}, \code{cor.matrix} is obtained by standardizing the matrix
\code{ t(contrasts) \%*\% solve(t(design) \%*\% design) \%*\% contrasts }
to a correlation matrix.
}
\seealso{
An overview of multiple testing functions is given in \link{08.Tests}.
}
\author{Gordon Smyth}
\examples{
tstat <- matrix(c(0,5,0, 0,2.5,0, -2,-2,2, 1,1,1), 4, 3, byrow=TRUE)
classifyTestsF(tstat)

# See also the examples for contrasts.fit and vennDiagram
}
\keyword{htest}
