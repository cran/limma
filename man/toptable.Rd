\name{toptable}
\alias{toptable}
\alias{topTable}
\title{Table of Top Genes from Linear Model Fit}
\description{
Extract a table of the top-ranked genes from a linear model fit.
}
\usage{
toptable(fit,coef=1,number=10,genelist=NULL,A=NULL,eb=NULL,adjust.method="holm",sort.by="B",resort.by=NULL,...)
topTable(fit,coef=1,number=10,genelist=NULL,adjust.method="holm",sort.by="B",resort.by=NULL)
}
\arguments{
  \item{fit}{list containing a linear model fit produced by \code{lmFit}, \code{lm.series}, \code{gls.series} or \code{rlm.series}.
  For \code{topTable}, \code{fit} should be an object of class \code{MArrayLM} as produced by \code{lmFit}.}
  \item{coef}{column number or column name specifying which coefficient or contrast of the linear model is of interest}
  \item{number}{how many genes to pick out}
  \item{genelist}{data frame or character vector containing gene information.
  If not specified, this will be taken from the \code{genes} component of \code{fit}.}
  \item{A}{matrix of A-values or vector of average A-values.}
  \item{eb}{output list from \code{ebayes(fit)}}
  \item{adjust.method}{method to use to adjust the p-values for multiple testing.  Options are \code{"bonferroni"}, \code{"holm"}, \code{"hochberg"}, \code{"hommel"}, \code{"fdr"} and \code{"none"}. If \code{"none"} then the p-values are not adjusted. A \code{NULL} value will result in the default adjustment method, which is \code{"holm"}.}
  \item{sort.by}{character string specifying statistic to rank genes by.  Possibilities are \code{"M"}, \code{"A"}, \code{"T"}, \code{"P"} or \code{"B"}.}
  \item{resort.by}{character string specifying statistic to sort the selected genes by in the output data.frame.  Possibilities are \code{"M"}, \code{"A"}, \code{"T"}, \code{"P"} or \code{"B"}.}
  \item{...}{any other arguments are passed to \code{ebayes} if \code{eb} is \code{NULL}}
}
\value{
  A dataframe with a row for the \code{number} top genes and the following columns:
  \item{genelist}{if genelist was included as input}
  \item{M}{estimate of the effect or the contrast, on the log2 scale}
  \item{t}{moderated t-statistic}
  \item{P.Value}{nominal P-value}
  \item{B}{log odds that the gene is differentially expressed}
}
\details{
This function summarizes a linear model fit object produced by \code{lmFit}, \code{lm.series}, \code{gls.series} or \code{rlm.series} by selecting the top-ranked genes for any given contrast.
\code{topTable()} assumes that the linear model fit has already been processed by \code{eBayes()}.

The p-values for the coefficient/contrast of interest are adjusted for multiple testing by a call to \code{\link[stats]{p.adjust}}.
The \code{"holm"} method is the default because it is conservative and valid for any type of dependence between the p-values.
In most microarray contexts however the less conservative Benjamini and Hochberg method \code{"fdr"} may be more suitable.
See \code{help("p.adjust")} for more information.
Note, if there is no good evidence for differential expression in the experiment, that it is quite possible for all the adjusted p-values to be large, even for all of them to be equal to one.
It is quite possible for all the adjusted p-values to be equal to one if the smallest p-value is no smaller than \code{1/ngenes} where \code{ngenes} is the number of genes with non-missing p-values.

The \code{sort.by} argument is used to select the top genes.
Normally the genes appear in order of selection in the output table.
If one wants the table to be in some order other than selection order, the \code{resort.by} argument may be used.
For example, \code{topTable(fit, sort.by="B", resort.by="M")} selects the top genes according to log-odds of differential expression and then orders the resulting genes by log-ratio.
\code{topTable(fit, sort.by="M", resort.by="M")} would select the genes by absolute log-ratio and then sort then by signed log-ratio from must positive to most negative.
}
\seealso{
\code{\link{ebayes}}, \code{\link[stats]{p.adjust}}, \code{\link{lm.series}}, \code{\link{gls.series}}, \code{\link{rlm.series}}.
}
\author{Gordon Smyth}
\examples{
#  Simulate gene expression data,
#  6 microarrays and 100 genes with first gene differentially expressed
M <- matrix(rnorm(100*6,sd=0.3),100,6)
M[1,1:3] <- M[1,1:3] + 2
#  Design matrix includes two treatments, one for first 3 and one for last 3 arrays
design <- cbind(First3Arrays=c(1,1,1,0,0,0),Last3Arrays=c(0,0,0,1,1,1))
fit <- lmFit(M, design=design)
fit <- eBayes(fit)
topTable(fit)
}
\keyword{htest}
