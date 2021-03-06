\name{toptable}
\alias{toptable}
\alias{topTable}
\alias{topTableF}
\title{Table of Top Genes from Linear Model Fit}
\description{
Extract a table of the top-ranked genes from a linear model fit.
}
\usage{
topTable(fit,coef=NULL,number=10,genelist=fit$genes,adjust.method="BH",sort.by="B",resort.by=NULL)
toptable(fit,coef=1,number=10,genelist=NULL,A=NULL,eb=NULL,adjust.method="BH",sort.by="B",resort.by=NULL,...)
topTableF(fit,number=10,genelist=fit$genes,adjust.method="BH")
}
\arguments{
  \item{fit}{list containing a linear model fit produced by \code{lmFit}, \code{lm.series}, \code{gls.series} or \code{mrlm}.
  For \code{topTable}, \code{fit} should be an object of class \code{MArrayLM} as produced by \code{lmFit} and \code{eBayes}.}
  \item{coef}{column number or column name specifying which coefficient or contrast of the linear model is of interest. Can also be a vector of column subscripts, in which case the gene ranking is by F-statistic for that set of contrasts.}
  \item{number}{how many genes to pick out}
  \item{genelist}{data frame or character vector containing gene information.
  For \code{topTable} only, this defaults to \code{fit$genes}.}
  \item{A}{matrix of A-values or vector of average A-values.
  For \code{topTable} only, this defaults to \code{fit$Amean}.}
  \item{eb}{output list from \code{ebayes(fit)}.
  If \code{NULL}, this will be automatically generated.}
  \item{adjust.method}{method used to adjust the p-values for multiple testing.  Options, in increasing conservatism, include \code{"none"}, \code{"BH"}, \code{"BY"} and \code{"holm"}.
  See \code{\link{p.adjust}} for the complete list of options. A \code{NULL} value will result in the default adjustment method, which is \code{"BH"}.}
  \item{sort.by}{character string specifying statistic to rank genes by.  Possibilities are \code{"logFC"}, \code{"A"}, \code{"T"}, \code{"t"}, \code{"P"}, \code{"p"} or \code{"B"}. \code{"M"} is allowed as a synonym for \code{"logFC"} for backward compatibility.}
  \item{resort.by}{character string specifying statistic to sort the selected genes by in the output data.frame.  Possibilities are \code{"logFC"}, \code{"A"}, \code{"T"}, \code{"t"}, \code{"P"}, \code{"p"} or \code{"B"}. \code{"M"} is allowed as a synonym for \code{"logFC"} for backward compatibility.}
  \item{...}{any other arguments are passed to \code{ebayes} if \code{eb} is \code{NULL}}
}
\value{
  A dataframe with a row for the \code{number} top genes and the following columns:
  \item{genelist}{if genelist was included as input}
  \item{logFC}{estimate of the log2-fold-change corresponding to the effect or contrast}
  \item{AveExpr}{average log2-expression for the probe over all arrays and channels, same as \code{Amean} in the \code{MarrayLM} object}
  \item{t}{moderated t-statistic}
  \item{P.Value}{raw p-value}
  \item{adj.P.Value}{adjusted p-value or q-value}
  \item{B}{log odds that the gene is differentially expressed}
}
\details{
Note that \code{toptable} is an earlier interface and is retained only for backward compatibility.

This function summarizes a linear model fit object produced by \code{lmFit}, \code{lm.series}, \code{gls.series} or \code{mrlm} by selecting the top-ranked genes for any given contrast.
\code{topTable()} assumes that the linear model fit has already been processed by \code{eBayes()}.

The p-values for the coefficient/contrast of interest are adjusted for multiple testing by a call to \code{\link[stats]{p.adjust}}.
The \code{"BH"} method, which controls the expected false discovery rate (FDR) below the specified value, is the default adjustment method because it is the most likely to be appropriate for microarray studies.
Note that the adjusted p-values from this method are bounds on the FDR rather than p-values in the usual sense.
Because they relate to FDRs rather than rejection probabilities, they are sometimes called q-values.
See \code{help("p.adjust")} for more information.

Note, if there is no good evidence for differential expression in the experiment, that it is quite possible for all the adjusted p-values to be large, even for all of them to be equal to one.
It is quite possible for all the adjusted p-values to be equal to one if the smallest p-value is no smaller than \code{1/ngenes} where \code{ngenes} is the number of genes with non-missing p-values.

The \code{sort.by} argument specifies the criterion used to select the top genes.
The choices are: \code{"logFC"} to sort by the (absolute) coefficient representing the log-fold-change; \code{"A"} to sort by average expression level (over all arrays) in descending order; \code{"T"} or \code{"t"} for absolute t-statistic; \code{"P"} or \code{"p"} for p-values; or \code{"B"} for the \code{lods} or B-statistic.

Normally the genes appear in order of selection in the output table.
If one wants the table to be in a different order, the \code{resort.by} argument may be used.
For example, \code{topTable(fit, sort.by="B", resort.by="logFC")} selects the top genes according to log-odds of differential expression and then orders the resulting genes by log-ratio in decreasing order.
Or \code{topTable(fit, sort.by="logFC", resort.by="logFC")} would select the genes by absolute log-ratio and then sort then by signed log-ratio from must positive to most negative.

\code{topTableF} ranks genes on the basis of the moderated F-statistic rather than t-statistics.
If \code{topTable} is called with \code{coef} has length greater than 1, then the specified columns will be extracted from \code{fit} and \code{topTableF} called on the result.
\code{topTable} with \code{coef=NULL} is the same as \code{topTableF}, unless the fitted model \code{fit} has only one column.
}
\note{
This is not the right function to use to create summary statistics for all the probes on an array.
Please consider using \code{write.fit} or \code{write} for this purpose, rather than using \code{topTable} with \code{number=nrow(fit)}.
}
\seealso{
An overview of linear model and testing functions is given in \link{06.LinearModels}.
See also \code{\link[stats]{p.adjust}} in the \code{stats} package.
}
\author{Gordon Smyth}
\examples{
#  See lmFit examples
}
\keyword{htest}
