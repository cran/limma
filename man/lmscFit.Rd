\name{lmscFit}
\alias{lmscFit}
\title{Fit Linear Model to Individual Channels of Two-Color Data}
\description{Fit a linear model to the individual log-intensities for each gene given a series of two-color arrays}
\usage{
lmscFit(object, design, correlation)
}
\arguments{
  \item{object}{an \code{\link[limma:malist]{MAList}} object or a list from which \code{M} and \code{A} values may be extracted}
  \item{design}{a numeric matrix containing the design matrix for linear model in terms of the individual channels.
  The number of rows should be twice the number of arrays.
  The number of columns will determine the number of coefficients estimated for each gene.}
  \item{correlation}{numeric value giving the intra-spot correlation}
}
\value{
An object of class \code{\link[limma:marraylm]{MArrayLM}}
}
\details{
For two color arrays, the channels measured on the same set of arrays are correlated.
The \code{M} and \code{A} however are uncorrelated for each gene.
This function fits a linear model to the set of M and A-values for each gene after re-scaling the M and A-values to have equal variances.
The input correlation determines the scaling required.
The input correlation is usually estimated using \code{\link{intraspotCorrelation}} before using \code{lmscFit}.

Missing values in \code{M} or \code{A} are not allowed.
}
\author{Gordon Smyth}
\seealso{
\code{\link[stats:lmfit]{lm.fit}}.

An overview of methods for single channel analysis in limma is given by \link{07.SingleChannel}.
}
\references{
Smyth, G. K. (2005). Individual channel analysis of two-colour microarray data. Invited Session IPM 11: Computational Tools For Microarray Analysis, 55th Session of the International Statistics Institute, Sydney, 12 April 2005. (Four-page paper distributed on CD at the conference.)
}
\examples{
library(sma)
# Subset of data from ApoAI case study in Limma User's Guide
data(MouseArray)
# Avoid non-positive intensities
RG <- backgroundCorrect(mouse.data,method="normexp")
MA <- normalizeWithinArrays(RG,mouse.setup)
MA <- normalizeBetweenArrays(MA,method="Aq")
# Randomly choose 500 genes for this example
i <- sample(1:nrow(MA),500)
MA <- MA[i,]
targets <- data.frame(Cy3=I(rep("Pool",6)),Cy5=I(c("WT","WT","WT","KO","KO","KO")))
targets.sc <- targetsA2C(targets)
targets.sc$Target <- factor(targets.sc$Target,levels=c("Pool","WT","KO"))
design <- model.matrix(~Target,data=targets.sc)
corfit <- intraspotCorrelation(MA,design)
fit <- lmscFit(MA,design,correlation=corfit$consensus)
cont.matrix <- cbind(KOvsWT=c(0,-1,1))
fit2 <- contrasts.fit(fit,cont.matrix)
fit2 <- eBayes(fit2)
topTable(fit2,adjust="fdr")
}
\keyword{models}
\keyword{regression}
