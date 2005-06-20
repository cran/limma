\name{auROC}
\alias{auROC}
\title{Area Under Receiver Operating Curve}
\description{
Compute exact area under the ROC for empirical data.
}
\usage{
auROC(truth, stat)
}
\arguments{
  \item{truth}{numeric vector of 0 and 1 indicating whether the null or alternative respectively is true for each case.
  If \code{stat} is missing then \code{truth} is assuming to be already sorted in decreasing test statistic order.}
  \item{stat}{numeric vector containing test statistics.
  Hypotheses are to be rejected if \code{stat} exceeds a given threshold.}
}
\details{
This function computes the exact area under an empirical ROC curve.
The number of true and false discoveries are determined by how well the true states represented by \code{truth} match up with the observed statistics given by \code{stat}.
}
\value{
Numeric vector giving area under the curve, 1 being perfect and 0 being the minimum, or \code{NULL} if \code{truth} has zero length.
}
\seealso{
See \link{08.Tests} for other functions for testing and processing p-values.

See also \code{\link[ROC]{AUC}} in the \code{ROC} package.
}
\examples{
auROC(c(1,1,0,0,0))
truth <- rbinom(30,size=1,prob=0.2)
stat <- rchisq(30,df=2)
auROC(truth,stat)
}
\author{Gordon Smyth}
\keyword{htest}
