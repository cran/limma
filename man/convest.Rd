\name{convest}
\alias{convest}
\title{Estimate Proportion of True Null Hypotheses}
\description{
Returns an estimate of the proportion of true null hypotheses using a convex decreasing density estimate on a vector of
p-values.
}
\usage{
convest(p, niter = 100, doplot = FALSE, doreport = FALSE)
}
\arguments{
  \item{p}{numeric vector of p-values, calculated using any test of your
    choice. Missing values are not allowed}
  \item{niter}{number of iterations to be used in fitting the convex,
    decreasing density for the p-values. Default is 100.}
  \item{doplot}{logical, should updated plots of fitted convex decreasing
    p-value density be produced at each iteration?
    Default is \code{FALSE}.}
  \item{doreport}{logical, should the estimated proportion be printed at each iteration?
  Default is \code{FALSE}.}
}
\details{
The proportion of true null hypotheses is often denoted pi0.
}
\value{
Numeric value in the interval [0,1] representing he estimated proportion of true null hypotheses.
}
\references{
Ferkingstad, E., Langaas, M., and Lindqvist, B. (2005).
Estimating the proportion of true null hypotheses, with application
to DNA microarray data.
\emph{Journal of the Royal Statistical Society Series} B, 67, 555-572.
Preprint at \url{http://www.math.ntnu.no/~mettela/SFG/research.imf}
}
\author{Egil Ferkingstad and Mette Langaas}
\seealso{
See \link{08.Tests} for other functions for producing or interpretting p-values.
}

\examples{
# First simulate data, use no.genes genes and no.ind individuals,
# with given value of pi0. Draw from normal distribution with mean=0
# (true null) and mean=mean.diff (false null).

no.genes <- 1000
no.ind <- 20
pi0 <- 0.9
mean.diff <- 1
n1 <- round(pi0*no.ind*no.genes)
n2 <- round((1-pi0)*no.ind*no.genes)
x <- matrix(c(rnorm(n1,mean=0),rnorm(n2,mean=mean.diff)),ncol=no.ind,byrow=TRUE)

# calculate p-values using your favorite method, e.g.
pvals <- ebayes(lm.series(x))$p.value

# run the convex decreasing density estimator to estimate pi0
convest(pvals,niter=100,doplot=interactive())
}
\keyword{htest}

