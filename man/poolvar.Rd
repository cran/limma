\name{poolVar}
\alias{poolVar}

\title{Pool Sample Variances with Unequal Variances}
\description{
Compute the Satterthwaite (1946) approximation to the distribution of a weighted sum of sample variances.
}

\usage{
poolVar(var, df=n-1, multiplier=1/n, n)
}

\arguments{
\item{var}{numeric vector of independent sample variances}
\item{df}{numeric vector of degrees of freedom for the sample variances}
\item{multiplier}{numeric vector giving multipliers for the sample variances}
\item{n}{numeric vector of sample sizes}
}

\details{
The sample variances \code{var} are assumed to follow scaled chi-square distributions.
A scaled chi-square approximation is found for the distribution of \code{sum(multiplier * var)} by equating first and second moments.
On output the sum to be approximated is equal to \code{multiplier * var} which follows approximately a scaled chisquare distribution on \code{df} degrees of freedom.
The approximation was proposed by Satterthwaite (1946).

If there are only two groups and the degrees of freedom are one less than the sample sizes then this gives the denominator of Welch's t-test for unequal variances.
}

\value{
A list with components
\item{var}{effective pooled sample variance}
\item{df}{effective pooled degrees of freedom}
\item{multiplier}{pooled multiplier}
}

\author{Gordon Smyth}

\references{
Welch, B. L. (1938). The significance of the difference between two means when the population variances are unequal.
\emph{Biometrika} \bold{29}, 350-362.

Satterthwaite, F. E. (1946). An approximate distribution of estimates of variance components.
\emph{Biometrics Bulletin} \bold{2}, 110-114.

Welch, B. L. (1947). The generalization of 'Student's' problem when several different population variances are involved.
\emph{Biometrika} \bold{34}, 28-35.

Welch, B. L. (1949). Further note on Mrs. Aspin's tables and on certain approximations to the tabled function. \emph{Biometrika} \bold{36}, 293-296.
}

\seealso{
\link{10.Other}
}

\examples{
#  Welch's t-test with unequal variances
x <- rnorm(10,mean=1,sd=2)
y <- rnorm(20,mean=2,sd=1)
s2 <- c(var(x),var(y))
n <- c(10,20)
out <- poolVar(var=s2,n=n)
tstat <- (mean(x)-mean(y)) / sqrt(out$var*out$multiplier)
pvalue <- 2*pt(-abs(tstat),df=out$df)
#  Equivalent to t.test(x,y)
}

\keyword{htest}
