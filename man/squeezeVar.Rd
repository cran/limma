\name{squeezeVar}
\alias{squeezeVar}

\title{Smooth Sample Variances}
\description{
Smooth a set of sample variances by computing empirical Bayes posterior means.
}

\usage{
squeezeVar(var, df)
}

\arguments{
\item{var}{numeric vector of independent sample variances}
\item{df}{numeric vector of degrees of freedom for the sample variances}
}

\details{
The sample variances \code{var} are assumed to follow scaled chi-squared distributions.
An inverse chi-squared prior is assumed for the true variances.
The scale and degrees of freedom for the prior distribution are estimated from the data.

The effect of this function is to smooth or shrink the variances towards a common value.
The smoothed variances have a smaller expected mean square error to the true variances than do the sample variances themselves.
}

\value{
A list with components
\item{var.post}{numeric vector of posterior variances}
\item{var.prior}{location of prior distribution}
\item{df.prior}{degrees of freedom of prior distribution}
}

\author{Gordon Smyth}

\references{
Smyth, G. K. (2004). Linear models and empirical Bayes methods for
     assessing differential expression in microarray experiments.
     \emph{Statistical Applications in Genetics and Molecular Biology}, \bold{3},
     No. 1, Article 3. \url{http://www.bepress.com/sagmb/vol3/iss1/art3}
}

\examples{
s2 <- rchisq(20,df=5)/5
squeezeVar(s2, df=5)
}

\keyword{htest}
