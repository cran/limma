\name{loessFit}
\alias{loessFit}
\title{Fast Simple Loess}
\description{
A fast version of locally weighted regression when there is only one x-variable and only the fitted values and residuals are required.
}
\usage{
loessFit(y, x, weights=NULL, span=0.3, bin=0.01/(2-is.null(weights)), iterations=4)
}
\arguments{
  \item{y}{numeric vector of response values.  Missing values are allowed.}
  \item{x}{numeric vector of predictor values  Missing values are allowed.}
  \item{weights}{numeric vector of non-negative weights.  Missing values are allowed.}
  \item{span}{numeric parameter between 0 and 1 specifying proportion of data to be used in the local regression moving window.
  Larger numbers give smoother fits.}
  \item{bin}{numeric value between 0 and 1 giving the proportion of the data which can be grouped in a single bin when doing local regression fit.
  \code{bin=0} forces an exact local regression fit with no interpolation.}
  \item{iterations}{number of iterations of loess fit}
}

\details{
This is a wrapper function to the Fortran and C code in the stats package which underlies the \code{lowess} and \code{loess} functions.
Its is to give a streamlined common interface to \code{lowess} and \code{loess} for use in \code{\link{normalizeWithinArrays}}.
When \code{weights} is null, this function is in effect a call to \code{lowess} in the stats package, with appropropriate choice of tuning parameters.
When \code{weights} is non-null, it is in effect a call to \code{loess}.
See the help pages for those functions for references and credits.

Note tha \code{lowess} is faster, needs less memory and is able to use a more accurate interpolation scheme than \code{loess}, so it is desirable to use \code{lowess} whenever \code{loess} is not needed to handle quantitative weights.

The arguments \code{span}, \code{cell} and \code{iterations} here have the same meaning as in \code{loess}.
\code{span} is equivalent to the argument \code{f} of \code{lowess} and \code{iterations} is equivalent to \code{iter+1}.
The parameter \code{bin} is intended to give a simple uniform interface to the \code{delta} argument of \code{lowess} and the \code{cell} argument of \code{loess}.
\code{bin} translates to \code{delta=bin*diff(range(x))} in a call to \code{lowess} or to \code{cell=bin/span} in a call to \code{loess}.

Unlike \code{lowess}, \code{loessFit} returns values in original rather than sorted order.
Also unlike \code{lowess}, \code{loessFit} allows missing values, the treatment being analogous to \code{na.exclude}.
Unlike \code{loess}, \code{loessFit} returns a linear regression fit if there are insufficient observations to estimate the loess curve.
}
\value{
A list with components
\item{fitted}{numeric vector of same length as \code{y} giving the loess fit}
\item{residuals}{numeric vector of same length as \code{x} giving residuals from the fit}
}

\author{Gordon Smyth, based on code from \code{lowess} and \code{loess} by BD Ripley}

\seealso{
See \code{\link[stats]{lowess}} and \code{\link[stats]{loess}} in the stats package.

See \link{05.Normalization} for an outline of the limma package normalization functions.
}

\examples{
y <- rnorm(1000)
x <- rnorm(1000)
w <- rep(1,1000)
# The following are equivalent apart from execution time
# and interpolation inaccuracies
system.time(fit <- loessFit(y,x)$fitted)
system.time(fit <- loessFit(y,x,w)$fitted)
system.time(fit <- fitted(loess(y~x,weights=w,span=0.3,family="symmetric",iterations=4)))
# The same but with sorted x-values
system.time(fit <- lowess(x,y,f=0.3)$y)
}

\keyword{models}
