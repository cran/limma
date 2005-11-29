\name{residuals.MArrayLM}
\title{Extract Residuals from MArrayLM Fit}
\usage{
\S3method{residuals}{MArrayLM}(object, y, \dots)
}
\alias{residuals.MArrayLM}
\arguments{
  \item{object}{a fitted model object inheriting from class \code{MarrayLM}.}
  \item{y}{a data object containing the response data used to compute the fit.
  This can be of any class for which \code{as.matrix} is defined, including \code{MAList}, \code{exprSet}, \code{marrayNorm} etc.}
  \item{\dots}{other arguments are not used}
}
\description{
This method extracts the residuals from all the probewise linear model fits
and returns them in a matrix.
}
\value{
  Numeric matrix of residuals.
}
\seealso{
  \code{\link[stats]{residuals}}.
}
\keyword{models}
\keyword{regression}
