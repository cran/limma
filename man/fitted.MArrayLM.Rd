\name{fitted.MArrayLM}
\alias{fitted.MArrayLM}
\concept{regression}
\title{Fitted Values Method for MArrayLM Fits}
\usage{
\S3method{fitted}{MArrayLM}(object, design = object$design, \dots)
}
\arguments{
  \item{object}{a fitted object of class inheriting from \code{"MArrayLM"}.}
  \item{design}{numeric design matrix.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\description{
  Obtains fitted values from a fitted microarray linear model object.
}
\value{
  A numeric matrix of fitted values.
}
\seealso{
  \code{\link[stats]{fitted}}
}

\author{Gordon Smyth}
\keyword{models}
\keyword{regression}
