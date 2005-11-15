\name{lm.series}
\alias{lm.series}
\title{Fit Linear Model to Microrray Data by Ordinary Least Squares}
\description{Fit a linear model genewise to expression data from a series of arrays.
This function uses ordinary least squares and is a utility function for \code{lmFit}.}
\usage{lm.series(M,design=NULL,ndups=1,spacing=1,weights=NULL)}
\arguments{
  \item{M}{numeric matrix containing log-ratio or log-expression values for a series of microarrays, rows correspond to genes and columns to arrays}
  \item{design}{numeric design matrix defining the linear model. The number of rows should agree with the number of columns of M. The number of columns will determine the number of coefficients estimated for each gene.}
  \item{ndups}{number of duplicate spots. Each gene is printed ndups times in adjacent spots on each array.}
  \item{spacing}{the spacing between the rows of \code{M} corresponding to duplicate spots, \code{spacing=1} for consecutive spots}
  \item{weights}{an optional numeric matrix of the same dimension as \code{M} containing weights for each spot. If it is of different dimension to \code{M}, it will be filled out to the same size.}
}
\value{
  A list with components
  \item{coefficients}{numeric matrix containing the estimated coefficients for each linear model. Same number of rows as \code{M}, same number of columns as \code{design}.}
  \item{stdev.unscaled}{numeric matrix conformal with \code{coef} containing the unscaled standard deviations for the coefficient estimators. The standard errors are given by \code{stdev.unscaled * sigma}.}
  \item{sigma}{numeric vector containing the residual standard deviation for each gene.}
  \item{df.residual}{numeric vector giving the degrees of freedom corresponding to \code{sigma}.}
  \item{qr}{QR-decomposition of \code{design}}
}
\details{
This is a utility function used by the higher level function \code{\link{lmFit}}.
Most users should not use this function directly but should use \code{\link{lmFit}} instead.

The linear model is fit for each gene by calling the function \code{lm.fit} or \code{lm.wfit} from the base library.
}
\author{Gordon Smyth}
\seealso{
\code{\link[base:lmfit]{lm.fit}}.

An overview of linear model functions in limma is given by \link{06.LinearModels}.
}
\keyword{models}
\keyword{regression}
