\name{bwss}
\alias{bwss}
\title{Between and within sums of squares}
\description{Sums of squares between and within groups. Allows for missing values.}
\usage{bwss(x,group)}
\arguments{
  \item{x}{a numeric vector giving the responses.}
  \item{group}{a vector or factor giving the grouping variable.}
}
\value{
  A list with components
  \item{bss}{sums of squares between the group means.}
  \item{wss}{sums of squares within the groups.}
  \item{bdf}{degrees of freedom corresponding to \code{bss}.}
  \item{wdf}{degrees of freedom corresponding to \code{wss}.}
}
\details{This is equivalent to one-way analysis of variance.}
\author{Gordon Smyth}
\seealso{\code{\link{bwss.matrix}}}
\keyword{models}
