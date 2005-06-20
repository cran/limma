\title{mdplot}
\name{mdplot}
\alias{mdplot}
\description{
Creates a mean-difference plot.
}
\usage{
mdplot(x, ...)
}
\arguments{
  \item{x}{numeric \code{matrix} with at least two columns}
  \item{...}{any other arguments are passed to \code{plot}}
}

\details{
Plots differences vs means for a set of bivariate values.
This is useful to contrast expression values for two microarrays.

Note that an MA-plot \code{\link{plotMA}} is a type of mean-difference plot.
}

\value{A plot is created on the current graphics device.}

\references{Chambers, J. M., Cleveland, W. S., Kleiner, B., and Tukey, P. A. (1983). Graphical Methods of Data Analysis. Wadsworth (pp. 48-57).

Cleveland, W. S., (1993). Visualizing Data. Hobart Press.

Bland, J. M., and Altman, D. G. (1986). Statistical methods for assessing agreement between two methods of clinical measurement. Lancet i, 307-310.

See also \url{http://www.statsci.org/micrarra/refs/maplots.html}
}
\author{Gordon Smyth}
\seealso{
An overview of diagnostic functions available in LIMMA is given in \link{09.Diagnostics}.
}
\keyword{hplot}
