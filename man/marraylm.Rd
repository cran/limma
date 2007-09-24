\name{MArrayLM-class}
\docType{class}
\alias{MArrayLM-class}
\title{Microarray Linear Model Fit - class}

\description{
A list-based class for storing the results of fitting gene-wise linear models to a batch of microarrays.
Objects are normally created by \code{\link{lmFit}}.
}

\section{Slots/Components}{
\code{MArrayLM} objects do not contain any slots (apart from \code{.Data}) but they should contain the following list components:
  \describe{
    \item{\code{coefficients}:}{\code{matrix} containing fitted coefficients or contrasts}
    \item{\code{stdev.unscaled}:}{\code{matrix} containing unscaled standard deviations of the coefficients or contrasts}
    \item{\code{sigma}:}{\code{numeric} vector containing residual standard deviations for each gene}
    \item{\code{df.residual}:}{\code{numeric} vector containing residual degrees of freedom for each gene}
  }
  Objects may also contain the following optional components:
  \describe{
    \item{\code{Amean}:}{\code{numeric} vector containing the average log-intensity for each probe over all the arrays in the original linear model fit.
    Note this vector does not change when a contrast is applied to the fit using \code{contrasts.fit}.}
    \item{\code{genes}:}{\code{data.frame} containing gene names and annotation}
    \item{\code{design}:}{design \code{matrix} of full column rank}
    \item{\code{contrasts}:}{\code{matrix} defining contrasts of coefficients for which results are desired}
    \item{\code{F}:}{\code{numeric} vector giving moderated F-statistics for testing all contrasts equal to zero}
    \item{\code{F.p.value}:}{\code{numeric} vector giving p-value corresponding to \code{F.stat}}
    \item{\code{s2.prior}:}{\code{numeric} value giving empirical Bayes estimated prior value for residual variances}
    \item{\code{df.prior}:}{\code{numeric} vector giving empirical Bayes estimated degrees of freedom associated with \code{s2.prior} for each gene}
    \item{\code{s2.post}:}{\code{numeric} vector giving posterior residual variances}
    \item{\code{t}:}{\code{matrix} containing empirical Bayes t-statistics}
    \item{\code{var.prior}:}{\code{numeric} vector giving empirical Bayes estimated variance for each true coefficient}
  }
}

\section{Methods}{
\code{RGList} objects will return dimensions and hence functions such as \code{\link[limma:dim]{dim}}, \code{\link[base:nrow]{nrow}} and \code{\link[base:nrow]{ncol}} are defined. 
\code{MArrayLM} objects inherit a \code{show} method from the virtual class \code{LargeDataObject}.

The functions \code{\link{ebayes}} and \code{\link{classifyTestsF}} accept \code{MArrayLM} objects as arguments.
}

\author{Gordon Smyth}

\seealso{
  \link{02.Classes} gives an overview of all the classes defined by this package.
}

\keyword{classes}
\keyword{regression}
