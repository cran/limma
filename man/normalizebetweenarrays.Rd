\name{normalizeBetweenArrays}
\alias{normalizeBetweenArrays}
\title{Normalize Between Arrays}

\description{
Normalizes expression intensities so that the intensities or log-ratios have similar distributions across a series of arrays.
}

\usage{
normalizeBetweenArrays(object, method="Aquantile", targets=NULL, ...)
}

\arguments{
  \item{object}{a \code{matrix}, \code{\link[limma:rglist]{RGList}} or \code{\link[limma:malist]{MAList}} object containing expression ratios for a series of arrays}
  \item{method}{character string specifying the normalization method to be used.
  Choices are \code{"none"}, \code{"scale"}, \code{"quantile"}, \code{"Aquantile"}, \code{"Gquantile"}, \code{"Rquantile"}, \code{"Tquantile"} or \code{"vsn"}.
  A partial string sufficient to uniquely identify the choice is permitted.}
  \item{targets}{vector, factor or matrix of length twice the number of arrays, used to indicate target groups if \code{method="Tquantile"}}
  \item{...}{other arguments are passed to \code{normalizeQuantiles} if one of the quantile methods are used or to \code{vsn} if \code{method="vsn"}}
}

\details{
\code{normalizeWithinArrays} normalizes expression values to make intensities consistent within each array.
\code{normalizeBetweenArrays} normalizes expression values to achieve consistency between arrays.
Normalization between arrays is usually, but not necessarily, applied after normalization within arrays.
An exception is \code{method="vsn"}, see below.

The scale normalization method was proposed by Yang et al (2001, 2002) and is further explained by Smyth and Speed (2003).
The idea is simply to scale the log-ratios to have the same median-abolute-deviation (MAD) across arrays.
This idea has also been implemented by the \code{maNormScale} function in the marrayNorm package.
The implementation here is slightly different in that the MAD scale estimator is replaced with the median-absolute-value and the A-values are normalized as well as the M-values.

Quantile normalization was proposed by Bolstad et al (2003) for Affymetrix-style single-channel arrays and by Yang and Thorne (2003) for two-color cDNA arrays.
\code{method="quantile"} ensures that the intensities have the same empirical distribution across arrays and across channels.
\code{method="Aquantile"} ensures that the A-values (average intensities) have the same empirical distribution across arrays leaving the M-values (log-ratios) unchanged.
These two methods are called "q" and "Aq" respectively in Yang and Thorne (2003).

\code{method="Tquantile"} performs quantile normalization separately for the groups indicated by \code{targets}.
\code{targets} may be a target matrix such as read by \code{readTargets} or can be a vector indicating green channel groups followed by red channel groups.

\code{method="Gquantile"} ensures that the green (first) channel has the same empirical distribution across arrays, leaving the M-values (log-ratios) unchanged.
This method might be used when the green channel is a common reference throughout the experiment.
In such a case the green channel represents the same target throughout, so it makes compelling sense to force the distribution of intensities to be same for the green channel on all the arrays, and to adjust to the red channel accordingly.
\code{method="Rquantile"} ensures that the red (second) channel has the same empirical distribution across arrays, leaving the M-values (log-ratios) unchanged.
Both \code{Gquantile} and \code{Rquantile} normalization have the implicit effect of changing the red and green log-intensities by equal amounts.

If \code{object} is a \code{matrix} then the scale, quantile or vsn normalization will be applied to the columns.
Applying \code{method="Aquantile"} when \code{object} is a \code{matrix} will produce an error.

\code{method="vsn"} uses the \code{\link[vsn]{vsn}} function from the vsn package.
For this option the input \code{object} should contain raw intensities, i.e., prior to background correction, log-transformation or any normalization.
Note that the normalized intensities are on the log-2 scale, not the log-e scale output by the \code{vsn} function in the vsn package.}

\value{
If \code{object} is a matrix then \code{normalizeBetweenArrays} produces a matrix of the same size.
Otherwise, \code{normalizeBetweenArrays} produces an \code{\link[limma:malist]{MAList}} object with M and A-values on the log-2 scale.
}

\author{Gordon Smyth}

\references{
Bolstad, B. M., Irizarry R. A., Astrand, M., and Speed, T. P. (2003), A comparison of normalization methods for high density oligonucleotide array data based on bias and variance. \emph{Bioinformatics} \bold{19}, 185-193.

Smyth, G. K., and Speed, T. P. (2003). Normalization of cDNA microarray data. In: \emph{METHODS: Selecting Candidate Genes from DNA Array Screens: Application to Neuroscience}, D. Carter (ed.). Methods Volume 31, Issue 4, December 2003, pages 265-273.

Yang, Y. H., Dudoit, S., Luu, P., and Speed, T. P. (2001). Normalization for cDNA microarray data. In \emph{Microarrays: Optical Technologies and Informatics}, M. L. Bittner, Y. Chen, A. N. Dorsel, and E. R. Dougherty (eds), Proceedings of SPIE, Volume 4266, pp. 141-152. 

Yang, Y. H., Dudoit, S., Luu, P., Lin, D. M., Peng, V., Ngai, J., and Speed, T. P. (2002). Normalization for cDNA microarray data: a robust composite method addressing single and multiple slide systematic variation. \emph{Nucleic Acids Research} \bold{30}(4):e15.

Yang, Y. H., and Thorne, N. P. (2003). Normalization for two-color cDNA microarray data.
In: D. R. Goldstein (ed.), \emph{Science and Statistics: A Festschrift for Terry Speed}, IMS Lecture Notes - Monograph Series, Volume 40, pp. 403-418.
}

\seealso{
  An overview of LIMMA functions for normalization is given in \link{05.Normalization}.

  See also \code{\link[marrayNorm]{maNormScale}} in the marrayNorm package, \code{\link[affy:normalize-methods]{normalize}} in the affy package and \code{\link[vsn]{vsn}} in the vsn package.
}

\examples{
library(sma)
data(MouseArray)
MA <- normalizeWithinArrays(mouse.data, mouse.setup)
plot.scale.box(MA$M)

#  Between array scale normalization as in Yang et al (2001):
MA <- normalizeBetweenArrays(MA,method="scale")
print(MA)
show(MA)
plot.scale.box(MA$M)

#  One can get the same results using the matrix method:
M <- normalizeBetweenArrays(MA$M,method="scale")
plot.scale.box(M)

#  MpAq normalization as in Yang and Thorne (2003):
MpAq <- normalizeWithinArrays(mouse.data, mouse.setup)
MpAq <- normalizeBetweenArrays(MpAq, method="Aq")
plotDensities(MpAq)
}

\keyword{models}
\keyword{multivariate}
