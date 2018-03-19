\name{cvtree}
\alias{cvtree}
\title{Extracts variable cluster IDs}
\description{Extracts cluster IDs for variables according to a dendrogram from object of class \code{\link{cvtree}}.}
\usage{
cvtree(object, k = 2, mincor = NULL, ...)
}
\arguments{
  \item{object}{Object of class \code{\link{corclust}}.}
  \item{k}{Number of clusters to be extracted from dendrogram.}
  \item{mincor}{Minimum within cluster correlation. Can be specified alternatively to \code{k}.}
  \item{...}{Currently not used.}  
}
\details{
Like in \code{\link{corclust}} for correlation comparison numerics and factors are considered separately. 
For factors Cramer's V statistic is used.
}
\value{Object of class \code{cvtree} with elements:
  \item{cluster}{Vector of cluster IDs.}
  \item{correlations}{Matrix of average within cluster correlations and average corrleation to all variables of the closest cluster as well as the ID of the closest cluster. For factor variables Cramer's V is computed.}
}

\references{Roever, C. and Szepannek, G. (2005): Application of a genetic algorithm
to variable selection in fuzzy clustering. In C. Weihs and W. Gaul (eds), Classification 
- The Ubiquitous Challenge, 674-681, Springer.}

\author{Gero Szepannek}

\seealso{See also \code{\link{corclust}},  \code{\link{plot.corclust}} and \code{\link{hclust}} for details on the clustering algorithm.}
\examples{
    data(B3)
    ccres <- corclust(B3)
    plot(ccres)
    cvtree(ccres, k = 3)
}

\keyword{cluster}
\keyword{multivariate}
\keyword{classif}
\keyword{manip}
\keyword{attribute}
