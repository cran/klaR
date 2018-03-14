\name{kmodes}
\alias{kmodes}
\alias{print.kmodes}
\title{
K-Modes Clustering
}
\description{
  Perform k-modes clustering on categorical data.
}
\usage{
kmodes(data, modes, iter.max = 10, weighted = FALSE, fast = TRUE)
}
\arguments{
  \item{data}{A matrix or data frame of categorical data. Objects have to be in rows, variables in columns.}
  \item{modes}{Either the number of modes or a set of initial
    (distinct) cluster modes.  If a number, a random set of (distinct)
    rows in \code{data} is chosen as the initial modes.}
  \item{iter.max}{The maximum number of iterations allowed.}
  \item{weighted}{Whether usual simple-matching distance between objects is used, or a weighted version of this distance.}
  \item{fast}{Logical Whether a fast version of the algorithm should be applied.}

}
\details{
    The \eqn{k}-modes algorithm (Huang, 1997) an extension of the k-means algorithm by MacQueen (1967).

  The data given by \code{data} is clustered by the \eqn{k}-modes method (Huang, 1997)
  which aims to partition the objects into \eqn{k} groups such that the
  distance from objects to the assigned cluster modes is minimized.
  
  By default simple-matching distance is used to determine the dissimilarity of two objects. It is computed by counting the number of mismatches in all variables.
  Alternative this distance is weighted by the frequencies of the categories in data (see Huang, 1997, for details).
 
  If an initial matrix of modes is supplied, it is possible that
  no object will be closest to one or more modes. In this case less cluster than supplied modes will be returned 
  and a warning is given.
  
  If called using \code{fast = TRUE} the reassignment of the data to clusters is done for the entire data set before recomputation of the modes is done. For computational reasons this option should be chosen unless moderate data sizes.   
  
  For clustering mixed type data it is referred to \code{\link[clustMixType]{kproto}}.

}
\value{
  An object of class \code{"kmodes"} which is a list with components:

  \item{cluster}{A vector of integers indicating the cluster to which each object is allocated.} 
  \item{size}{The number of objects in each cluster.}
  \item{modes}{A matrix of cluster modes.}
  \item{withindiff}{The within-cluster simple-matching distance for each cluster.}
  \item{iterations}{The number of iterations the algorithm has run.}
  \item{weighted}{Whether weighted distances were used or not.}
}
\author{Christian Neumann, \email{christian2.neumann@tu-dortmund.de}, Gero Szepannek, \email{gero.szepannek@web.de}}
\references{
  Huang, Z. (1997)  A Fast Clustering Algorithm to Cluster Very Large Categorical Data Sets in Data Mining.
  in \emph{KDD: Techniques and Applications} (H. Lu, H. Motoda and H. Luu, Eds.), pp. 21-34, World Scientific, Singapore.

  MacQueen, J. (1967)  Some methods for classification and analysis of
  multivariate observations. In \emph{Proceedings of the Fifth Berkeley
    Symposium on  Mathematical Statistics and  Probability},
  eds L. M. Le Cam \& J. Neyman,
  \bold{1}, pp. 281-297. Berkeley, CA: University of California Press.
}
\examples{
### a 5-dimensional toy-example:

## generate data set with two groups of data:
set.seed(1)
x <- rbind(matrix(rbinom(250, 2, 0.25), ncol = 5),
           matrix(rbinom(250, 2, 0.75), ncol = 5))
colnames(x) <- c("a", "b", "c", "d", "e")

## run algorithm on x:
(cl <- kmodes(x, 2))

## and visualize with some jitter:
plot(jitter(x), col = cl$cluster)
points(cl$modes, col = 1:5, pch = 8)
}
\keyword{multivariate}
\keyword{cluster}
\keyword{category}

