\name{EDAM}
\alias{EDAM}

\title{Computation of an Eight Direction Arranged Map}
\description{
  Produces an object of class \code{EDAM} which is a two dimensional representation 
  of data in a rectangular, equally spaced 
    grid as known from Self-Organizing Maps.
}
\usage{
EDAM(EV0, nzx = 0, iter.max = 10, random = TRUE, standardize = FALSE, 
    wghts = 0, classes = 0, sa = TRUE, temp.in = 0.5, temp.fin = 1e-07, 
    temp.gamma = 0)
}

\arguments{
  \item{EV0}{ either a symmetric dissimilarity matrix or a matrix of arbitrary dimensions whose 
                n rows correspond to cases and whose k columns correspond to variables. }
  \item{nzx}{ an integer specifying the number of vertical bars in the grid. By default, 
                \code{nzx} is chosen automatically, so that the grid
                gets closest do a square. If n is no multiple of \code{nzx}, all surplus objects are skipped.}
  \item{iter.max}{an integer giving the maxmimum number of iterations to perform for the same neighborhood size.}
  \item{random}{logical. If \code{TRUE}, the initital order is drawn from a uniform distribution.}
  \item{standardize}{logical. If \code{TRUE}, the measurements in \code{EV0} are standardized before 
        calculating Euclidean distances. 
        Measurements are standardized for each variable by dividing by the variable's 
        standard deviation. Meaningless if \code{EV0} is a dissimilarity matrix. }
  \item{wghts}{ an optional vector of length k giving relative weights of the variables in 
        computing Euclidean distances. Meaningless if \code{EV0} is a dissimilarity matrix. }
  \item{classes}{ an optional vector of length n specifying the membership to classes for all objects.}
  \item{sa}{logical. If \code{TRUE}, the optimization is obtained by Simulated Annealing.}
  \item{temp.in}{numeric giving the initial temperature, if \code{sa} is set to \code{TRUE}.}
  \item{temp.fin}{ numeric giving the final temperature, if \code{sa} is set to \code{TRUE}. 
        Meaningless if \code{temp.gamma} is greater than 0.}
  \item{temp.gamma}{ numeric giving the relative change of the temperature from one iteration to the other, 
        if \code{sa} is set to \code{TRUE}.}
}
\details{
  The data given by \code{EV0} is visualized by the EDAM-algorithm. This method approximates the best visualization where
  goodness is measured by \code{S}, a transformation of the criterion \code{stress} as i.e. 
  known from \code{\link[MASS]{sammon}}. 
  The target space of the visualization is restricted to a grid so the problem has a discrete solution space. 
  Originally this restriction was made to make the results
  comparable to those of Kohonen Self-Organizing Maps. But it turns out that also for reasons of a clear arrangement the
  representation in a grid can be more favorable than in the hole plane.
  
  During the computation of EDAM 3 values indicating its progress are given online. The first is the number of the actual
  iteration, the second the maximum number of overall performed iterations. The latter may reduce during computation,
  since the neighborhood reduces in case of convergence before the last iteration. 
  The last number gives the actual criterion S.
  The default plot method \code{plot.edam} for objects of class \code{EDAM} is \code{\link{shardsplot}}.
}
\value{
  EDAM returns an object of \code{class} \code{EDAM}, which is a list containing the following components:
  \item{preimages}{the re-ordered data; the position of the i-th object is where \code{Z} equals i.}
  \item{Z}{a matrix representing the positions of the \code{preimages} in the grid by their numbers.}
  \item{Z.old.terms}{a matrix representing the positions of the data in original order in the grid by their numbers.}
  \item{cl.ord}{a vector giving the re-ordered classes. All elements equal 1 if argument \code{classes} is undefined.}
  \item{S}{the criterion of the map} 
}
\references{
  Raabe, N. (2003).
  \emph{Vergleich von Kohonen Self-Organizing-Maps mit einem nichtsimultanen 
  Klassifikations- und Visualisierungsverfahren}.
  Diploma Thesis, Department of Statistics, University of Dortmund.
  \url{http://www.statistik.tu-dortmund.de/de/content/einrichtungen/lehrstuehle/personen/raabe/Diplomarbeit.pdf}.
}
\author{Nils Raabe}


\seealso{
    \code{\link{shardsplot}}, \code{\link{TopoS}}
}
\examples{
# Compute an Eight Directions Arranged Map for a random sample 
# of the iris data.
data(iris)
set.seed(1234)
iris.sample <- sample(150, 42)
irisEDAM <- EDAM(iris[iris.sample, 1:4], classes = iris[iris.sample, 5], 
    standardize = TRUE, iter.max = 3)
plot(irisEDAM, vertices = FALSE)
legend(3, 5, col = rainbow(3), legend = levels(iris[,5]), pch = 16)
print(irisEDAM)

# Construct clusters within the phases of the german business data 
# and visualize the centroids by EDAM.
data(B3)
phasemat <- lapply(1:4, function(x) B3[B3[,1] == x, 2:14])
subclasses <- lapply(phasemat, 
    function(x) cutree(hclust(dist(x)), k = round(nrow(x) / 4.47)))
centroids <- lapply(1:4, 
    function(y) apply(phasemat[[y]], 2, 
        function(x) by(x, subclasses[[y]], mean)))
centmat <- matrix(unlist(sapply(centroids, t)), ncol = 13, 
    byrow = TRUE, dimnames = list(NULL, colnames(centroids[[1]])))
centclasses <- unlist(lapply(1:4, 
    function(x) rep(x, unlist(lapply(centroids, nrow))[x])))
B3EDAM <- EDAM(centmat, classes = centclasses, standardize = TRUE, 
    iter.max = 6, rand = FALSE)
plot(B3EDAM, standardize = TRUE)
opar <- par(xpd = NA)
legend(4, 5.1, col = rainbow(4), pch = 16, xjust = 0.5, yjust = 0,
    ncol = 2, legend = c("upswing", "upper turning point", 
                         "downswing", "lower turning point"))
print(B3EDAM)
par(opar)
}
\keyword{multivariate}
\concept{Eight Direction Arranged Maps}
\concept{EDAM}
\concept{Cluster analysis}
