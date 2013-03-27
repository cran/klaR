\name{TopoS}
\alias{TopoS}

\title{Computation of criterion S of a visualization}
\description{
  \code{TopoS} computes one version of the criterion \code{stress} as i.e. 
    known from \code{\link[MASS]{sammon}} for a given visualization.
}
\usage{
TopoS(EV.dist, Cells.dist)
}

\arguments{
  \item{EV.dist}{a symmetric distance matrix consisting of distances in the origin space (\code{dx})}
  \item{Cells.dist}{a symmetric distance matrix consisting of distances in the target space (\code{dy})}
}
\details{
  TopoS first performs a homogeneous linear regression where \code{dx} is predicted by \code{dy}. 
  Then the residual sum of squares of this regression is computed and divided by the sum of squared \code{dx}. 
  Finally one minus the square root of the latter result defines \code{S}.
}
\value{
  \code{TopoS} returns a numeric between zero and one which is the criterion \code{S} of the given visualization.
}

\author{Nils Raabe}


\seealso{\code{\link{EDAM}}}
\examples{
# Compute S for the MDS visualization of the german business data
data(B3)
plot(cmdscale(dist(B3[, 2:14])), col = rainbow(4)[B3[, 1]], pch = 16)
TopoS(dist(B3[, 2:14]), dist(cmdscale(dist(B3[, 2:14]))))
}
\keyword{internal}
\concept{S criterion}
