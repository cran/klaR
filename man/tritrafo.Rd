\name{tritrafo}
\alias{tritrafo}
\encoding{latin1}
\title{Barycentric plots}
\description{
Function to carry out the transformation into 2D space
for \code{\link{triplot}}, \code{\link{trilines}} etc.
}
\usage{
tritrafo(x, y = NULL, z = NULL, check = TRUE, tolerance = 0.0001)
}
\arguments{
  \item{x}{Vector of fractions of first component
OR 3-column matrix containing all three components (omitting \code{y} and \code{z})
OR 3-element vector (for all three components, omitting \code{y} and \code{z}).}
  \item{y}{(optional) vector of fractions of second component.}
  \item{z}{(optional) vector of fractions of third component.}
  \item{check}{if \code{TRUE}, it is checked whether \code{x+y+z=1} and \code{x,y,z>=0} for all cases.}
  \item{tolerance}{tolerance for above sum check.}
}
\details{
Projects the mixture given by \code{x}, \code{y}, and \code{z} 
with \code{x}, \code{y}, \code{z} between one and zero and \code{x+y+z=1} into
a two-dimensional space.

For further details see \code{\link{triplot}}.
}
\value{
A matrix with two columns corresponding to the two dimensions.
}
\author{Christian Röver, \email{roever@statistik.uni-dortmund.de}}
\seealso{\code{\link{triplot}}, \code{\link{tripoints}}, \code{\link{trilines}}, \code{\link{trigrid}}}
\examples{
tritrafo(0.1, 0.2, 0.7)
tritrafo(0.1, 0.2, 0.6) # warning 

triplot()
points(tritrafo(0.1, 0.2, 0.7), col="red")
tripoints(0.1, 0.2, 0.7, col="green")  # the same

tritrafo(c(0.1,0.2), c(0.3,0.4), c(0.6,0.4))
tritrafo(diag(3))

point <- c(0.25,0.6,0.15)
triplot(point, pch=16)
text(tritrafo(point), "(0.25, 0.60, 0.15)", adj=c(0.5,2)) # add a label 
}
\keyword{dplot}
\concept{Visualizing Classification Performance Measures}
\concept{Barycentric plots}
