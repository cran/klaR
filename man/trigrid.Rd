\name{trigrid}
\alias{trigrid}
\encoding{latin1}
\title{Barycentric plots}
\description{Function to add a grid to an existing (barycentric) plot.}
\usage{
trigrid(x = seq(0.1, 0.9, by = 0.1), y = NULL, z = NULL, 
    lty = "dashed", col = "grey", ...)
}

\arguments{
  \item{x}{Values along which to draw grid lines for first dimension 
    (or all dimensions if \code{y} and \code{z} omitted). 
  For NO grid lines in some dimensions just supply an \code{NA}.}
  \item{y}{Grid lines for second dimension.}
  \item{z}{Grid lines for third dimension.}
  \item{lty}{Line type (see \code{\link{par}}).}
  \item{col}{Line colour (see \code{\link{par}}).}
  \item{\dots}{Further graphical parameters passed to \code{\link{trilines}}.}
}
\details{
Grid lines illustrate the set of points for which one of the dimensions is held constant; 
e.g. horizontal lines contain all points with a certain value y for the second dimension, 
connecting the two extreme points (0,y,1-y) and (1-y,y,0).

Grids may be designed more flexible than with \code{triplot}'s \code{grid} option.}
\author{Christian Röver, \email{roever@statistik.uni-dortmund.de}}
\seealso{\code{\link{triplot}}, \code{\link{trilines}}, \code{\link{triframe}}, \code{\link{centerlines}}}
\examples{
triplot(grid = FALSE)
trigrid(c(1/3, 0.5)) # same grid for all 3 dimensions 

triplot(grid = c(1/3, 0.5))  # (same effect) 

triplot(grid = FALSE)
# different grids for all dimensions:
trigrid(x = 1/3, y = 0.5, z = seq(0.2, 0.8, by=0.2))  

triplot(grid = FALSE)
# grid for third dimension only:
trigrid(x = NA, y = NA, z = c(0.1, 0.2, 0.4, 0.8))  
}

\keyword{aplot}
\concept{Visualizing Classification Performance Measures}
\concept{Barycentric plots}
