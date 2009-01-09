\name{dkernel}
\alias{dkernel}
\title{Estimate density of a given kernel}
\description{
Given an estimated kernel density this function estimates the density of a new vector.
}
\usage{
dkernel(x, kernel = density(x), interpolate = FALSE, ...)
}
\arguments{
  \item{x}{vector of which the density should be estimated}
  \item{kernel}{object of \code{class} \code{\link{density}}}
  \item{interpolate}{Interpolate or use \code{\link{density}} of nearest point?}
  \item{\dots}{currently not used.}
}
\value{
    Denstiy of \code{x} in \code{kernel}.
}
\author{Karsten Luebke}
\seealso{\code{\link{density}}, \code{\link{NaiveBayes}}}
\examples{
kern <- density(rnorm(50))
x <- seq(-3, 3, len = 100)
y <- dkernel(x, kern)
plot(x, y, type = "l")
}
\keyword{distribution}
\keyword{nonparametric}
\concept{Kernel density estimation}
