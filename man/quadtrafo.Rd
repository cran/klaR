\name{quadtrafo}
\alias{quadtrafo}
\alias{quadpoints}
\alias{quadlines}

\title{Transforming of 4 dimensional values in a barycentric coordinate system.}
\description{
Transforming of 4 dimensional values in a barycentric coordinate system.
}
\usage{
quadtrafo(e, f = NULL, g = NULL, h = NULL)
        
quadlines(e, f = NULL, g = NULL, h = NULL, sp, ...)

quadpoints(e, f = NULL, g = NULL, h = NULL, sp, ...)
}
\arguments{
  \item{e}{either a matrix with 4 columns represanting the membership values or a 
            vector with the membership values of the first class}
  \item{f}{vector with the membership values of the second class}
  \item{g}{vector with the membership values of the third class}
  \item{h}{vector with the membership values of the forth class}
  \item{sp}{\code{scatterplot3d} object to which \code{points} or \code{lines} should be plotted}
  \item{\dots}{further arguments passed to the underlyind \code{plot} functions}
}

\value{
For \code{quadtrafo} the 3 dimensional values in the barycentrix coordinate system. 
}

\author{Karsten Luebke, \email{karsten.luebke@fom.de}}

\seealso{\code{\link{triplot}}, \code{\link{quadplot}}}
\examples{
library(MASS)
data(B3)
posterior <- predict(lda(PHASEN ~ ., data = B3))$post
quadtrafo(posterior)
}
\keyword{internal}
\concept{Visualizing Classification Performance Measures}
\concept{Barycentric plots}
