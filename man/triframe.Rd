\name{triframe}
\alias{triframe}
\encoding{latin1}
\title{Barycentric plots}
\description{Function to add a frame to an existing (barycentric) plot.}
\usage{
triframe(label = 1:3, label.col = 1, cex = 1, ...)
}
\arguments{
  \item{label}{labels for the three corners of the plot.}
  \item{label.col}{text color for labels.}
  \item{cex}{Magnification factor for label text relative to the default.}
  \item{\dots}{Further graphical parameters passed to \code{\link{trilines}}.}
}
\author{Christian Röver, \email{roever@statistik.uni-dortmund.de}}
\seealso{\code{\link{triplot}}, \code{\link{trilines}}, \code{\link{trigrid}}, \code{\link{centerlines}}}
\examples{
triplot(grid = TRUE, frame = FALSE)     # plot without frame 
some.triangle <- rbind(c(0, 0.65, 0.35), c(0.53, 0.47, 0), 
                       c(0.72, 0, 0.28))[c(1:3, 1), ]
trilines(some.triangle, col = "red", pch = 16, type = "b")
triframe(label = c("left", "top", "right"), col = "blue", 
         label.col = "green3")          # frame on top of points 
}
\keyword{aplot}
\concept{Visualizing Classification Performance Measures}
\concept{Barycentric plots}
