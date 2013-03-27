\name{quadplot}
\alias{quadplot}
\title{Plotting of 4 dimensional membership representation simplex}
\description{
For a 4 class discrimination problem the membership values of each class are visualized in
a 3 dimensional barycentric coordinate system.
}
\usage{
quadplot(e = NULL, f = NULL, g = NULL, h = NULL, angle = 75, 
    scale.y = 0.6, label = 1:4, labelcol = rainbow(4), 
    labelpch = 19, labelcex = 1.5, main = "", s3d.control = list(), 
    simplex.control = list(), legend.control = list(), ...)
}
\arguments{
  \item{e}{either a matrix with 4 columns represanting the membership values or a 
            vector with the membership values of the first class}
  \item{f}{vector with the membership values of the second class}
  \item{g}{vector with the membership values of the third class}
  \item{h}{vector with the membership values of the forth class}
  \item{angle}{angle between x and y axis }
  \item{scale.y}{scale of y axis related to x- and z axis}
  \item{label}{label for the classes }
  \item{labelcol}{colors to use for the labels}
  \item{labelpch}{\code{pch} for the labels}
  \item{labelcex}{\code{cex} for the labels}
  \item{main}{main title of the plot}
  \item{s3d.control}{a \emph{list} with further arguments passed to the underlying 
    \code{\link[scatterplot3d]{scatterplot3d}} function call that sets up the plot}
  \item{simplex.control}{a \emph{list} with further arguments passed to the underlying 
    function call that draws the barycentric coordinate system}
  \item{legend.control}{a \emph{list} with further arguments passed to the underlying 
    function call that adds the \code{\link{legend}}}
  \item{\dots}{further arguments passed to the underlying \code{plot} function that draws the data points}
}
\details{
The membership values are calculated with \code{\link{quadtrafo}} and plotted 
with \code{\link[scatterplot3d]{scatterplot3d}}.
}
\value{
A \code{\link[scatterplot3d]{scatterplot3d}} object.
}

\references{Garczarek, Ursula Maria (2002): Classification rules in standardized partition spaces.
            Dissertation, University of Dortmund. 
            URL \url{http://hdl.handle.net/2003/2789}}

\author{Karsten Luebke, \email{karsten.luebke@fom.de}, and Uwe Ligges}

\seealso{\code{\link{triplot}}, \code{\link[scatterplot3d]{scatterplot3d}}}
\examples{
library("MASS")
data(B3)
opar <- par(mfrow = c(1, 2), pty = "s")
posterior <- predict(lda(PHASEN ~ ., data = B3))$post
s3d <- quadplot(posterior, col = rainbow(4)[B3$PHASEN], 
        labelpch = 22:25, labelcex = 0.8,
        pch = (22:25)[apply(posterior, 1, which.max)], 
        main = "LDA posterior assignments")
quadlines(centerlines(4), sp = s3d, lty = "dashed")

posterior <- predict(qda(PHASEN ~ ., data = B3))$post
s3d <- quadplot(posterior, col = rainbow(4)[B3$PHASEN], 
        labelpch = 22:25, labelcex = 0.8,
        pch = (22:25)[apply(posterior, 1, which.max)],
        main = "QDA posterior assignments")
quadlines(centerlines(4), sp = s3d, lty = "dashed")
par(opar)
}
\keyword{classif}
\keyword{dplot}
\concept{Visualizing Classification Performance Measures}
\concept{Barycentric plots}
