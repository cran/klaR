\name{triplot}
\alias{triplot}
\encoding{latin1}
\title{Barycentric plots}
\description{
Function to produce triangular (barycentric) plots 
illustrating proportions of 3 components, 
e.g. discrete 3D-distributions or mixture fractions that sum up to 1.
}
\usage{
triplot(x = NULL, y = NULL, z = NULL, main = "", frame = TRUE, 
    label = 1:3, grid = seq(0.1, 0.9, by = 0.1), center = FALSE, 
    set.par = TRUE, ...)
}
\arguments{
  \item{x}{Vector of fractions of first component
    OR 3-column matrix containing all three components (omitting \code{y} and \code{z})
    OR 3-element vector (for all three components, omitting \code{y} and \code{z}).}
  \item{y}{(Optional) vector of fractions of second component.}
  \item{z}{(Optional) vector of fractions of third component.}
  \item{main}{Main title}
  \item{frame}{Controls whether a frame (triangle) and labels are drawn.}
  \item{label}{(Character) vector of labels for the three corners.}
  \item{grid}{Values along which grid lines are to be drawn (or \code{FALSE} for no grid at all). 
    Default is steps of 10 percent.}
  \item{center}{Controls whether or not to draw centerlines at which there is a 
    \sQuote{tie} between any two dimensions (see also \code{\link{centerlines}}).}
  \item{set.par}{Controls whether graphical parameter \code{mar} is set so 
    the plot fills the window (see \code{\link{par}}).}
  \item{\dots}{Further graphical parameters passed to \code{\link{trilines}}.}
}
\details{
The barycentric plot illustrates the set of points (x,y,z) with x,y,z between 0 and 1 and x+y+z=1; 
that is, the triangle spanned by (1,0,0), (0,1,0) and (0,0,1) in 3-dimensional space.
The three dimensions x, y and z correspond to lower left, upper and lower right corner of the plot.
The greater the share of x in the proportion, the closer the point is to the lower left corner;
Points on the opposite (upper right) side have a zero x-fraction.
The grid lines show the points at which one dimension is held constant, horizontal lines for 
example contain points with a constant second dimension.
}
\author{Christian Röver, \email{roever@statistik.uni-dortmund.de}}
\seealso{
\code{\link{tripoints}}, \code{\link{trilines}}, \code{\link{triperplines}}, \code{\link{trigrid}}, 
\code{\link{triframe}} for points, lines and layout, \code{\link{tritrafo}} for placing labels,
and \code{\link{quadplot}} for the same in 4 dimensions.
}
\examples{
# illustrating probabilities: 
triplot(label = c("1, 2 or 3", "4 or 5", "6"), 
    main = "die rolls: probabilities", pch = 17)
triperplines(1/2, 1/3, 1/6)

# expected... 
triplot(1/2, 1/3, 1/6, label = c("1, 2 or 3", "4 or 5", "6"), 
    main = "die rolls: expected and observed frequencies", pch = 17)
# ... and observed frequencies. 
dierolls <- matrix(sample(1:3, size = 50*20, prob = c(1/2, 1/3, 1/6), 
                          replace = TRUE), ncol = 50)
frequencies <- t(apply(dierolls, 1, 
    function(x)(summary(factor(x, levels = 1:3)))) / 50)
tripoints(frequencies)

# LDA classification posterior: 
data(iris)
require(MASS)
pred <- predict(lda(Species ~ ., data = iris),iris)
plotchar <- rep(1,150)
plotchar[pred$class != iris$Species] <- 19
triplot(pred$posterior, label = colnames(pred$posterior), 
        main = "LDA posterior assignments", center = TRUE, 
        pch = plotchar, col = rep(c("blue", "green3", "red"), rep(50, 3)), 
        grid = TRUE)
legend(x = -0.6, y = 0.7, col = c("blue", "green3", "red"), 
    pch = 15, legend = colnames(pred$posterior))
}
\keyword{aplot}
\concept{Visualizing Classification Performance Measures}
\concept{Barycentric plots}
