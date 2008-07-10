\name{shardsplot}
\alias{shardsplot}
\alias{plot.EDAM}

\title{Plotting Eight Direction Arranged Maps or Self-Organizing Maps}
\description{
  Plotting method for objects of \code{class} \code{\link{EDAM}} or \code{\link[som]{som}}. 
}
\usage{
shardsplot(object, plot.type = c("eight", "four", "points", "n"), 
    expand = 1, stck = TRUE, grd = FALSE, standardize = FALSE, 
    data.or = NA, label = FALSE, plot = TRUE, classes = 0, 
    vertices = TRUE, classcolors = "rainbow", wghts = 0, 
    xlab = "Dimension 1", ylab = "Dimension 2", xaxs = "i", 
    yaxs = "i", ...)

\method{plot}{EDAM}(...)
}

\arguments{
  \item{object}{an object of class \code{EDAM} or \code{som}.}
  \item{plot.type}{a character giving the shape of the shards. 
    Available are \dQuote{\code{eight}} and \dQuote{\code{four}} for octagons resp. rectangles, 
    and \dQuote{\code{points}} for points. If \code{plot.type} is \dQuote{\code{n}}, 
    no shards are plotted at all.}
  \item{expand}{a numeric giving the relative expansion of the axes. 
    A value greater than one implies smaller shards. Varying \code{expand}
    can be sensible for visual reasons.}
  \item{stck}{logical. If \code{TRUE} the cells are varied continously corresponding to 
    the differences of direct neighbors in the origin space. 
    Within this variation the relative order of the cells is always preserved.}
  \item{grd}{logical. If \code{TRUE} (which automatically sets \code{stck} to \code{TRUE}), 
    the variation of cells is restricted to their original discrete values.}
  \item{standardize}{logical. If \code{TRUE}, then the measurements in \code{object$preimages} 
    are standardized before calculating Euclidean distances. 
        Measurements are standardized for each variable by dividing by the variable's 
        standard deviation. Meaningless if \code{object$preimages} is a dissimilarity matrix. }
  \item{data.or}{original data and classes where the first k columns are variables and the (k+1)-th column are the classes.
    If defined and class of \code{object} is \code{som}, \code{data.or} is used to assign a class to each codebook. There
    a codebook receives the class, from which the majority of its assigned objects origins.}
  \item{label}{logical. If \code{TRUE}, the shards are labeled by the rownames of the preimages.}
  \item{plot}{logical. If \code{FALSE}, all graphical output is suppressed.}
  \item{classes}{a vector giving alternative classes for objects of class \code{EDAM}; \code{classes} have to be given in 
    the original order of the data to which \code{\link{EDAM}} was applied.}
  \item{vertices}{logical. If \code{TRUE} the grid is drawn.}
  \item{classcolors}{colors to represent the classes, or a character giving the \emph{colorscale} for the classes. 
    Since now available scales are \code{rainbow}, \code{topo} and \code{gray}.}
  \item{wghts}{ an optional vector of length k giving relative weights of the variables 
    in computing Euclidean distances. Meaningless if \code{object$preimages} is a dissimilarity matrix.}
  \item{xaxs}{see \code{\link{par}}}
  \item{yaxs}{see \code{\link{par}}}
  \item{xlab}{see \code{\link{par}}}
  \item{ylab}{see \code{\link{par}}}
  \item{\dots}{further plotting parameters.}
  }
  
\details{
  If \code{plot.type} is \dQuote{\code{four}} or \dQuote{\code{eight}}, the shape of each shard depends 
  on the relative distances of the actual object 
  or codebook to its up to eight neighbours. If \code{plot.type} is \dQuote{\code{eight}}, \code{shardsplot}
  corresponds to the representation method
  suggested by Cottrell and de Bodt (1996) for Kohonen Self-Organizing Maps. 
  If \code{plot.type} is \dQuote{\code{points}}, \code{shardsplot} reduces to a usual scatter plot.
}
\value{
  The following list is (invisibly) returned:
  \item{Cells.ex}{the images of the visualized data}
  \item{S}{the criterion of the visualization}
}
\references{Cottrell, M., and de Bodt, E. (1996).
  A Kohonen Map Representation to Avoid Misleading Interpretations.
  \emph{Proceedings of the European Symposium on Atrificial Neural Networks}, D-Facto, pp. 103--110.
}
\author{Nils Raabe}

\seealso{
    \code{\link{EDAM}}, \code{\link{TopoS}}, \code{\link[som]{som}}
}
\examples{
# Compute clusters and an Eight Directions Arranged Map for the 
# country data. Plotting the result.
data(countries)
logcount <- log(countries[,2:7])
sdlogcount <- apply(logcount, 2, sd)
logstand <- t((t(logcount) / sdlogcount) * c(1,2,6,5,5,3))
cclasses <- cutree(hclust(dist(logstand)), k = 6)
countryEDAM <- EDAM(logstand, classes = cclasses, sa = FALSE, 
    iter.max = 10, random = FALSE)
plot(countryEDAM, vertices = FALSE, label = TRUE, stck = FALSE)

# Compute and plot a Self-Organizing Map for the iris data
data(iris)
library(som)
irissom <- som(iris[,1:4], xdim = 6, ydim = 14)
shardsplot(irissom, data.or = iris, vertices = FALSE)
opar <- par(xpd = NA)
legend(7.5, 6.1, col = rainbow(3), xjust = 0.5, yjust = 0,
    legend = levels(iris[, 5]), pch = 16, horiz = TRUE)
par(opar)    
}
\keyword{hplot}
\concept{Vizualizing Eight Direction Arranged Maps}
\concept{EDAM}
\concept{SOM}
\concept{Cluster analysis}
