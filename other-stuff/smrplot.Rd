\name{smrplot}
\alias{smrplot}
\title{Plotting of membership representation simplex}
\description{
For a 3 or 4 class discrimination problem the membership values of each class are visualized in
a barycentric coordinate system. 
}
\usage{
smrplot(posterior, trueclass = NULL, center = FALSE, 
    col = rainbow(ncol(posterior)), pchest = rep(21, ncol(posterior)), 
    pchwrong = rep(20, ncol(posterior)), 
    pchclass = rep(19, ncol(posterior)), ...)
}
\arguments{
  \item{posterior}{(scaled) posteriori probabilities of different classes.}
  \item{trueclass}{vector of true classes (if known).}
  \item{center}{should the partitions to center be plotted?}
  \item{col}{colors to use for different classes.}
  \item{pchest}{pch symbol to use for estimated membership values.}
  \item{pchwrong}{pch symbol to use for misclassifications.}
  \item{pchclass}{pch symbol to mark the classes.}
  \item{...}{further graphical parameters to be passed to the underlying plot functions}
}

\references{Garczarek, Ursula Maria (2002): Classification rules in standardized partition spaces.
            Dissertation, University of Dortmund. 
            URL \url{http://eldorado.uni-dortmund.de:8080/FB5/ls7/forschung/2002/Garczarek}}
            
\author{Karsten Luebke, \email{luebke@statistik.uni-dortmund.de}}
\examples{
library(MASS)
data(iris)
smrplot(predict(lda(Species ~ ., data = iris))$posterior, 
    trueclass = iris$Species, center = TRUE)
}
\keyword{classif}
\keyword{dplot}
