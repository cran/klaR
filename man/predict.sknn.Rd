\name{predict.sknn}
\alias{predict.sknn}
\title{Simple k Nearest Neighbours Classification}
\description{
Classifies new observations using the sknn learned by 
the \code{sknn}-function.
}
\usage{
\method{predict}{sknn}(object, newdata,...)
}
\arguments{
\item{object}{Object of class \code{\link{sknn}}.}
\item{newdata}{Data frame (or matrix) of cases to be classified.}
\item{...}{\dots}
}
\value{
A list with elements \sQuote{\code{class}} and \sQuote{\code{posterior}}.
}
\author{Karsten Luebke, \email{karsten.luebke@fom.de}}
\seealso{
  \code{\link{sknn}}, \code{\link[class]{knn}}
}
\examples{
data(iris)
x <- sknn(Species ~ ., data = iris)
predict(x, iris)
x <- sknn(Species ~ ., gamma = 10, kn = 10, data = iris)
predict(x, iris)
}
\keyword{classif}
\concept{k nearest neighbour classification}
\concept{KNN}
