\name{predict.svmlight}
\alias{predict.svmlight}
\title{Interface to SVMlight}
\description{
Classifies new observations using the SVM learned by 
the \code{svmlight}-function.
}
\usage{
\method{predict}{svmlight}(object, newdata, scal = TRUE, ...)
}
\arguments{
\item{object}{Object of class \code{\link[klaR]{svmlight}}.}
\item{newdata}{Data frame (or matrix) of cases to be classified.}
\item{scal}{Logical, whether to scale membership values via \code{\link{e.scal}}.}
\item{...}{\dots}
}
\value{
A list with elements \sQuote{\code{class}} and \sQuote{\code{posterior}} (scaled, if \code{scal = TRUE}).
}
\author{Karsten Luebke, \email{luebke@statistik.uni-dortmund.de}}
\seealso{
  \code{\link[klaR]{svmlight}}, \code{\link[e1071]{svm}}
}
\examples{
\dontrun{
data(iris)
x <- svmlight(Species ~ ., data = iris)
predict(x, iris)
}
}
\keyword{classif}
