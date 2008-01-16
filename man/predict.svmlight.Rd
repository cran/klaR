\name{predict.svmlight}
\alias{predict.svmlight}
\title{Interface to SVMlight}
\description{
Predicts new observations using the SVM learned by 
the \code{svmlight}-function.
}
\usage{
\method{predict}{svmlight}(object, newdata, scal = TRUE, ...)
}
\arguments{
\item{object}{Object of class \code{\link{svmlight}}.}
\item{newdata}{Data frame (or matrix) of cases to be predicted.}
\item{scal}{Logical, whether to scale membership values via \code{\link{e.scal}}.}
\item{...}{\dots}
}
\value{
If a classification is learned (\code{type="C"}) in \code{\link{svmlight}} a 
list with elements \sQuote{\code{class}} and \sQuote{\code{posterior}} (scaled, if \code{scal = TRUE}).

If a Regression is learned (\code{type="R"}) in \code{\link{svmlight}} the predicted values.
}
\author{Karsten Luebke, \email{luebke@statistik.tu-dortmund.de}}
\seealso{
  \code{\link{svmlight}}, \code{\link[e1071]{svm}}
}
\examples{
\dontrun{
data(iris)
x <- svmlight(Species ~ ., data = iris)
predict(x, iris)
}
}
\keyword{classif}
\concept{Support vector machines}
\concept{SVMlight}
\concept{SVM}
\concept{Classification}
