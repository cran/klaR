\name{predict.rda}
\alias{predict.rda}
\encoding{latin1}
\title{Regularized Discriminant Analysis (RDA)}
\description{
Classifies new observations using parameters determined by 
the \code{rda}-function.
}
\usage{
\method{predict}{rda}(object, newdata, posterior = TRUE, 
    aslist = TRUE,...)
}
\arguments{
\item{object}{Object of class \code{rda}.}
\item{newdata}{Data frame (or matrix) of cases to be classified.}
\item{posterior}{Logical; indicates whether a matrix of 
posterior probabilites over all classes for each observation 
shall be returned in addition to classifications.}
\item{aslist}{Logical; if \code{TRUE}, a list
containing classifications and posterior probabilities is returned, 
otherwise a vector with an attribute \sQuote{\code{posterior}}.}
\item{...}{}
}
\value{
Depends on the value of argument \sQuote{\code{aslist}}:

Either a vector (of class \code{factor}) of classifications 
that (optionally) has an attribute \sQuote{\code{posterior}} 
containing the posterior probability matrix, or

A list with elements \sQuote{\code{class}} and \sQuote{\code{posterior}}.
}
\author{Christian Röver, \email{roever@statistik.uni-dortmund.de}}
\seealso{
  \code{\link[klaR]{rda}}
}
\examples{
data(iris)
x <- rda(Species ~ ., data = iris, gamma = 0.05, lambda = 0.2)
predict(x, iris[, 1:4])
}
\keyword{multivariate}
\concept{Regularized Discriminant Analysis}
\concept{Linear Discriminant Analysis}
\concept{Quadratic Discriminant Analysis}
\concept{robust against multicollinearity}
