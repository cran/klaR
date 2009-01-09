\name{predict.NaiveBayes}
\alias{predict.NaiveBayes}
\title{Naive Bayes Classifier}
\description{
  Computes the conditional a-posterior probabilities of a categorical
  class variable given independent predictor variables using
  the Bayes rule.
}
\usage{
\method{predict}{NaiveBayes}(object, newdata, threshold = 0.001, ...)
}
\arguments{
  \item{object}{An object of class \code{"naiveBayes"}.}
  \item{newdata}{A dataframe with new predictors.}
  \item{threshold}{Value replacing cells with 0 probabilities.}
  \item{\dots}{passed to \code{\link{dkernel}} function if neccessary.}
}
\value{
A list with the conditional a-posterior
probabilities for each class and the estimated class are returned.
}
\details{
  This implementation of Naive Bayes as well as this help is based on the code by 
  David Meyer in the package e1071 but extended for kernel estimated densities.
  The standard naive Bayes classifier (at least this implementation)
  assumes independence of the predictor
  variables. For attributes with missing values, the
  corresponding table entries are omitted for prediction.
}
\author{Karsten Luebke}
\seealso{\code{\link{NaiveBayes}},\code{\link{dkernel}}\code{\link[e1071]{naiveBayes}},\code{\link[MASS]{qda}}}
\examples{

data(iris)
m <- NaiveBayes(Species ~ ., data = iris)
predict(m)
}
\keyword{classif}
\keyword{category}
\concept{Naive Bayes Classification}
\concept{Kernel estimated densities}
