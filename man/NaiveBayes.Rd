\name{NaiveBayes}
\alias{NaiveBayes}
\alias{NaiveBayes.default}
\alias{NaiveBayes.formula}
\title{Naive Bayes Classifier}
\description{
  Computes the conditional a-posterior probabilities of a categorical
  class variable given independent predictor variables using
  the Bayes rule.
}
\usage{
\method{NaiveBayes}{formula}(formula, data, ..., subset, na.action = na.pass)
\method{NaiveBayes}{default}(x, grouping, prior, usekernel = FALSE, fL = 0, ...)
}
\arguments{
  \item{x}{a numeric matrix, or a data frame of categorical and/or
    numeric variables.}
  \item{grouping}{class vector (a factor).}
  \item{formula}{a formula of the form \code{class ~ x1 + x2 +
      \dots}. Interactions are not allowed.}
  \item{data}{a data frame of predictors (caegorical and/or
    numeric).}
  \item{prior}{the prior probabilities of class membership. If unspecified, 
    the class proportions for the training set are used. If present, 
    the probabilities should be specified in the order of the factor levels.}
  \item{usekernel}{if \code{TRUE} a kernel density estimate (\code{\link{density}})
    is used for denstity estimation. If \code{FALSE} a normal density is estimated.}
  \item{fL}{Factor for Laplace correction, default factor is 0, i.e. no correction.}
  \item{\dots}{arguments passed to \code{\link{density}}.}
  \item{subset}{for data given in a data frame, an index vector
    specifying the cases to be used in the
    training sample.  (NOTE: If given, this argument must be
    named.)}
  \item{na.action}{a function to specify the action to be taken if \code{NA}s are
          found. The default action is not to count them for the
          computation of the probability factors. An
          alternative is na.omit, which leads to rejection of cases
          with missing values on any required variable.  (NOTE: If
          given, this argument must be named.)}
}
\value{
  An object of class \code{"NaiveBayes"} including components:

  \item{apriori}{Class distribution for the dependent variable.}
  \item{tables}{A list of tables, one for each predictor variable. For each
    categorical variable a table giving, for each attribute level, the conditional
    probabilities given the target class. For each numeric variable, a
    table giving, for each target class, mean and standard deviation of
    the (sub-)variable or a object of \code{class} \code{\link{density}}.}
}
\details{
 This implementation of Naive Bayes as well as this help is based on the code by 
  David Meyer in the package e1071 but extended for kernel estimated densities and user 
  specified \code{prior} probabilities. 
  The standard naive Bayes classifier (at least this implementation)
  assumes independence of the predictor
  variables.
  }
\author{Karsten Luebke, \email{karsten.luebke@fom.de}}
\seealso{\code{\link{predict.NaiveBayes}},\code{\link{plot.NaiveBayes}},\code{\link[e1071]{naiveBayes}},\code{\link[MASS]{qda}}}
\examples{
data(iris)
m <- NaiveBayes(Species ~ ., data = iris)
}
\keyword{classif}
\keyword{category}
\concept{Naive Bayes Classification}
\concept{Kernel estimated densities}
