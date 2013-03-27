\name{sknn}
\alias{sknn}
\alias{sknn.default}
\alias{sknn.formula}
\alias{sknn.matrix}
\alias{sknn.data.frame}
\title{Simple k nearest Neighbours}
\description{
Function for simple knn classification.
}
\usage{
sknn(x, ...)

\method{sknn}{default}(x, grouping, kn = 3, gamma=0, ...)
\method{sknn}{data.frame}(x, ...)
\method{sknn}{matrix}(x, grouping, ..., subset, na.action = na.fail)
\method{sknn}{formula}(formula, data = NULL, ..., subset, na.action = na.fail)
}

\arguments{
  \item{x}{matrix or data frame containing the explanatory variables 
    (required, if \code{formula} is not given).}
  \item{grouping}{factor specifying the class for each observation 
    (required, if \code{formula} is not given).}
  \item{formula}{formula of the form \code{groups ~ x1 + x2 + ...}. 
    That is, the response is the grouping factor and the right hand side specifies the (non-factor) discriminators.}
  \item{data}{Data frame from which variables specified in \code{formula} are preferentially to be taken.}
  \item{kn}{Number of nearest neighbours to use.}
  \item{gamma}{gamma parameter for rbf in knn. If \code{gamma=0} ordinary knn classification is used.}
  \item{subset}{An index vector specifying the cases to be used in the
    training sample. (Note: If given, this argument must be named.)}
  \item{na.action}{specify the action to be taken if \code{NA}s are
    found. The default action is for the procedure to fail. An
    alternative is \code{\link{na.omit}}, which leads to rejection of cases with
    missing values on any required variable. (Note: If given, this
    argument must be named.) }
  \item{...}{currently unused}
}

\details{
If \code{gamma>0} an gaussian like density is used to weight the classes of the \code{kn} nearest neighbors. 
\code{weight=exp(-gamma*distance)}. This is similar to an rbf kernel. 
If the distances are large it may be useful to \code{\link{scale}} the data first.
}
\value{
  A list containing the function call. 
}
\author{Karsten Luebke, \email{karsten.luebke@fom.de}}

\seealso{\code{\link{predict.sknn}}, \code{\link[class]{knn}}}
\examples{
data(iris)
x <- sknn(Species ~ ., data = iris)
x <- sknn(Species ~ ., gamma = 4, data = iris)
}
\keyword{classif}
\concept{k nearest neighbour classification}
\concept{KNN}
