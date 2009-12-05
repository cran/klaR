\name{nm}
\alias{nm}
\alias{nm.default}
\alias{nm.formula}
\alias{nm.matrix}
\alias{nm.data.frame}
\title{Nearest Mean Classification}
\description{
Function for nearest mean classification.
}
\usage{
nm(x, ...)

\method{nm}{default}(x, grouping, gamma = 0, ...)
\method{nm}{data.frame}(x, ...)
\method{nm}{matrix}(x, grouping, ..., subset, na.action = na.fail)
\method{nm}{formula}(formula, data = NULL, ..., subset, na.action = na.fail)
}

\arguments{
  \item{x}{matrix or data frame containing the explanatory variables 
    (required, if \code{formula} is not given).}
  \item{grouping}{factor specifying the class for each observation 
    (required, if \code{formula} is not given).}
  \item{formula}{formula of the form \code{groups ~ x1 + x2 + ...}. 
    That is, the response is the grouping factor and the right hand side specifies the (non-factor) discriminators.}
  \item{data}{Data frame from which variables specified in \code{formula} are preferentially to be taken.}
  \item{gamma}{gamma parameter for rbf weight of the distance to mean. If \code{gamma=0} the posterior is 1 for the
    nearest class (mean) and 0 else.}
  \item{subset}{An index vector specifying the cases to be used in the
    training sample. (Note: If given, this argument must be named.)}
  \item{na.action}{specify the action to be taken if \code{NA}s are
    found. The default action is for the procedure to fail. An
    alternative is \code{\link{na.omit}}, which leads to rejection of cases with
    missing values on any required variable. (Note: If given, this
    argument must be named.) }
  \item{...}{}
}

\details{
\code{nm} is calling \code{\link{sknn}} with the class means as observations.
If \code{gamma>0} a gaussian like density is used to weight the distance to the class means
\code{weight=exp(-gamma*distance)}. This is similar to an rbf kernel. 
If the distances are large it may be useful to \code{\link{scale}} the data first.
}
\value{
  A list containing the function call and the class means (\code{learn})). 
}
\author{Karsten Luebke, \email{karsten.luebke@fom.de}}

\seealso{\code{\link{sknn}}, \code{\link{rda}}, \code{\link[class]{knn}}}
\examples{
data(B3)
x <- nm(PHASEN ~ ., data = B3)
x$learn
x <- nm(PHASEN ~ ., data = B3, gamma = 0.1)
predict(x)$post
}
\keyword{classif}
\concept{Nearest Mean Classification}
