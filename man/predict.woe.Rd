\name{predict.woe}
\alias{predict.woe}
\title{Weights of evidence}
\description{Applies weight of evidence transform of factor variables for binary classification based on a model of class \code{\link{woe}}.}
\usage{
\method{predict}{woe}(object, newdata, replace = TRUE, ...)
}
\arguments{
  \item{object}{Object resulting from a call of \code{woe}.}
  \item{newdata}{A matrix or data frame where WOE transform should be applied of the same dimension as the data used for training the \code{woe} object.}
  \item{replace}{Logical flag specifying whether the original factor variables should be kept in the output.}
  \item{\dots}{Currently not used.}        
}

\value{
Data frame including the transformed numeric \code{woe} variables. 
}

\references{
Good, I. (1950): \emph{Probability and the Weighting of Evidences.} Charles Griffin, London.
}

\author{Gero Szepannek}
\seealso{\code{\link{woe}}, \code{\link{plot.woe}}}

\examples{
# see examples in ?woe
}
\keyword{classif}
\keyword{multivariate}
