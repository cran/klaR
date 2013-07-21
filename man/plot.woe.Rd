\name{plot.woe}
\alias{plot.woe}
\title{Plot information values}
\description{
Barplot of information values to compare dicriminator of the transformed variables.}
\usage{
\method{plot}{woe}(x, type = c("IV", "woes"), ...)
}
\arguments{
  \item{x}{An object of class \code{woe}.}
  \item{type}{Character to specify the plot type, see below. Either \code{"IV"} (default) or \code{"woes"}.}
  \item{\dots}{Further arguments to be passed to the barplot function.}
}
\details{
For \code{type=="IV"} a barplot of information values for all transformed variables. 
A thumb rule of interpretation is that Values above 0.3 are considered as strongly discrimative where values below 0.02 are considered to characterize unpredictive variables. 
For \code{type=="woes"} for each variable the relative frequencies of all transformed levels are plotted.
}
\value{No value is returned.}

\references{
Good, I. (1950): \emph{Probability and the Weighting of Evidences.} Charles Griffin, London.

Kullback, S. (1959): \emph{Information Theory and Statistics.} Wiley, New York.
}

\author{Gero Szepannek}
\seealso{\code{\link{woe}}, \code{\link{predict.woe}}}

\examples{
# see examples in ?woe
}

\keyword{classif}
\keyword{multivariate}
