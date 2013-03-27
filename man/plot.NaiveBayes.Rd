\name{plot.NaiveBayes}
\alias{plot.NaiveBayes}
\title{Naive Bayes Plot}
\description{
  Visualizes the marginal probabilities of predictor variables 
  given the class.
}
\usage{
\method{plot}{NaiveBayes}(x, vars, n = 1000, legendplot = TRUE, lty, col,
    ylab = "Density", main = "Naive Bayes Plot", ...)
}

\arguments{
  \item{x}{an object of class \code{\link{NaiveBayes}}}
  \item{vars}{variables to be plotted. If missing, all predictor variables are plotted.}
  \item{n}{number of points used to plot the density line.}
  \item{legendplot}{logical, whether to print a \code{\link{legend}}}
  \item{lty}{line type for different classes, defaults to the first \code{length(x$apriori)} 
    colors of the current palette in use.}
  \item{col}{color for different classes, defaults to \code{rainbow(length(x$apriori))}.}
  \item{ylab}{label for y-axis.}
  \item{main}{title of the plots.}
  \item{\dots}{furhter arguments passed to the underlying plot functions.}
}
\details{
 For metric variables the estimated density is plotted. For categorial variables \code{\link{mosaicplot} is called.}
}

\author{Karsten Luebke, \email{karsten.luebke@fom.de}}

\seealso{\code{\link{NaiveBayes}}}
\examples{
data(iris)
mN <- NaiveBayes(Species ~ ., data = iris)
plot(mN)

mK <- NaiveBayes(Species ~ ., data = iris, usekernel = TRUE)
plot(mK)
}
\keyword{classif}
\keyword{dplot}
\concept{Visualizing Naive Bayes Classification}
