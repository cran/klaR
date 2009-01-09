\name{plineplot}
\alias{plineplot}

\title{Plotting marginal posterior class probabilities}
\description{
For a given variable the posteriori probabilities of the classes given by a 
classification method are plotted. The variable need not be used for the actual
classifcation.
}
\usage{
plineplot(formula, data, method, x, col.wrong = "red", 
          ylim = c(0, 1), loo = FALSE, mfrow, ...)
}

\arguments{
  \item{formula}{formula of the form \code{groups ~ x1 + x2 + ...}. 
    That is, the response is the grouping factor and the right hand side specifies the (non-factor) discriminators.}
  \item{data}{ Data frame from which variables specified in formula are preferentially to be taken.}
  \item{method}{character, name of classification function 
    (e.g. \dQuote{\code{\link[MASS]{lda}}}).}
  \item{x}{variable that should be plotted. See examples.}
  \item{col.wrong}{color to use for missclassified objects.}
  \item{ylim}{\code{ylim} for the plot.}
  \item{loo}{logical, whether leave-one-out estimate is used for prediction}
  \item{mfrow}{number of rows and columns in the graphics device, see \code{\link{par}}. 
    If missing, number of rows equals number of classes, and 1 column.}
  \item{\dots}{further arguments passed to the underlying classification method or plot functions.}
}

\value{
The actual error rate.
}
\author{Karsten Luebke}
\seealso{\code{\link{partimat}}}
\examples{
library(MASS)

# The name of the variable can be used for x
data(B3)
plineplot(PHASEN ~ ., data = B3, method = "lda", 
    x = "EWAJW", xlab = "EWAJW")

# The plotted variable need not be in the data
data(iris)
iris2 <- iris[ , c(1,3,5)]
plineplot(Species ~ ., data = iris2, method = "lda", 
    x = iris[ , 4], xlab = "Petal.Width")
}
\keyword{classif}
\keyword{dplot}
\concept{Vizualizing classification results}
