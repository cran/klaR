\name{classscatter}
\alias{classscatter}

\title{Classification scatterplot matrix}
\description{
Function to plot a scatterplot matrix with a classification result.
}
\usage{
classscatter(formula, data, method, col.correct = "black", 
    col.wrong = "red", gs = NULL, ...)
}

\arguments{
   \item{formula}{ formula of the form \code{groups ~ x1 + x2 + ...}. 
    That is, the response is the grouping factor and the right hand side specifies the (non-factor) discriminators.  }
  \item{data}{ Data frame from which variables specified in formula are preferentially to be taken.}
  \item{method}{character, name of classification function 
    (e.g. \dQuote{\code{\link[MASS]{lda}}}).}
  \item{col.correct}{color to use for correct classified objects.}
  \item{col.wrong}{color to use for missclassified objects.}
  \item{gs}{group symbol (plot character), must have the same length as the data. 
    If \code{NULL}, \code{as.character(groups)} is the default.}
  \item{\dots}{further arguments passed to the underlying classification method or plot functions.}
}

\value{
 The actual error rate.
}

\author{Karsten Luebke}

\seealso{\code{\link{plot}}}
\examples{
data(B3)
library(MASS)
classscatter(PHASEN ~ BSP91JW + EWAJW + LSTKJW, 
    data = B3, method = "lda")
}
\keyword{dplot}
\keyword{classif}
\concept{Visualizing classification results}
