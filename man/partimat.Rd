\name{partimat}
\alias{partimat}
\alias{partimat.default}
\alias{partimat.formula}
\alias{partimat.data.frame}
\alias{partimat.matrix}

\title{Plotting the 2-d partitions of classification methods}
\description{
Provides a multiple figure array which shows the classification of observations based on 
classification methods (e.g. \code{lda}, \code{qda}) for every combination of two variables. 
Moreover, the classification borders are displayed and the apparent error rates are given in each title.
}
\usage{
partimat(x,...)

\method{partimat}{default}(x, grouping, method = "lda", prec = 100, 
    nplots.vert, nplots.hor, main = "Partition Plot", name, mar, 
    plot.matrix = FALSE, plot.control = list(), ...)
\method{partimat}{data.frame}(x, ...)
\method{partimat}{matrix}(x, grouping, ..., subset, na.action = na.fail)
\method{partimat}{formula}(formula, data = NULL, ..., subset, na.action = na.fail)
}


\arguments{
  \item{x}{matrix or data frame containing the explanatory variables (required, if \code{formula} is not given).}
  \item{grouping}{factor specifying the class for each observation (required, if \code{formula} is not given).}
  \item{formula}{formula of the form \code{groups ~ x1 + x2 + ...}. 
    That is, the response is the grouping factor and the right hand side specifies the (non-factor) discriminators. }
  \item{method}{the method the classification is based on, currently supported are:
  \code{\link[MASS]{lda}}, \code{\link[MASS]{qda}}, \code{\link[rpart]{rpart}}, \code{\link[e1071]{naiveBayes}}, 
  \code{\link{rda}}, \code{\link{sknn}} and \code{\link{svmlight}}}.
  \item{prec}{precision used to draw the classification borders (the higher the more precise; default: 100).}
  \item{data}{Data frame from which variables specified in formula are preferentially to be taken.}
  \item{nplots.vert}{number of rows in the multiple figure array}
  \item{nplots.hor}{number of columns in the multiple figure array}
  \item{subset}{index vector specifying the cases to be used in the
    training sample. (Note: If given, this argument must be named.)}
  \item{na.action}{specify the action to be taken if \code{NA}s are
    found. The default action is for the procedure to fail. An
    alternative is \code{\link{na.omit}}, which leads to rejection of cases with
    missing values on any required variable. (Note: If given, this argument must be named.) }
  \item{main}{title}
  \item{name}{Variable names to be printed at the axis / into the diagonal.}
  \item{mar}{numerical vector of the form \code{c(bottom, left, top, right)}
    which gives the lines of margin to be specified on the four sides of the plot.
    Defaults are \code{rep(0, 4)} if \code{plot.matrix = TRUE}, \code{c(5, 4, 2, 1) + 0.1} otherwise.}
  \item{plot.matrix}{logical; if \code{TRUE}, like a scatterplot matrix; 
    if \code{FALSE} (default) uses less space and arranges the plots \dQuote{optimal} 
    (using a fuzzy algorithm) in an array by plotting each pair of variables once.}
  \item{plot.control}{A list containing further arguments passed to the underlying 
    plot functions (and to \code{\link{drawparti}}).}
  \item{...}{Further arguments passed to the classification \code{method} (through \code{\link{drawparti}}).}
}
\seealso{for much more fine tuning see \code{\link{drawparti}}}
\author{Karsten Luebke, \email{karsten.luebke@fom.de}, Uwe Ligges, Irina Czogiel}
\note{Warnings such as  \sQuote{parameter \dQuote{xyz} couldn't be set in high-level plot function} are expected,
    if making use of \code{...}.}
\examples{
library(MASS)
data(iris)
partimat(Species ~ ., data = iris, method = "lda")
\dontrun{
partimat(Species ~ ., data = iris, method = "lda", 
    plot.matrix = TRUE, imageplot = FALSE) # takes some time ...
}
}
\keyword{classif}
\keyword{dplot}
\concept{Vizualizing classification results}
