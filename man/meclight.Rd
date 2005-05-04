\name{meclight.default}
\alias{meclight}
\alias{meclight.default}
\alias{meclight.formula}
\alias{meclight.matrix}
\alias{meclight.data.frame}
\alias{print.meclight}
\title{Minimal Error Classification}
\description{Computer intensive method for linear dimension reduction that minimizes the classification error directly.}
\usage{
meclight(x, ...)

\method{meclight}{default}(x, grouping, r = 1, fold = 10, ...)
\method{meclight}{formula}(formula, data = NULL, ..., subset, na.action = na.fail)
\method{meclight}{data.frame}(x, ...)
\method{meclight}{matrix}(x, grouping, ..., subset, na.action = na.fail)
}
\arguments{
  \item{x}{(required if no formula is given as the principal argument.) A matrix or data frame 
  containing the explanatory variables.}
  \item{grouping}{(required if no formula principal argument is given.) A factor specifying the class for each observation.}  
  \item{r}{Dimension of projected subspace.}
  \item{fold}{Number of Bootstrap samples.}
  \item{formula}{A formula of the form \code{groups ~ x1 + x2 + ...}. That is, the response is the grouping factor and 
  the right hand side specifies the (non-factor) discriminators. }
  \item{data}{Data frame from which variables specified in formula are preferentially to be taken.}
  \item{subset}{An index vector specifying the cases to be used in the training sample. 
  (NOTE: If given, this argument must be named.)}
  \item{na.action}{A function to specify the action to be taken if NAs are found. 
      The default action is for the procedure to fail. 
      An alternative is \code{na.omit}, 
      which leads to rejection of cases with missing values on any required variable. 
      (NOTE: If given, this argument must be named.) }
  \item{\dots}{Further arguments passed to \code{\link[MASS]{lda}}.}
}
\details{Computer intensive method for linear dimension reduction that minimizes the classification error in the projected
subspace directly. Classification is done by \code{\link[MASS]{lda}}. In contrast to the reference function minimization is
done by Nelder-Mead in \code{\link{optim}}.}
\value{
  \item{method.model}{An object of class \sQuote{lda}.}
  \item{Proj.matrix}{Projection matrix.}
  \item{B.error}{Estimated bootstrap error rate.}
  \item{B.impro}{Improvement in \code{\link[MASS]{lda}} error rate.}
  }
\references{Roehl, M.C., Weihs, C., and Theis, W. (2002): 
Direct Minimization in Multivariate Classification. \emph{Computational Statistics}, 17, 29-46.}
\author{Maria Eveslage, Karsten Luebke}
\seealso{\code{\link{predict.meclight}}}
\examples{
data(iris)
meclight.obj <- meclight(Species ~ ., data = iris)
meclight.obj
}
\keyword{classif}
\concept{Linear Discriminant Analysis}
\concept{Linear Dimension Reduction}
\concept{Classification}
