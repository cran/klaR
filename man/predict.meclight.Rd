\name{predict.meclight}
\alias{predict.meclight}
\title{Prediction of Minimal Error Classification}
\description{
  Classify multivariate observations in conjunction with \code{\link{meclight}} and
  \code{\link[MASS]{lda}}. 
}
\usage{
\method{predict}{meclight}(object, newdata,...)
}
\arguments{
  \item{object}{Object of class \code{\link{meclight}}.}
  \item{newdata}{Data frame of cases to be classified or, 
    if object has a formula, a data frame with columns of the same names
    as the variables used. A vector will be interpreted as a row
    vector.}
  \item{...}{currently ignored}
}
\details{
Classify multivariate observations in conjunction with \code{\link{meclight}} and
\code{\link[MASS]{lda}}. 
}
\value{
  \item{class}{The estimated class (\code{factor}).}
  \item{posterior}{Posterior probabilities for the classes.}
}
\references{Roehl, M.C., Weihs, C., and Theis, W. (2002): 
Direct Minimization in Multivariate Classification. \emph{Computational Statistics}, 17, 29-46.}
\author{Karsten Luebke, \email{karsten.luebke@fom.de}}
\seealso{\code{\link{meclight}}}
\examples{
data(iris)
meclight.obj <- meclight(Species ~ ., data = iris)
predict(meclight.obj, iris)
}
\keyword{classif}
\concept{Linear Discriminant Analysis}
\concept{Linear Dimension Reduction}
\concept{Classification}
