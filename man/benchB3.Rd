\name{benchB3}
\alias{benchB3}
\title{Benchmarking on B3 data}
\description{
Evaluates the performance of a classification method on the \code{\link{B3}} data. 
}
\usage{
benchB3(method, prior = rep(1/4, 4), sv = "4", scale = FALSE, ...)
}
\arguments{
  \item{method}{classification method to use}
  \item{prior}{prior probabilities of classes}
  \item{sv}{class of the start of a business cycle}
  \item{scale}{logical, whether to use \code{\link{scale}} first}
  \item{\dots}{furhter arguments passed to \code{method}}
}
\details{
The performance of classification methods on cyclic data can be measured by a special form 
of cross-validation: Leave-One-Cycle-Out. That means that a complete cycle is used as test data 
and the others are used as training data. This is repeated for all complete cycles in the data. 
}
\value{
  A list with elements
  \item{MODEL}{list with the model returned by \code{method} of the training data}
  \item{error}{vector of test error rates in cycles}
  \item{l1co.error}{leave-one-cycle-out error rate}
}

\author{Karsten Luebke}

\seealso{\code{\link{B3}}}
\examples{
perLDA <- benchB3("lda")
\dontrun{
## due to parameter optimization rda takes a while 
perRDA <- benchB3("rda")
library(rpart)
## rpart will not work with prior argument:
perRpart <- benchB3("rpart", prior = NULL)
}
}
\keyword{classif}
\concept{Benchmark}
