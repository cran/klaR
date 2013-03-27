\name{e.scal}
\alias{e.scal}
\title{Function to calculate e- or softmax scaled membership values}
\description{
Calculates the e- or softmax scaled membership values of an argmax based classification rule.
}
\usage{
e.scal(x, k = 1, tc = NULL)
}
\arguments{
  \item{x}{matrix of membership values}
  \item{k}{parameter for e-scaling (1 for softmax)}
  \item{tc}{vector of true classes (required if \code{k} has to be optimized)}
}
\details{
For any  membership vector y \eqn{\exp(y\cdot k) / \sum\exp(y\cdot k)}{exp(y*k) / sum(exp(y*k)} is calculated. 
If \code{k=1}, the classical softmax scaling is used. If the true classes are given, \code{k} is optimized 
so that the apparent error rate is minimized.
}
\value{
A list containing elements
  \item{sv}{Scaled values}
  \item{k}{Optimal \code{k}}
}
\references{Garczarek, Ursula Maria (2002): Classification rules in standardized partition spaces.
            Dissertation, University of Dortmund. 
            URL \url{http://hdl.handle.net/2003/2789}}

\author{Karsten Luebke, \email{karsten.luebke@fom.de}}

\examples{
library(MASS)
data(iris)
ldaobj <- lda(Species ~ ., data = iris)
ldapred <- predict(ldaobj)$posterior
e.scal(ldapred)
e.scal(ldapred, tc = iris$Species)
}
\keyword{classif}
\concept{Softmax scaling}
\concept{E scaling}
\concept{Classification}
\concept{Membership values}
