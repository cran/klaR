\name{cond.index}
\alias{cond.index}

\title{Calculation of Condition Indices for Linear Regression}
\description{
Diagnosis of collinearity in X
}
\usage{
cond.index(formula, data, ...)
}

\arguments{

  \item{formula}{formula of the form \sQuote{\code{groups ~ x1 + x2 + ...}}}
  \item{data}{data frame (or matrix) containing the explanatory variables}
  \item{...}{further arguments to be passed to \code{lm}}
}
\details{
Collinearities can inflate the variance of the estimated regression coefficients and numerical stability. The condition indices are calculated by the eigenvalues of the crossproduct matrix of the scaled but uncentered explanatory variables.
Indices > 30 may indicate collinearity.
}
\value{
  A vector of the condition indices.
}
\references{Belsley, D. , Kuh, E. and Welsch, R. E. (1979), \emph{Regression Diagnostics: Identifying Influential Data and Sources of Collinearity}, 
            John Wiley (New York) }
\author{ Andrea Preusser, Karsten Luebke (\email{karsten.luebke@fom.de})}
\seealso{ \code{\link{stepclass}}, \code{\link[stats]{manova}} }
\examples{
data(Boston)
condition_medv <- cond.index(medv ~ ., data = Boston)
condition_medv
}
\keyword{multivariate}
\concept{Collinearity diagnosis in regression}

