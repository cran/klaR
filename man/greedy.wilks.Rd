\name{greedy.wilks}
\alias{greedy.wilks}
\alias{greedy.wilks.default}
\alias{greedy.wilks.formula}
\alias{print.greedy.wilks}
\title{Stepwise forward variable selection for classification}
\description{
 Performs a stepwise forward variable/model selection using the Wilk's Lambda criterion.
}
\usage{
greedy.wilks(X, ...)
\method{greedy.wilks}{default}(X, grouping, niveau = 0.2, ...)
\method{greedy.wilks}{formula}(formula, data = NULL, ...) 
}

\arguments{
  \item{X}{matrix or data frame (rows=cases, columns=variables) }
  \item{grouping}{class indicator vector}
  \item{formula}{formula of the form \sQuote{\code{groups ~ x1 + x2 + ...}}}
  \item{data}{data frame (or matrix) containing the explanatory variables}
  \item{niveau}{level for the approximate F-test decision}
  \item{...}{further arguments to be passed to the default method, e.g. \code{niveau}}
}
\details{
 A stepwise forward variable selection is performed. The initial model is defined by starting with the
 variable which separates the groups most. The model is then extended by including further variables depending 
 on the Wilk's lambda criterion: Select the one which minimizes the Wilk's lambda of the model including the
 variable if its p-value still shows statistical significance.
}
\value{
  A list of two components, a \code{formula} of the form \sQuote{\code{response ~ list + of + selected + variables}},
  and a data.frame \code{results} containing the following variables:
      \item{vars }{the names of the variables in the final model in the order of selection.}
      \item{Wilks.lambda }{the appropriate Wilks' lambda for the selected variables.}
      \item{F.statistics.overall}{the approximated F-statistic for the so far selected model.}
      \item{p.value.overall}{the appropriate p-value of the F-statistic.}
      \item{F.statistics.diff}{the approximated F-statistic of the partial Wilks's lambda (for comparing the model 
                               including the new variable with the model not including it).}
      \item{p.value.diff}{the appropriate p-value of the F-statistic of the partial Wilk's lambda.}
}
\references{Mardia, K. V. , Kent, J. T. and Bibby, J. M. (1979), \emph{Multivariate analysis}, 
            Academic Press (New York; London) }
\author{ Andrea Preusser, Karsten Luebke }
\seealso{ \code{\link{stepclass}}, \code{\link[stats]{manova}} }
\examples{
data(B3)
gw_obj <- greedy.wilks(PHASEN ~ ., data = B3, niveau = 0.1)
gw_obj
## now you can say stuff like
## lda(gw_obj$formula, data = B3)
}
\keyword{multivariate}
\concept{Stepwise variable selection in classification}
\concept{Wilk's lambda}
