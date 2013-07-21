\name{woe}
\alias{woe}
\alias{woe.default}
\alias{woe.formula}
\alias{print.woe}
\title{Weights of evidence}
\description{Computes weight of evidence transform of factor variables for binary classification.}
\usage{
woe(x, ...)
\method{woe}{default}(x, grouping, zeroadj = 0, ids = NULL, 
                      appont = TRUE, ...)
\method{woe}{formula}(formula, data = NULL, ...)
}
\arguments{
  \item{x}{A matrix or data frame containing the explanatory variables.}
  \item{grouping}{A factor specifying the binary class for each observation.}
  \item{formula}{A formula of the form \code{grouping ~ x1 + x2 + ...} That is, the response is the grouping factor and the right hand side specifies the discriminators.}
  \item{data}{Data frame from which variables specified in formula are to be taken.}
  \item{zeroadj}{Additive constant to be added for a level with 0 observations in a class.}  
  \item{ids}{Vector of either indices or variable names that specifies the variables to be transformed.}
  \item{appont}{Application on training data: logical indicating whether the transformed values for the training data should be returned by recursive calling of \code{predict.woe}.}
  \item{\dots}{For \code{woe.formula}: Further arguments passed to function \code{woe.default} such as \code{ids}. For \code{woe.default}: \code{replace = FALSE} can be passed to recursive call of \code{predict.woe} if \code{appont}.}
}
\details{
To each factor level \eqn{x} a numeric value \eqn{WOE(x) = ln(f(x|1)/f(x|2))} is assigned where 1 and 2 denote the class labels. The WOE transform is motivated for subsequent modelling by logistic regression. Note that the frequencies of the classes should be investigated before. Information values heuristically quantify the discriminatory power of a variable by \eqn{IV = (f(x|1)-f(x|2)) ln(f(x|1)/f(x|2))}.
}
\value{
  Returns an object of class \emph{woe} that can be applied to new data. 
  \item{woe}{WOE coefficients for factor2numeric transformation of each (specified) variable.}
  \item{IV}{Vector of information values of all transformed variables.}
  \item{newx}{Data frame of transformed data if \code{appont}.}
}

\references{
Good, I. (1950): \emph{Probability and the Weighting of Evidences.} Charles Griffin, London.

Kullback, S. (1959): \emph{Information Theory and Statistics.} Wiley, New York.
}

\author{Gero Szepannek}
\seealso{\code{\link{predict.woe}}, \code{\link{plot.woe}}}

\examples{
## load German credit data
library("evtree")
data("GermanCredit")

## training/validation split
train <- sample(nrow(GermanCredit), round(0.6*nrow(GermanCredit)))

woemodel <- woe(credit_risk~., data = GermanCredit[train,], zeroadj=0.5, applyontrain = TRUE)
woemodel

## plot variable information values and woes
plot(woemodel)
plot(woemodel, type = "woes")

## apply woes 
traindata <- predict(woemodel, GermanCredit[train,], replace = TRUE)
str(traindata)

## fit logistic regression model
glmodel     <- glm(credit_risk~., traindata, family=binomial)
summary(glmodel)
pred.trn <- predict(glmodel, traindata, type = "response")

## predict validation data
validata <- predict(woemodel, GermanCredit[-train,], replace = TRUE)
pred.val <- predict(glmodel, validata, type = "response")

}
\keyword{classif}
\keyword{multivariate}
