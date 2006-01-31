\name{loclda}
\alias{loclda}
\encoding{latin1}
\alias{loclda.default}
\alias{loclda.formula}
\alias{loclda.matrix}
\alias{loclda.data.frame}
\alias{print.loclda}
\title{Localized Linear Discriminant Analysis (LocLDA)}
\description{
A localized version of Linear Discriminant Analysis.
}
\usage{
loclda(x, ...)

\method{loclda}{formula}(formula, data, ..., subset, na.action)

\method{loclda}{default}(x, grouping, weight.func = function(x) 1/exp(x), 
    k = nrow(x), weighted.apriori = TRUE, ...)

\method{loclda}{data.frame} (x, ...)

\method{loclda}{matrix}(x, grouping, ..., subset, na.action)
}
\arguments{
\item{formula}{Formula of the form \sQuote{\code{groups ~ x1 + x2 + ...}}.}
\item{data}{Data frame from which variables specified in \code{formula} are to be taken.}
\item{x}{Matrix or data frame containing the explanatory variables
    (required, if \code{formula} is not given).}
\item{grouping}{(required if no \code{formula} principal argument is given.) 
    A factor specifying the class for each observation.}
\item{weight.func}{Function used to compute local weights. 
    Must be finite over the interval [0,1]. See Details below.}
\item{k}{Number of nearest neighbours used to construct localized classification rules. 
    See Details below.}
\item{weighted.apriori}{Logical: if \code{TRUE}, class prior probabilities are computed 
    using local weights (see Details below). If \code{FALSE}, equal priors for all classes 
        actually occurring in the train data are used.}
\item{subset}{An index vector specifying the cases to be used in the training sample.}
\item{na.action}{A function to specify the action to be taken if \code{NA}s are found. 
    The default action is for the procedure to fail. An alternative is \code{na.omit} 
        which leads to rejection of cases with missing values on any required variable.}
\item{\dots}{Further arguments to be passed to \code{loclda.default}.}
}
\details{
This is an approach to apply the concept of localization described by Tutz and Binder (2005) 
    to Linear Discriminant Analysis. The function \code{loclda} generates an object of class \code{loclda} 
    (see Value below). As localization makes it necessary to build an 
    individual decision rule for each test observation, 
    this rule construction has to be handled by \code{\link[klaR]{predict.loclda}}. 
    For convenience, the rule building procedure is still described here.

To classify a test observation \eqn{x_s}{x_s}, only the \code{k} nearest neighbours of 
    \eqn{x_s}{x_s} within the train data are used. Each of these k train observations 
        \eqn{x_i, i = 1,\dots,k}{x_i, i=1,...,k}, is assigned a weight \eqn{w_i}{w_i} according to 
        \deqn{w_i = K\left(\frac{||x_i-x_s||}{d_k}\right), i=1,\dots,k}{w_i := K ( ||x_i - x_s|| / d_k ), i=1,...,k,} 
        where K is the weighting function given by \code{weight.func}, \eqn{||x_i-x_s||}{||x_i - x_s||} 
        is the euclidian distance of \eqn{x_i}{x_i} and \eqn{x_s}{x_s} 
        and \eqn{d_k}{d_k} is the euclidian distance of \eqn{x_s}{x_s} 
        to its \eqn{k}{k}-th nearest neighbour. 
        With these weights for each class \eqn{A_g, g=1,\dots,G}{A_g, g=1,...,G}, 
        its weighted empirical mean \eqn{\hat{\mu}_g}{mu_g_hat} and weighted empirical 
        covariance matrix are computed. The estimated pooled (weighted) covariance matrix 
        \eqn{\hat{\Sigma}}{Sigma_hat} is then calculated from the individual weighted 
        empirical class covariance matrices. If \code{weighted.apriori} is \code{TRUE} (the default), 
        prior class probabilities are estimated according to: 
        \deqn{prior_g := \frac{\sum_{i=1}^k \left(w_i \cdot I (x_i \in A_g)\right)}{\sum_{i=1}^k \left( w_i \right)}}{prior_g := [  Sum_{i=1,..,k} ( w_i * I(x_i in A_g) )  ]   /   [  Sum_{i=1,...,k} ( w_i )  ], g = 1,...,G,} 
        where I is the indicator function. If \code{FALSE}, equal priors for all classes are used. 
        In analogy to Linear Discriminant Analysis, the decision rule for \eqn{x_s}{x_s} is 
        \deqn{\hat{A} := argmax_{g \in 1,\dots,G} (posterior_g)}{A_hat := argmax_{g in 1,...,G} (posterior_g),} 
        where \deqn{posterior_g := prior_g \cdot \exp{\left( (-\frac{1}{2}) t(x_s-\hat{\mu}_g)\hat{\Sigma}^{-1}(x_s-\hat{\mu}_g)\right)} }{posterior_g := prior_g * exp [ (-1/2) * t( x_s - mu_g_hat ) * Sigma_hat^(-1) * ( x_s - mu_g_hat ) ] .}
If \eqn{posterior_g < 10^{(-150)} \forall g \in \{1,\dots,G\}}{posterior_g < 1e-150 for all g in 1,...,G}, 
    \eqn{posterior_g}{posterior_g} is set to \eqn{\frac{1}{G}}{1/G} for all \eqn{g \in 1,\dots,G}{g in 1,...,G} 
    and the test observation \eqn{x_s}{x_s} is simply assigned to the class whose weighted mean has the lowest 
    euclidian distance to \eqn{x_s}{x_s}.
}
\value{
A list of class \code{loclda} containing the following components:
\item{call}{The (matched) function call.}
\item{learn}{Matrix containing the values of the explanatory variables for all train observations.}
\item{grouping}{Factor specifying the class for each train observation.}
\item{weight.func}{Value of the argument \code{weight.func}.}
\item{k}{Value of the argument \code{k}.}
\item{weighted.apriori}{Value of the argument \code{weighted.apriori}.}
}
\references{
  Tutz, G. and Binder, H. (2005): Localized classification. \emph{Statistics and Computing} 15, 155-166.
}
\author{Marc Zentgraf (\email{marc-zentgraf@gmx.de}) and Karsten Luebke (\email{luebke@statistik.uni-dortmund.de})}
\seealso{
  \code{\link{predict.loclda}},
  \code{\link[MASS]{lda}}
}
\examples{
benchB3("lda")$l1co.error
benchB3("loclda")$l1co.error
}
\keyword{multivariate}
\concept{Localized Linear Discriminant Analysis}
\concept{Linear Discriminant Analysis}
\concept{Localization}
