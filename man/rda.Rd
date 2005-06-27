\name{rda}
\alias{rda}
\encoding{latin1}
\alias{rda.default}
\alias{rda.formula}
\alias{plot.rda}
\alias{print.rda}
\title{Regularized Discriminant Analysis (RDA)}
\description{
Builds a classification rule using regularized group covariance 
matrices that are supposed to be more robust against 
multicollinearity in the data.
}
\usage{
rda(x, ...)

\method{rda}{default}(x, grouping = NULL, prior = NULL, gamma = NA, 
    lambda = NA, regularization = c(gamma = gamma, lambda = lambda), 
    crossval = TRUE, fold = 10, train.fraction = 0.5, 
    estimate.error = TRUE, output = FALSE, startsimplex = NULL, 
    max.iter = 100, trafo = TRUE, simAnn = FALSE, schedule = 2, 
    T.start = 0.1, halflife = 50, zero.temp = 0.01, alpha = 2, 
    K = 100, ...)
\method{rda}{formula}(formula, data, ...)
}
\arguments{
\item{x}{Matrix or data frame containing the explanatory variables 
    (required, if \code{formula} is not given).}
\item{formula}{Formula of the form \sQuote{\code{groups ~ x1 + x2 + ...}}.}
\item{data}{A data frame (or matrix) containing the explanatory 
    variables.}
\item{grouping}{(Optional) a vector specifying the class for 
    each observation; if not specified, the first column of 
    \sQuote{\code{data}} is taken.}
\item{prior}{(Optional) prior probabilities for the classes. 
    Default: proportional to training sample sizes. 
    \dQuote{\code{prior=1}} indicates equally likely classes.}
\item{gamma, lambda, regularization}{
    One or both of the rda-parameters may be fixed manually. 
    Unspecified parameters are determined by minimizing the 
    estimated error rate (see below).}
\item{crossval}{Logical. If \code{TRUE}, in the optimization 
    step the error rate is estimated by Cross-Validation, 
    otherwise by drawing several training- and test-samples.}
\item{fold}{The number of Cross-Validation- or Bootstrap-samples
    to be drawn.}
\item{train.fraction}{In case of Bootstrapping: the fraction of 
    the data to be used for training in each Bootstrap-sample; 
    the remainder is used to estimate the misclassification rate.}
\item{estimate.error}{Logical. If \code{TRUE}, the apparent 
    error rate for the final parameter set is estimated.}
\item{output}{Logical flag to indicate whether text output 
    during computation is desired.}
\item{startsimplex}{(Optional) a starting simplex for the 
    Nelder-Mead-minimization.}
\item{max.iter}{Maximum number of iterations for Nelder-Mead.}
\item{trafo}{Logical; indicates whether minimization is carrried 
    out using transformed parameters.}
\item{simAnn}{Logical; indicates whether Simulated Annealing 
    shall be used.}
\item{schedule}{Annealing schedule 1 or 2 (exponential or polynomial).}
\item{T.start}{Starting temperature for Simulated Annealing.}
\item{halflife}{Number of iterations until temperature is reduced to a half  (schedule 1).}
\item{zero.temp}{Temperature at which it is set to zero  (schedule 1).} 
\item{alpha}{Power of temperature reduction (linear, quadratic, cubic,...)  (schedule 2).}
\item{K}{Number of iterations until temperature = 0  (schedule 2).}
\item{...}{}
}
\details{
J.H. Friedman (see references below) suggested a method to fix 
almost singular covariance matrices in discriminant analysis. 
Basically, individual covariances as in QDA are used, but 
depending on two parameters (\eqn{\gamma}{gamma} and 
\eqn{\lambda}{lambda}), these can be shifted towards a 
diagonal matrix and/or the pooled covariance 
matrix. For (\eqn{\gamma=0}{gamma=0}, \eqn{\lambda=0}{lambda=0}) it equals QDA, 
for (\eqn{\gamma=0}{gamma=0}, \eqn{\lambda=1}{lambda=1}) it equals LDA.

You may fix these parameters at certain values or leave it to 
the function to try to find \dQuote{optimal} values. If one 
parameter is given, the other one is determined using the 
R-function \sQuote{\code{\link{optimize}}}. If no parameter is 
given, both are determined numerically by a 
Nelder-Mead-(Simplex-)algorithm with the option of using 
Simulated Annealing.
The goal function to be minimized is the (estimated) 
misclassification rate; the misclassification rate is estimated 
either by Cross-Validation or by repeatedly dividing the data 
into training- and test-sets (Boostrapping).

\emph{Warning}: If these sets are small, optimization is expected 
to produce almost random results. We recommend to adjust the 
parameters manually in such a case.
In all other cases it is recommended to run the optimization 
several times in order to see whether stable results are gained.

Since the Nelder-Mead-algorithm is actually intended for 
\emph{continuous} functions while the observed error rate 
by its nature is \emph{discrete}, a greater number of 
Boostrap-samples might improve the optimization by increasing 
the smoothness of the response surface (and, of course, by 
reducing variance and bias). 
If a set of parameters leads to singular covariance 
matrices, a penalty term is added to the misclassification rate 
which will hopefully help to maneuver back out of singularity
(so do not worry about error rates greater than one during 
optimization).
}

\section{More details}{
The explicit defintion of \eqn{\gamma}{gamma},
\eqn{\lambda}{lambda} and the resulting covariance estimates
is as follows:

The pooled covariance estimate \eqn{\hat{\Sigma}}{SigmaHat} is 
given as well as the individual covariance estimates 
\eqn{\hat{\Sigma}_k}{SigmaHat_k} for each group.

First, using \eqn{\lambda}{lambda}, a convex combination of 
these two is computed:
\deqn{\hat{\Sigma}_k (\lambda) := (1-\lambda) \hat{\Sigma}_k + \lambda \hat{\Sigma}.}{%
SigmaHat_k(lambda) := (1-lambda)*SigmaHat_k + lambda*SigmaHat.}
Then, another convex combination is constructed using the 
above estimate and a (scaled) identity matrix:
\deqn{\hat{\Sigma}_k (\lambda,\gamma) = (1-\gamma)\hat{\Sigma}_k(\lambda)+
\gamma\frac{1}{d}\mathrm{tr}[\hat{\Sigma}_k(\lambda)]\mathrm{I}.}{%
SigmaHat_k(lambda,gamma) := (1-gamma) SigmaHat_k(lambda) 
+ gamma 1/d trace(SigmaHat_k(lambda)) I.}
The factor 
\eqn{\frac{1}{d}\mathrm{tr}[\hat{\Sigma}_k(\lambda)]}{1/d trace(SigmaHat_k(lambda))}
in front of the identity matrix I is the mean of the diagonal 
elements of 
\eqn{\hat{\Sigma}_k(\lambda)}{SigmaHat_k(lambda)}, so it is 
the mean variance of all \eqn{d} variables assuming the group 
covariance \eqn{\hat{\Sigma}_k(\lambda)}{SigmaHat_k(lambda)}.

For the four extremes of (\eqn{\gamma}{gamma},\eqn{\lambda}{lambda}) 
the covariance structure reduces to special cases:
\itemize{
\item (\eqn{\gamma=0}{gamma=0}, \eqn{\lambda=0}{lambda=0}): 
QDA - individual covariance for each group.
\item (\eqn{\gamma=0}{gamma=0}, \eqn{\lambda=1}{lambda=1}): 
LDA - a common covariance matrix.
\item (\eqn{\gamma=1}{gamma=1}, \eqn{\lambda=0}{lambda=0}): 
Conditional independent variables - similar to Naive Bayes,
but variable variances within group (main diagonal elements) 
are equal.
\item (\eqn{\gamma=1}{gamma=1}, \eqn{\lambda=1}{lambda=1}):
Classification using euclidean distance - as in previous case, 
but variances are the same for all groups. Objects are assigned 
to group with nearest mean.
}
}
\value{
A list of class \code{rda} containing the following 
components:
\item{call}{The (matched) function call.}
\item{regularization}{vector containing the two regularization 
parameters (gamma, lambda)}
\item{classes}{the names of the classes}
\item{prior}{the prior probabilities for the classes}
\item{error.rate}{apparent error rate (if computation 
was not suppressed), and, if any optimization took place, the
final (cross-validated or bootstrapped) error rate estimate as 
well.}
\item{means}{Group means.}
\item{covariances}{Array of group covariances.}
\item{covpooled}{Pooled covariance.}
\item{converged}{(Logical) indicator of convergence (only for 
Nelder-Mead).}
\item{iter}{Number of iterations actually performed (only for 
Nelder-Mead).}
}
\references{
  Friedman, J.H. (1989): Regularized Discriminant Analysis.
  In: \emph{Journal of the American Statistical Association} 84, 
  165-175.
 
  Press, W.H., Flannery, B.P., Teukolsky, S.A., Vetterling, W.T. (1992): 
  \emph{Numerical Recipes in C}.  Cambridge: Cambridge University Press.
}
\author{Christian Röver, \email{roever@statistik.uni-dortmund.de}}
\seealso{
  \code{\link[klaR]{predict.rda}},
  \code{\link[MASS]{lda}}, \code{\link[MASS]{qda}}
}
\examples{
data(iris)
x <- rda(Species ~ ., data = iris, gamma = 0.05, lambda = 0.2)
predict(x, iris)
}
\keyword{multivariate}
\concept{Regularized Discriminant Analysis}
\concept{Linear Discriminant Analysis}
\concept{Quadratic Discriminant Analysis}
\concept{robust against multicollinearity}
