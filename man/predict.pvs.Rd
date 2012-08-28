\name{predict.pvs}
\alias{predict.pvs}
\title{predict method for pvs objects}
\description{
Prediction of class membership and posterior probabilities using pairwise variable selection.
}
\usage{
\method{predict}{pvs}(object, newdata, quick = FALSE, detail = FALSE, ...)
}
\arguments{
\item{object}{an object of class \sQuote{\code{pvs}}, as that created by the function \dQuote{\code{\link[klaR]{pvs}}} }
\item{newdata}{a data frame or matrix containing new data. If not given the same datas as used for training the \sQuote{\code{pvs}}-model are used. }
\item{quick}{indicator (logical), whether a quick, but less accurate computation of posterior probabalities should be used or not.}
\item{detail}{indicator (logical), whether the returned object includes additional information about the posterior probabilities for each date in each submodel.}
\item{\dots}{Further arguments are passed to underlying \code{predict} calls.}
}
\details{
If \dQuote{\code{quick=FALSE}} the posterior probabilites for each case are computed using the pairwise coupling algorithm presented by Hastie, Tibshirani (1998). 
If \dQuote{\code{quick=FALSE}} a much quicker solution is used, which leads to less accurate posterior probabalities. 
In almost all cases it doesn't has a negative effect on the classification result.
}
\value{
a list with components:
\item{class}{the predicted classes}
\item{posterior}{posterior probabilities for the classes}
\item{details}{(only if \dQuote{\code{details=TRUE}}. A list containing matrices of posterior probabalities computated by the classification method for each case and classpair. }
}
\author{Gero Szepannek, \email{szepannek@statistik.tu-dortmund.de}, Christian Neumann}
\references{Szepannek, G. and Weihs, C. (2006)  Variable Selection for Classification of More than Two 
Classes Where the Data are Sparse. In \emph{From Data and Information Analysis to Kwnowledge Engineering.},
  eds Spiliopolou, M., Kruse, R., Borgelt, C., Nuernberger, A. and Gaul, W. pp. 700-708. Springer, Heidelberg.}
\seealso{
For more details and examples how to use this predict method, see \code{\link[klaR]{pvs}}.
}
\keyword{classif}
\concept{Pairwise variable selection for classification}
