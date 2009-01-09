\name{predict.locpvs}
\alias{predict.locpvs}
\title{predict method for locpvs objects}
\description{
Prediction of class membership and posterior probabilities in local models using pairwise variable selection.
}
\usage{
\method{predict}{locpvs}(object,newdata, quick = FALSE, return.subclass.prediction = TRUE, ...)
}
\arguments{
\item{object}{an object of class \sQuote{\code{locpvs}}, as that created by the function \dQuote{\code{\link[klaR]{locpvs}}}}
\item{newdata}{a data frame or matrix containing new data. If not given the same datas as used for training the \sQuote{\code{pvs}}-model are used. }
\item{quick}{indicator (logical), whether a quick, but less accurate computation of posterior probabalities should be used or not.}
\item{return.subclass.prediction}{indicator (logical), whether the returned object includes posterior probabilities for each date in each subclass }
\item{\dots}{Further arguments are passed to underlying \code{predict} calls.}
}
\details{
Posterior probabilities are predicted as if object is a standard \sQuote{\code{pvs}}-model with the subclasses as classes. Then the posterior probabalities are summed over all subclasses for each class. The class with the highest value becomes the prediction.

If \dQuote{\code{quick=FALSE}} the posterior probabilites for each case are computed using the pairwise coupling algorithm presented by Hastie, Tibshirani (1998). If \dQuote{\code{quick=FALSE}} a much quicker solution is used, which leads to less accurate posterior probabalities. In almost all cases it doesn't has a negative effect on the classification result.
}
\value{
a list with components:
\item{class}{the predicted (upper) classes}
\item{posterior}{posterior probabilities for the (upper) classes}
\item{subclass.posteriors}{(only if \dQuote{\code{return.subclass.prediction=TRUE}}. A matrix containing posterior probabalities for the subclasses. }
}
\author{Gero Szepannek, \email{szepannek@statistik.tu-dortmund.de}, Christian Neumann}
\references{Szepannek, G. and Weihs, C. (2006) Local Modelling in Classification on Different Feature Subspaces. 
In \emph{Advances in Data Mining.}, ed Perner, P., LNAI 4065, pp. 226-234. Springer, Heidelberg.}
\seealso{
\code{\link[klaR]{locpvs}} for learning \sQuote{\code{locpvs}}-models, \code{\link[klaR]{pvs}} for pairwise variable selection without modeling subclasses, \code{\link[klaR]{predict.pvs}} for predicting \sQuote{\code{pvs}}-models 
}
\examples{

## this example might be a bit artificial, but it sufficiently shows how locpvs has to be used

## learn a locpvs-model on the Vehicle dataset

library("mlbench")
data("Vehicle")

subclass <- Vehicle$Class # use four car-types in dataset as subclasses
## aggregate "bus" and "van" to upper-class "big" and "saab" and "opel" to upper-class "small":
subclass_class <- matrix(c("bus","van","saab","opel","big","big","small","small"), ncol=2) 

## learn now a locpvs-model for the subclasses:
model <- locpvs(Vehicle[,1:18], subclass, subclass_class)
model # short summary, showing the class-pairs of the submodels 
# together with the selected variables and the relation of sub- to upperclasses

## predict:
pred <- predict(model,Vehicle[,1:18])

## now you can look at the predicted classes:
pred$class
## or at the posterior probabilities:
pred$posterior
## or at the posterior probabilities for the subclasses:
pred$subclass.posteriors

}
\keyword{classif}
\concept{Pairwise variable selection for classification}
