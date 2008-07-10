\name{locpvs}
\alias{locpvs}
\title{Pairwise variable selection for classification in local models}
\description{
Performs pairwise variable selection on subclasses.
}
\usage{
locpvs(x, subclasses, subclass.labels, prior=NULL, method="lda", 
    vs.method = c("ks.test", "stepclass", "greedy.wilks"),
    niveau=0.05, fold=10, impr=0.1, direct="backward", out=FALSE, ...)
}
\arguments{
\item{x}{matrix or data frame containing the explanatory variables. x must consist of numerical data only. }
\item{subclasses}{vector indicating the subclasses (a factor) }
\item{subclass.labels}{must be a matrix with 2 coloumns, where the first coloumn specifies the subclass and the second coloumn the according upper class }
\item{prior}{prior probabilites for the classes. If not specified the prior probabilities will be set according to proportion in \dQuote{subclasses}. If specified the order of prior 
        probabilities must be the same as in \dQuote{subclasses}. }
\item{method}{character, name of classification function (e.g. \dQuote{\code{\link[MASS]{lda}}} (default)).}
\item{vs.method}{character, name of variable selection method. Must be one of \dQuote{\code{\link[ctest]{ks.test}}} (default),
    \dQuote{\code{\link[klaR]{stepclass}}} or \dQuote{\code{\link[klaR]{greedy.wilks}}}. }
\item{niveau}{used niveau for \dQuote{\code{\link[ctest]{ks.test}}}}
\item{fold}{parameter for cross-validation, if \dQuote{\code{\link[klaR]{stepclass}}} is chosen \sQuote{\code{vs.method}}}                                            
\item{impr}{least improvement of performance measure desired to include or exclude any variable (<=1), if \dQuote{\code{\link[klaR]{stepclass}}} is chosen \sQuote{\code{vs.method}} }  
\item{direct}{direction of variable selection, if \dQuote{\code{\link[klaR]{stepclass}}} is chosen \sQuote{\code{vs.method}}. 
        Must be one if \dQuote{\code{forward}}, \dQuote{\code{backward}} (default) or \dQuote{\code{both}}. }
\item{out}{indicator (logical) for textoutput during computation (slows down computation!), if \dQuote{\code{\link[klaR]{stepclass}}} is chosen \sQuote{\code{vs.method}} }
\item{...}{further parameters passed to classification function (\sQuote{\code{method}}) or variable selection method (\sQuote{\code{vs.method}}) }
}
\details{
A call on \code{\link[klaR]{pvs}} is performed using \dQuote{subclasses} as grouping variable. See \code{\link[klaR]{pvs}} for further details.
}
\value{
An object of class \sQuote{\code{locpvs}} containing the following components:
\item{pvs.result}{the complete output of the call to \code{\link[klaR]{pvs}} (see \code{\link[klaR]{pvs}} for further details}
\item{subclass.labels}{the subclass.labels as specified in function call}
}
\author{Gero Szepannek, \email{szepannek@statistik.tu-dortmund.de}, Christian Neumann}
\references{Szepannek, G. and Weihs, C. (2006) Local Modelling in Classification on Different Feature Subspaces. 
In \emph{Advances in Data Mining.}, ed Perner, P., LNAI 4065, pp. 226-234. Springer, Heidelberg.}
\seealso{
\code{\link[klaR]{predict.locpvs}} for predicting \sQuote{\code{locpvs}} models and \code{\link[klaR]{pvs}}
}
\examples{

## this example might be a bit artificial, but it sufficiently shows how locpvs has to be used

## learn a locpvs-model on the Vehicle dataset

library("mlbench")
data("Vehicle")

subclass <- Vehicle$Class # use four car-types in dataset as subclasses
## aggregate "bus" and "van" to upper-class "big" and "saab" and "opel" to upper-class "small"
subclass_class <- matrix(c("bus","van","saab","opel","big","big","small","small"),ncol=2) 

## learn now a locpvs-model for the subclasses:
model <- locpvs(Vehicle[,1:18], subclass, subclass_class) 
model # short summary, showing the class-pairs of the submodels 
# together with the selected variables and the relation of sub- to upperclasses

## predict:
pred <- predict(model, Vehicle[,1:18])

## now you can look at the predicted classes:
pred$class
## or at the posterior probabilities:
pred$posterior
## or at the posterior probabilities for the subclasses:
pred$subclass.posteriors

}
\keyword{classif}
\concept{Pairwise variable selection for classification}
