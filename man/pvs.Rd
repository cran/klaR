\name{pvs}
\alias{pvs}
\alias{pvs.default}
\alias{pvs.formula}
\alias{print.pvs}
\title{Pairwise variable selection for classification}
\description{
Pairwise variable selection for numerical data, allowing the use of different classifiers and different variable selection methods.
}
\usage{
pvs(x, ...)

\method{pvs}{default}(x, grouping, prior=NULL, method="lda", 
    vs.method=c("ks.test","stepclass","greedy.wilks"), niveau=0.05, 
    fold=10, impr=0.1, direct="backward", out=FALSE, ...)
    
\method{pvs}{formula}(formula, data = NULL, ...)
}
\arguments{
\item{x}{matrix or data frame containing the explanatory variables 
    (required, if \code{formula} is not given). x must consist of numerical data only. }
\item{formula}{A formula of the form \code{groups ~ x1 + x2 + ...}. 
    That is, the response is the grouping factor (the classes) and the right hand side 
    specifies the (numerical) discriminators. 
    Interaction terms are not supported.}
\item{data}{data matrix (rows=cases, columns=variables)}
\item{grouping}{class indicator vector (a factor)}
\item{prior}{prior probabilites for the classes. If not specified the prior probabilities will be set according to proportion in \dQuote{grouping}. If specified the order of prior 
        probabilities must be the same as in \dQuote{grouping}. }
\item{method}{character, name of classification function (e.g. \dQuote{\code{\link[MASS]{lda}}} (default)).}
\item{vs.method}{character, name of variable selection method. Must be one of \dQuote{\code{\link[stats]{ks.test}}} (default),
    \dQuote{\code{\link[klaR]{stepclass}}} or \dQuote{\code{\link[klaR]{greedy.wilks}}}. }
\item{niveau}{used niveau for \dQuote{\code{\link[stats]{ks.test}}}}
\item{fold}{parameter for cross-validation, if \dQuote{\code{\link[klaR]{stepclass}}} is chosen \sQuote{\code{vs.method}}}                                            
\item{impr}{least improvement of performance measure desired to include or exclude any variable (<=1), if \dQuote{\code{\link[klaR]{stepclass}}} is chosen \sQuote{\code{vs.method}} }  
\item{direct}{direction of variable selection, if \dQuote{\code{\link[klaR]{stepclass}}} is chosen \sQuote{\code{vs.method}}. 
        Must be one if \dQuote{\code{forward}}, \dQuote{\code{backward}} (default) or \dQuote{\code{both}}. }
\item{out}{indicator (logical) for textoutput during computation (slows down computation!), if \dQuote{\code{\link[klaR]{stepclass}}} is chosen \sQuote{\code{vs.method}} }
\item{...}{further parameters passed to classification function (\sQuote{\code{method}}) or variable selection method (\sQuote{\code{vs.method}}) }
}
\details{
The classification \dQuote{method} (e.g. \sQuote{\code{\link[MASS]{lda}}}) must have its own 
\sQuote{\code{predict}} method (like \sQuote{\code{\link[MASS]{predict.lda}}} for \sQuote{\code{lda}}) 
returns a list with an element \sQuote{\code{posterior}} containing the posterior probabilties. It must be able to deal with matrices as in \code{method(x, grouping, ...)}. 
Examples of such classification methods are \sQuote{\code{\link[MASS]{lda}}}, \sQuote{\code{\link[MASS]{qda}}}, \sQuote{\code{\link[klaR]{rda}}}, 
\sQuote{\code{\link[klaR]{NaiveBayes}}} or \sQuote{\code{\link[klaR]{sknn}}}.\\
For the classification methods \dQuote{\code{\link[e1071]{svm}}} and \dQuote{\code{\link[randomForest]{randomForest}}} there are special routines implemented, to make them work with \sQuote{\code{pvs}} method even though their \sQuote{\code{predict}} methods don't provide the demanded posteriors. However those two classfiers can not be used together with variable selection method \dQuote{\code{\link[klaR]{stepclass}}}.

\sQuote{\code{pvs}} performs a variable selection using the selection method chosen in \sQuote{\code{vs.method}} for each pair of classes in \sQuote{\code{x}}. 
Then for each pair of classes a submodel using \sQuote{\code{method}} is trained (using only the earlier selected variables for this class-pair).

If \sQuote{\code{method}} is \dQuote{\code{\link[stats]{ks.test}}}, then for each variable the empirical distribution functions of the cases of both classes are compared via \dQuote{\code{\link[stats]{ks.test}}}. Only variables with a p-values below \sQuote{\code{niveau}} are used for training the submodel for this pair of classes.

If \sQuote{\code{method}} is \dQuote{\code{\link[klaR]{stepclass}}} the variable selection is performed using the \dQuote{\code{\link[klaR]{stepclass}}} method.

If \sQuote{\code{method}} is \dQuote{\code{\link[klaR]{greedy.wilks}}} the variable selection is performed using Wilk's lambda criterion.
}
\value{
An object of class \sQuote{\code{pvs}} containing the following components:
\item{classes}{the classes in grouping}
\item{prior}{used prior probabilities}
\item{method}{name of used classification function}
\item{vs.method}{name of used function for variable selection}
\item{submodels}{containing a list of submodels. For each pair of classes there is a list element being another list of 3 containing the class-pair of this submodel, the selected variables 
        for the subspace of classes and the result of the trained classification function.}
\item{call}{the (matched) function call}
}
\author{Gero Szepannek, \email{szepannek@statistik.tu-dortmund.de}, Christian Neumann}
\references{
    \itemize{
        \item{}{Szepannek, G. and Weihs, C. (2006)  Variable Selection for Classification of More than Two 
            Classes Where the Data are Sparse. In \emph{From Data and Information Analysis to Kwnowledge Engineering.},
          eds Spiliopolou, M., Kruse, R., Borgelt, C., Nuernberger, A. and Gaul, W. pp. 700-708. Springer, Heidelberg.}
        \item{}{Szepannek, G. (2008): Different Subspace Classification - Datenanalyse, -interpretation, -visualisierung und 
          Vorhersage in hochdimensionalen Raeumen, ISBN 978-3-8364-6302-7, vdm, Saarbruecken.}
    }
}
\seealso{
\code{\link[klaR]{predict.pvs}} for predicting \sQuote{\code{pvs}} models and \code{\link[klaR]{locpvs}} for pairwisevariable selection in local models of several subclasses
}
\examples{
## Example 1: learn an "lda" model on the waveform data using pairwise variable 
## selection (pvs) using "ks.test" and compare it to using lda without pvs 

library("mlbench")
trainset <- mlbench.waveform(300) 
pvsmodel <- pvs(trainset$x, trainset$classes, niveau=0.05) # default: using method="lda"
## short summary, showing the class-pairs of the submodels and the selected variables
pvsmodel
 
testset <-  mlbench.waveform(500)
## prediction of the test data set: 
prediction <- predict(pvsmodel, testset$x)

## calculating the test error rate
1-sum(testset$classes==prediction$class)/length(testset$classes)
## Bayes error is 0.149

## comparison to performance of simple lda
ldamodel <- lda(trainset$x, trainset$classes)
LDAprediction <- predict(ldamodel, testset$x)

## test error rate
1-sum(testset$classes==LDAprediction$class)/length(testset$classes)


## Example 2: learn a "qda" model with pvs on half of the Satellite dataset, 
## using "ks.test" 

\donttest{
library("mlbench")
data("Satellite")

## takes few seconds as exact KS tests are calculated here:
model <- pvs(classes ~ ., Satellite[1:3218,], method="qda", vs.method="ks.test")
## short summary, showing the class-pairs of the submodels and the selected variables
model 

## now predict on the rest of the data set:
## pred <- predict(model,Satellite[3219:6435,]) # takes some time
pred <- predict(model,Satellite[3219:6435,], quick=TRUE) # that's much quicker

## now you can look at the predicted classes:
pred$class 
## or the posterior probabilities:
pred$posterior
}
}
\keyword{classif}
\keyword{multivariate}
\concept{Pairwise variable selection for classification}
