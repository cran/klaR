\name{svmlight}
\alias{svmlight}
\alias{svmlight.default}
\alias{svmlight.formula}
\alias{svmlight.matrix}
\alias{svmlight.data.frame}
\alias{svmlight.file}
\title{Interface to SVMlight}
\description{
Function to call SVMlight from R for classification. Multiple group
classification is done with the one-against-rest partition of data.  
}
\usage{
svmlight(x, ...)

\method{svmlight}{default}(x, grouping, temp.dir = NULL, pathsvm = NULL, del = TRUE,
    svm.options = NULL, prior = NULL, out = FALSE,...)
\method{svmlight}{data.frame}(x, ...)
\method{svmlight}{matrix}(x, grouping, ..., subset, na.action = na.fail)
\method{svmlight}{formula}(formula, data = NULL, ..., subset, na.action = na.fail)
}

\arguments{
  \item{x}{matrix or data frame containing the explanatory variables 
    (required, if \code{formula} is not given).}
  \item{grouping}{factor specifying the class for each observation 
    (required, if \code{formula} is not given).}
  \item{formula}{formula of the form \code{groups ~ x1 + x2 + ...}. 
    That is, the response is the grouping factor and the right hand side specifies the (non-factor) discriminators.}
  \item{data}{Data frame from which variables specified in \code{formula} are preferentially to be taken.}
  \item{temp.dir}{directory for temporary files.}
  \item{pathsvm}{Path to SVMlight binaries (required, if path os unknown by the OS).}
  \item{del}{Logical: whether to delete temporary files}
  \item{svm.options}{Optional parameters to SVMlight.
  
    For further details see: \dQuote{How to use} on \url{http://svmlight.joachims.org/}.
  }
  \item{prior}{A Priori probabilities of classes.}
  \item{out}{Logical: whether SVMlight output ahouild be printed on console 
    (only for Windows OS.)}
  \item{subset}{An index vector specifying the cases to be used in the
    training sample. (Note: If given, this argument must be named.)}
  \item{na.action}{specify the action to be taken if \code{NA}s are
    found. The default action is for the procedure to fail. An
    alternative is \code{\link{na.omit}}, which leads to rejection of cases with
    missing values on any required variable. (Note: If given, this
    argument must be named.) }
  \item{...}{}
  
}
\details{
  Function to call SVMlight from R for classification. Multiple group
  classification is done with the one-against-rest partition of data.
  SVMlight is an implementation of Vapnik's Support Vector Machine. It
  is written in C by Thorsten Joachims. On the homepage (see below) the
  source-code and several binaries for SVMlight are available. If more
  then two classes are given the SVM is learned by the one-against-rest
  scheme. That means that each class is trained against the other K-1
  classes. The class with the highest decision function in the SVM
  wins. So K SVMs have to be learned.
}
\value{
  A list containing the function call and the result of SVMlight. 
}
\references{\url{http://svmlight.joachims.org/}}
\author{Karsten Luebke, \email{luebke@statistik.uni-dortmund.de}}
\seealso{\code{\link{predict.svmlight}},\code{\link[e1071]{svm}},
}
\examples{
\dontrun{
## Only works if the svmlight binaries are in the path.
data(iris)
x <- svmlight(Species ~ ., data=iris)
}
}
\keyword{classif}
