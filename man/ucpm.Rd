\name{ucpm}
\alias{ucpm}
\title{Uschi's classification performance measures}
\description{
Function to calculate the Correctness Rate, the Accuracy, the Ability to Seperate and the Confidence of 
a classification rule.
}
\usage{
ucpm(m, tc, ec = NULL)
}
\arguments{
  \item{m}{matrix of (scaled) membership values}
  \item{tc}{vector of true classes}
  \item{ec}{vector of estimated classes (only required if scaled membership values are used)}
}
\details{ 
\itemize{
\item The \emph{correctness rate} is the estimator for the correctness of a classification rule (1-error rate). 

\item The \emph{accuracy} is based on the euclidean distances between (scaled) membership vectors and the vectors 
representing the true class corner. These distances are standardized so that a measure of 1 is achieved
if all vectors lie in the correct corners and 0 if they all lie in the center. 

\item Analougously, the \emph{ability to seperate} is based on the distances between (scaled) membership 
vectors and the vector representing the corresponding assigned class corner. 

\item The \emph{confidence} is the mean of the membership values of the assigned classes. 
}

}
\value{
A list with elements:
  \item{CR}{Correctness Rate}
  \item{AC}{Accuracy}
  \item{AS}{Ability to Seperate}
  \item{CF}{Confidence}
  \item{CFvec}{Confidence for each (true) class}
}

\references{Garczarek, Ursula Maria (2002): Classification rules in standardized partition spaces.
            Dissertation, University of Dortmund. 
            URL \url{http://hdl.handle.net/2003/2789}}
            
\author{Karsten Luebke, \email{luebke@statistik.uni-dortmund.de}}
\examples{
library(MASS)
data(iris)
ucpm(predict(lda(Species ~ ., data = iris))$posterior, iris$Species)
}
\keyword{classif}
\concept{Classification Performance Measures}
\concept{Correctness Rate}
\concept{Accuracy}
\concept{Confidence}
\concept{Ability to Seperate}
