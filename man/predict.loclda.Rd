\name{predict.loclda}
\alias{predict.loclda}
\encoding{latin1}
\title{Localized Linear Discriminant Analysis (LocLDA)}
\description{
Classifies new observations using parameters determined by
the \code{loclda}-function.
}
\usage{
\method{predict}{loclda}(object, newdata, ...)
}
\arguments{
\item{object}{Object of class \code{loclda}.}
\item{newdata}{Data frame of cases to be classified.}
\item{\dots}{Further arguments are ignored.}
}
\value{
A list with components:
\item{class}{Vector (of class \code{factor}) of classifications.}
\item{posterior}{Posterior probabilities for the classes. 
    For details of computation see \code{\link{loclda}} 
    (+ normalization so posterior-values add up to 1 for each observation).}
\item{all.zero}{Vector (of class \code{integer}) indicating for which rows of
     \code{newdata} all corresponding posterior-values are \eqn{< 10^{-150}}{< 1e-150} before normalization. 
     Those observations are assigned to the class to whose (locally weighted) centroid they have the 
     lowest euclidian distance.}
}
\author{Marc Zentgraf (\email{marc-zentgraf@gmx.de}) and Karsten Luebke (\email{karsten.luebke@fom.de})}
\seealso{
  \code{\link{loclda}}
}
\examples{
data(B3)
x <- loclda(PHASEN ~ ., data = B3, subset = 1:80)
predict(x, B3[-(1:80),])
}
\keyword{multivariate}
\concept{Localized Linear Discriminant Analysis}
\concept{Linear Discriminant Analysis}
\concept{Localization}

