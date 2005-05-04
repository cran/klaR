\name{errormatrix}
\alias{errormatrix}
\encoding{latin1}
\title{Tabulation of prediction errors by classes}
\description{
Cross-tabulates true and predicted classes 
with the option to show relative frequencies.
}
\usage{
errormatrix(true, predicted, relative = FALSE)
}
\arguments{
\item{true}{Vector of true classes.}
\item{predicted}{Vector of predicted classes.}
\item{relative}{Logical. If \code{TRUE} rows are normalized
to show relative frequencies (see below).}
}
\details{
Given vectors of true and predicted classes, a (symmetric) 
table of misclassifications is constructed.

Element [i,j] shows the number of objects of class i that 
were classified as class j; so the main diagonal shows the 
correct classifications. The last row and column show the
corresponding sums of misclassifications, the lower right 
element is the total sum of misclassifications.

If \sQuote{\code{relative}} is \code{TRUE}, the \emph{rows} are 
normalized so they show relative frequencies instead. The 
lower right element now shows the total error rate, and the
remaining last row sums up to one, so it shows \dQuote{where the
misclassifications went}.
}
\note{
Concerning the case that \sQuote{\code{relative}} is \code{TRUE}:

If a prior distribution over the classes is given, the 
misclassification rate that is returned as the lower right 
element (which is only the fraction of misclassified 
\emph{data}) is not an estimator for the expected 
misclassification rate.

In that case you have to multiply the individual error rates for 
each class (returned in the last column) with the corresponding
prior probabilities and sum these up (see example below).

Both error rate estimates are equal, if the fractions of classes 
in the data are equal to the prior probabilities.
}
\value{
A (named) matrix.
}
\author{Christian Röver, \email{roever@statistik.uni-dortmund.de}}
\seealso{\code{\link{table}}}
\examples{
data(iris)
library(MASS)
x <- lda(Species ~ Sepal.Length + Sepal.Width, data=iris)
y <- predict(x, iris)

# absolute numbers: 
errormatrix(iris$Species, y$class)

# relative frequencies: 
errormatrix(iris$Species, y$class, relative = TRUE)

# percentages: 
round(100 * errormatrix(iris$Species, y$class, relative = TRUE), 0)

# expected error rate in case of class prior: 
indiv.rates <- errormatrix(iris$Species, y$class, relative = TRUE)[1:3, 4]
prior <- c("setosa" = 0.2, "versicolor" = 0.3, "virginica" = 0.5)
total.rate <- t(indiv.rates) \%*\% prior
total.rate
}
\keyword{multivariate}
\concept{Visualizing classification results}
\concept{Confusion matrix}
