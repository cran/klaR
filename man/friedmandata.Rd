\name{friedman.data}
\alias{friedman.data}
\title{Friedman's classification benchmark data}
\description{Function to generate 3-class classification benchmarking data
  as introduced by J.H. Friedman (1989)}
\usage{
friedman.data(setting = 1, p = 6, samplesize = 40, asmatrix = FALSE)
}
\arguments{
  \item{setting}{the problem setting (integer 1,2,...,6).}
  \item{p}{number of variables (6, 10, 20 or 40).}
  \item{samplesize}{sample size (number of observations, >=6).}
  \item{asmatrix}{if \code{TRUE}, results are returned as a matrix,
  otherwise as a data frame (default).}
}
\details{
  When J.H. Friedman introduced the Regularized Discriminant Analysis
  (\code{\link[klaR]{rda}}) in 1989, he used artificially generated data
  to test the procedure and to examine its performance in comparison to
  Linear and Quadratic Discriminant Analysis
  (see also \code{\link[MASS]{lda}} and \code{\link[MASS]{qda}}).

  6 different settings were considered to demonstrate potential strengths
  and weaknesses of the new method:
  \enumerate{
    \item equal spherical covariance matrices,
    \item unequal spherical covariance matrices,
    \item equal, highly ellipsoidal covariance matrices with mean
        differences in low-variance subspace,
    \item equal, highly ellipsoidal covariance matrices with mean
        differences in high-variance subspace,
    \item unequal, highly ellipsoidal covariance matrices with zero mean
        differences and
    \item unequal, highly ellipsoidal covariance matrices with nonzero mean
        differences.
  }

  For each of the 6 settings data was generated with 6, 10, 20 and 40
  variables.

  Classification performance was then measured by repeatedly creating
  training-datasets of 40 observations and estimating the
  misclassification rates by test sets of 100 observations.

  The number of classes is always 3, class labels are assigned randomly
  (with equal probabilities) to observations, so the contributions of
  classes to the data differs from dataset to dataset. To make sure
  covariances can be estimated at all, there are always at least two
  observations from each class in a dataset.
}
\value{Depending on \code{asmatrix} either a data frame or a matrix with
  \code{samplesize} rows and \code{p+1} columns, the first column containing
  the class labels, the remaining columns being the variables.
}
\references{Friedman, J.H. (1989): Regularized Discriminant Analysis.
  In: \emph{Journal of the American Statistical Association} 84, 165-175.}
\author{Christian Röver, \email{roever@statistik.uni-dortmund.de}}
\seealso{
  \code{\link[klaR]{rda}}
}
\examples{
# Reproduce the 1st setting with 6 variables.
# Error rate should be somewhat near 9 percent.
training <- friedman.data(1, 6, 40)
x <- rda(class ~ ., data = training, gamma = 0.74, lambda = 0.77)
test <- friedman.data(1, 6, 100)
y <- predict(x, test[,-1])
errormatrix(test[,1], y$class)
}
\keyword{multivariate}
\concept{Friedman's classification benchmark data}
\concept{Regularized Discriminant Analysis}
