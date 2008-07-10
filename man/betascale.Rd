\name{betascale}
\alias{betascale}
\title{Scale membership values according to a beta scaling}
\description{
Performs the scaling for beta scaling learned by \code{\link{b.scal}}.
}
\usage{
betascale(betaobj, member)
}

\arguments{
  \item{betaobj}{A model learned by \code{\link{b.scal}}.}
  \item{member}{Membership values to be scaled.}
}
\details{
See \code{\link{b.scal}}.
}
\value{
A matrix with the scaled membership values.
}
\seealso{\code{\link{b.scal}}, \code{\link{e.scal}}}
\examples{
library(MASS)
data(B3)
pB3 <- predict(lda(PHASEN ~ ., data = B3))$posterior
pbB3 <- b.scal(pB3, B3$PHASEN)
betascale(pbB3)
}
\keyword{classif}
\concept{Beta scaling}
\concept{Classification}
\concept{Posterior probabilities}
