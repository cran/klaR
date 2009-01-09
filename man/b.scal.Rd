\name{b.scal}
\alias{b.scal}

\title{Calculation of beta scaling parameters}
\description{
Calculates the scaling parameter for \code{\link{betascale}}.
}
\usage{
b.scal(member, grouping, dis = FALSE, eps = 1e-04)
}
\arguments{
  \item{member}{Membership values of an argmax classification method. 
                Eg. posterior probabilities of \code{\link[MASS]{lda}}. 
                Row-wise values must sum up to 1 and must be in the interval [0,1].}
  \item{grouping}{Class vector.}
  \item{dis}{Logical, whether to optimize the dispersion parameter in \code{\link{pbeta}}.}
  \item{eps}{Minimum variation of membership values. If variance is smaller than \code{eps},
     the values are treated as one point.}
}
\details{
With \code{\link{betascale}} and \code{b.scal}, membership values of an argmax classifier 
are scaled in such a way, that the mean membership value of those values which are assigned 
to each class reflect the mean correctness rate of that values.
This is done via \code{\link{qbeta}} and \code{\link{pbeta}} with the appropriate shape parameters. 
If \code{dis} is \code{TRUE}, it is tried that the variation of membership values 
is optimal for the accuracy relative to the correctness rate. 
If the variation of the membership values is less than \code{eps},
they are treated as one point and shifted towards the correctness rate.
}
\value{
  A list containing
  \item{model}{Estimated parameters for \code{\link{betascale}}.}
  \item{eps}{Value of \code{eps} from the call.}
  \item{member}{Scaled membership values.}
}
\references{Garczarek, Ursula Maria (2002): Classification rules in standardized partition spaces.
            Dissertation, University of Dortmund. 
            URL \url{http://hdl.handle.net/2003/2789}}
\author{Karsten Luebke, Uwe Ligges}
\seealso{\code{\link{betascale}}, \code{\link{e.scal}}}
\examples{
library(MASS)
data(B3)
pB3 <- predict(lda(PHASEN ~ ., data = B3))$posterior
pbB3 <- b.scal(pB3, B3$PHASEN, dis = TRUE)
ucpm(pB3, B3$PHASEN)
ucpm(pbB3$member, B3$PHASEN)
}
\keyword{classif}
\concept{Beta scaling}
\concept{Classification}
\concept{Posterior probabilities}
