\name{hmm.sop}
\alias{hmm.sop}
\title{Calculation of HMM Sum of Path}
\description{
A Hidden Markov Model for the classification of states in a time series. 
Based on the transition probabilities and the so called emission probabilities
(\eqn{p(class|x)}{p(class|x)}) the \sQuote{prior probablilities} of states (classes) in time period \emph{t} 
given all past information in time period \emph{t} are calculated.
}
\usage{
hmm.sop(sv, trans.matrix, prob.matrix)
}
\arguments{
  \item{sv}{state at time 0}
  \item{trans.matrix}{matrix of transition probabilities}
  \item{prob.matrix}{matrix of \eqn{p(class|x)}{p(class|x)}}
}

\value{
Returns the \sQuote{prior probablilities} of states. 
}

\references{Garczarek, Ursula Maria (2002): Classification rules in standardized partition spaces.
            Dissertation, University of Dortmund. 
            URL \url{http://hdl.handle.net/2003/2789}}

\author{Daniel Fischer, Reinald Oetsch}


\seealso{\code{\link{calc.trans}}}
\examples{
library(MASS)
data(B3)
trans.matrix <- calc.trans(B3$PHASEN)

# Calculate posterior prob. for the classes via lda
prob.matrix <- predict(lda(PHASEN ~ ., data = B3))$post
errormatrix(B3$PHASEN, apply(prob.matrix, 1, which.max))
prior.prob <- hmm.sop("2", trans.matrix, prob.matrix)
errormatrix(B3$PHASEN, apply(prior.prob, 1, which.max))
}
\keyword{ts}
\keyword{classif}
\concept{Hidden Markov Model}
\concept{Classification in time series}
\concept{HMM}
