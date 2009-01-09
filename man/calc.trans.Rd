\name{calc.trans}
\alias{calc.trans}
\title{Calculation of transition probabilities}
\description{
Function to estimate the probabilities of a time series to stay or change the state.
}
\usage{
calc.trans(x)
}
\arguments{
  \item{x}{(factor) vector of states}
}
\details{
To estimate the transition probabilities the empirical frequencies are counted.
}
\value{
The transition probabilities matrix. 
\code{x[i,j]} is the probability to change from state \code{i} to state \code{j}.
}
\author{Karsten Luebke}

\examples{
data(B3)
calc.trans(B3$PHASEN)
}
\keyword{ts}
