\name{rerange}
\alias{rerange}
\title{Linear transformation of data}
\description{
 The function performs a linear transformation of the data, such that
 afterwards \code{range(data)=c(theMin,theMax)}.
}
\usage{
rerange(data, min.goal = 0, max.goal = 1, min.data =
    min(data), max.data = max(data), center = NA)
}
\arguments{
  \item{data}{ vector with the \code{data} to transform }
  \item{min.goal}{ new minimum value }
  \item{max.goal}{ new maximum value }
  \item{min.data}{ old minimum value }
  \item{max.data}{ old maximum value }
  \item{center}{which old value should become the new center \code{(max.goal + min.goal) / 2)}}
}
\value{
  vector with the transformed data
}
\author{Dominik Reusser}
\examples{
   klaR:::rerange(data=1:20)
   klaR:::rerange(data=1:30, center=5)
}
\keyword{internal}
