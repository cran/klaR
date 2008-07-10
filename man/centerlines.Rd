\name{centerlines}
\alias{centerlines}
\title{Lines from classborders to the center}
\description{
Function which constructs the lines from the borders between 
two classes to the center. To be used in connection with
\code{\link{triplot}} and \code{\link{quadplot}}.
}
\usage{
centerlines(n)
}
\arguments{
  \item{n}{number of classes. Meaningful are \code{3} or \code{4}.}
}

\value{
a matrix with \code{n}-columns. 
}
\author{Karsten Luebke}

\seealso{\code{\link{triplot}}, \code{\link{quadplot}}}
\examples{
centerlines(3)
centerlines(4)
}
\keyword{classif}
\keyword{dplot}
\concept{Visualizing Classification Performance Measures}
\concept{Barycentric plots}
