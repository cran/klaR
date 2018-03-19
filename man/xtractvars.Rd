\name{xtractvars}
\alias{xtractvars}
\title{Variable clustering based variable selection}
\description{Applies variable selection to data based on variable clusterings as resulting from \code{\link{corclust}} or  \code{\link[ClustVarLV]{CLV}}.}  
\usage{
xtractvars(object, data, thres = 0.5)
}
\arguments{
  \item{object}{Object of class \code{\link{cvtree}} applied to a \code{\link{corclust}} object or the \code{summary()} of a \code{clv} object as created by \code{\link[ClustVarLV]{CLV}}.}
  \item{data}{Data where variables are to be selected. Coloumn names must be identical to those used in corclust model.}
  \item{thres}{Maximum accepted average within cluster correlation for selection of a variable.}
}
\details{
Of each cluster the first variable is selected as well as all other variables with an average within cluster correlation below \code{thres}. 
}
\value{
The data is returned where unselected coloumns are removed.
}

\references{Roever, C. and Szepannek, G. (2005): Application of a genetic algorithm
to variable selection in fuzzy clustering. In C. Weihs and W. Gaul (eds), Classification 
- The Ubiquitous Challenge, 674-681, Springer.}


\author{Gero Szepannek}

\seealso{See also \code{\link{corclust}}, \code{\link{cvtree}} and  \code{\link[ClustVarLV]{CLV}}.}
\examples{
    data(B3)
    ccres <- corclust(B3)
    plot(ccres)
    cvtres <- cvtree(ccres, k = 3)
    newdata <- xtractvars(cvtres, B3, thres = 0.5) 
}

\keyword{cluster}
\keyword{multivariate}
\keyword{classif}
\keyword{manip}
\keyword{attribute}
