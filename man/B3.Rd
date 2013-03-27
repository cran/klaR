\name{B3}
\alias{B3}
\encoding{latin1}
\docType{data}
\title{West German Business Cycles 1955-1994}
\description{West German Business Cycles 1955-1994}
\usage{data(B3)}
\format{
  A data frame with 157 observations on the following 14 variables.
  \describe{
    \item{PHASEN}{a factor with levels \code{1} (upswing), 
        \code{2} (upper turning points), \code{3} (downswing), and \code{4} (lower turning points).}
    \item{BSP91JW}{GNP (y)}
    \item{CP91JW}{Private Consumption (y)}
    \item{DEFRATE}{Government deficit (percent of GNP)}
    \item{EWAJW}{Wage and salary earners (y)}
    \item{EXIMRATE}{Net exports as (percent of GNP)}
    \item{GM1JW}{Money supply M1 (y)}
    \item{IAU91JW}{Investment in equipment (y)}
    \item{IB91JW}{Investment in construction (y)}
    \item{LSTKJW}{Unit labor cost (y)}
    \item{PBSPJW}{GNP price deflator (y)}
    \item{PCPJW}{Consumer price index (y)}
    \item{ZINSK}{Short term interest rate (nominal)}
    \item{ZINSLR}{Long term interest rate (real)}
  }
  where (y) stands for \dQuote{yearly growth rates}.
  
  Note that years and corresponding year quarters are given in the row names of the data frame,
  e.g. \dQuote{1988,3} for the third quarter in 1988. 
}
\details{
The West German Business Cycles data (1955-1994) is analyzed by the project \emph{B3}
of the SFB475 (Collaborative Research Centre 
\dQuote{Reduction of Complexity for Multivariate Data Structures}),
supported by the Deutsche Forschungsgemeinschaft.
}
\source{
RWI (Rheinisch Westfälisches Institut für Wirtschaftsforschung), Essen, Germany.
}
\seealso{For benchmarking on this data see also \code{\link{benchB3}}}
\references{
Heilemann, U. and Münch, H.J. (1996): 
    West German Business Cycles 1963-1994: A Multivariate Discriminant Analysis.
    \emph{CIRET--Conference in Singapore, CIRET--Studien} 50.
}
\examples{
data(B3)
summary(B3)
}
\keyword{datasets}
