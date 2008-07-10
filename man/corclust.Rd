\name{corclust}
\alias{corclust}
\title{Function to identify groups of highly correlated variables 
for removing correlated features from the data for further analysis.}
\description{
A hierarchical clustering of variables using \code{hclust} is performed using 
1 - the absolute correlation as a distance measure between tow variables.  

}
\usage{
corclust(x, cl = NULL, mincor = NULL, prnt = FALSE, method = "complete")
}
\arguments{
  \item{x}{Either a data frame or a matrix consisting of numerical attributes.}
  \item{cl}{Optional vector of ty factor indicating class levels, if class specific correlations should to be considered.}
  \item{mincor}{Optional vector of degrees of correlation within a cluster of variables 
                that will be indicated in the plot by a line.}
  \item{prnt}{Logical indicating whether the matrix of distances should be printed.}
  \item{method}{Linkage to be used for clustering. Default is \code{complete} linkage.}
}
\details{
The main output consists in the tree visualization of the clustered variables. 
Each cluster consists of a set of correlated variables according to the chosen clustering criterion. 
The default criterion is \sQuote{\code{complete}}. This choice is meaningful as it represents the 
{\emph minimum absolute correlation} between all variables of a cluster.\\
Further proceeding would consist in chosing one variable of each cluster to obtain a 
subset of rather uncorrelated variables for further analysis.\\
If an additional class vector \code{cl} is given to the function for any two variables their minimum 
correlation over all classes is used.  
}
\value{
  \item{min.abs.cor}{Matrix of distances used for clustering containing 
                1 - the absolute correlation between any two variables.}
  \item{clustering}{Result object of the hierarchical clustering.}
}
\author{Gero Szepannek}

\seealso{See also \code{\link{hclust}}, for details on the clustering algorithm.}
\examples{
    data(iris)
    classes <- iris$Species
    variables <- iris[,1:4]
    corclust(variables, classes, mincor = 0.6) 
}

\keyword{cluster}
\keyword{multivariate}
\keyword{classif}
\keyword{manip}
\keyword{attribute}
