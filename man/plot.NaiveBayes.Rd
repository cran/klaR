\name{plot.NaiveBayes}
\alias{plot.NaiveBayes}
\title{Naive Bayes Plot}
\description{
  Visualizes the marginal probabilities of predictor variables 
  given the class.
}
\usage{
plot.NaiveBayes(x, vars, n = 1000, legendplot = TRUE, 
    lty = 1:length(x$apriori), col = rainbow(length(x$apriori)), 
    ylab = "Density", main = "Naive Bayes Plot", ...)
    }

\arguments{
  \item{x}{an object of \code{class} \code{\link{NaiveBayes}}}
  \item{vars}{variables to be plotted. If \code{missing} all predictor variables are plotted.}
  \item{n}{Number points used to plot the density line.}
  \item{legendplot}{print a \code{\link{legend}}?}
  \item{lty}{\code{lty} of the different classes.}
  \item{col}{\code{color} for the different classes.}
  \item{ylab}{label for y-axis.}
  \item{main}{title of the plots.}
  \item{\dots}{furhter arguments passed to the underlying plot functions.}
}
\details{
 For metric variables the estimated density is plotted. For categorial variables \code{\link{mosaicplot} is called.}
}

\author{Karsten Luebke}

\seealso{\code{\link{NaiveBayes}}}
\examples{
data(iris)
mN<-NaiveBayes(Species~.,data=iris)
plot(mN)

mK<-NaiveBayes(Species~.,data=iris,usekernel=TRUE)
plot(mK)
}
\keyword{classif}
\keyword{dplot}
