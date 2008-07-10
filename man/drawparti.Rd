\name{drawparti}
\alias{drawparti}

\title{Plotting the 2-d partitions of classification methods}
\description{
Plot showing the classification of observations based on 
classification methods (e.g. \code{lda}, \code{qda}) for two variables. 
Moreover, the classification borders are displayed and the apparent error rates are given in each title.
}
\usage{
drawparti(grouping, x, y, method = "lda", prec = 100, xlab = NULL, 
    ylab = NULL, col.correct = "black", col.wrong = "red", 
    col.mean = "black", col.contour = "darkgrey", 
    gs = as.character(grouping), pch.mean = 19, cex.mean = 1.3, 
    print.err = 0.7, legend.err = FALSE, legend.bg = "white",
    imageplot = TRUE, image.colors = cm.colors(nc), 
    plot.control = list(), ...)        
}

\arguments{
  \item{grouping}{factor specifying the class for each observation.}
  \item{x}{first explanatory vector.}
  \item{y}{second explanatory vector.}
  \item{method}{the method the classification is based on, currently supported are:
      \code{\link[MASS]{lda}}, \code{\link[MASS]{qda}}, \code{\link[rpart]{rpart}}, \code{\link[e1071]{naiveBayes}}, 
      \code{\link{rda}}, \code{\link{sknn}} and \code{\link{svmlight}}}.
  \item{prec}{precision used to draw the classification borders (the higher the more precise; default: 100).}
  \item{xlab}{a title for the x axis.}
  \item{ylab}{a title for the y axis.}
  \item{col.correct}{color for correct classified objects.}
  \item{col.wrong}{color for wrong classified objects.}
  \item{col.mean}{color for class means (only for methods \code{lda} and \code{qda}).}
  \item{col.contour}{color of the contour lines (if \code{imageplot = FALSE}).}
  \item{gs}{group symbol (plot character), must have the same length as \code{grouping}.}            
  \item{pch.mean}{plot character for class means (only for methods \code{lda} and \code{qda}).}
  \item{cex.mean}{character expansion for class means (only for methods \code{lda} and \code{qda}).}
  \item{print.err}{character expansion for text specifying the apparent error rate. 
    If \code{print.err = 0}, nothing is printed.}
  \item{legend.err}{logical; whether to plot the apparent error rate above the plot (if \code{FALSE}),
    or into a legend into the upper right corner of the plot (if \code{TRUE}). 
    This argument is ignored, if \code{print.err = 0}, i.e. if no error rate is printed.}
  \item{legend.bg}{Backgound colour to use for the legend.}
  \item{imageplot}{logical; whether to use an \code{\link{image}} plot or \code{\link{contour}} lines.}
  \item{image.colors}{colors used for the \code{imageplot}, if TRUE.}
  \item{plot.control}{A list containing further arguments passed to the underlying plot functions.}
  \item{...}{Further arguments passed to the classification \code{method}.}
}
\seealso{\code{\link{partimat}}}
\author{Karsten Luebke, \email{luebke@statistik.tu-dortmund.de}, Uwe Ligges, Irina Czogiel}
\keyword{classif}
\keyword{dplot}
\concept{Vizualizing classification results}
