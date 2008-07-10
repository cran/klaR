\name{countries}
\alias{countries}
\docType{data}
\title{Socioeconomic data for the most populous countries.}
\description{
Socioeconomic data for the most populous countries.
}
\usage{data(countries)}
\format{
  A data frame with 42 observations on the following 7 variables.
  \describe{
    \item{Country}{name of the country.}
    \item{Popul}{population.}
    \item{PopDens}{population density.}
    \item{GDPpp}{GDP per inhabitant.}
    \item{LifeEx}{mean life expectation}
    \item{InfMor}{infant mortality}
    \item{Illit}{illiteracy rate}
  }
}

\source{
CIA World Factbook \url{http://www.cia.gov/cia/publications/factbook/}
}

\examples{
data(countries)
summary(countries)
}
\keyword{datasets}
