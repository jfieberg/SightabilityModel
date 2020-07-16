\name{obs.m}
\alias{obs.m}
\docType{data}
\title{
MN moose survey data
}
\description{
Operational survey data for moose in MN (during years 2004-2007).  Each record corresponds to an 
independently sighted group of moose, with variables that capture
individual covariates (used in the detection model) as well as plot-level information (stratum identifier, sampling probability, etc). 
}
\usage{data(obs.m)}
\format{
  A data frame with 805 observations on the following 11 variables.
  \describe{
    \item{\code{year}}{year of survey}
    \item{\code{stratum}}{stratum identifier}
    \item{\code{subunit}}{sample plot ID}
    \item{\code{total}}{number of moose observed}
    \item{\code{cows}}{number of cows observed}
    \item{\code{calves}}{number of calves observed}
    \item{\code{bulls}}{number of bulls observed}
    \item{\code{unclass}}{number of unclassified animals observed (could not identify sex/age class)}
    \item{\code{voc}}{ measurement of visual obstruction}
    \item{\code{grpsize}}{group size (cluster size)}
  }
}
\references{
 Giudice, J H. and Fieberg, J. and Lenarz, M. S.  2012. Spending Degrees of Freedom in a Poor Economy: A Case Study of Building a Sightability Model for Moose in Northeastern Minnesota.
   Journal of Wildlife Management 76(1):75-87.
}
\examples{
data(obs.m) 
obs.m[1:5, ]
}
\keyword{datasets}
