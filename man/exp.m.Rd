\name{exp.m}
\alias{exp.m}
\docType{data}
\title{
Experimental (test trials) data set used to estimate detection probabilities for moose in MN
}
\description{
Experimental (test trials) data set used to estimate detection probabilities for moose in MN
}
\usage{data(exp.m)}
\format{
  A data frame with 124 observations on the following 4 variables.
  \describe{
    \item{\code{year}}{year of the experimental survey (test trial)} 
    \item{\code{observed}}{Boolean variable (=1 if moose was observed and 0 otherwise)}
    \item{\code{voc}}{measurement of visual obstruction}
    \item{\code{grpsize}}{group size (number of observed moose in each independently sighted group)}
  }
}
\references{
  Giudice, J H. and Fieberg, J. and Lenarz, M. S.  2012. Spending Degrees of Freedom in a Poor Economy: A Case Study of Building a Sightability Model for Moose in Northeastern Minnesota.
   Journal of Wildlife Management 76(1):75-87.}
\examples{
data(exp.m)
exp.m[1:5,]
}
\keyword{datasets}
