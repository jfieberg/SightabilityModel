% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{exp.m}
\alias{exp.m}
\title{Experimental (test trials) data set used to estimate detection probabilities
for moose in MN}
\format{
A data frame with 124 observations on the following 4 variables.
\describe{ 
  \item{year}{year of the experimental survey (test trial)}
  \item{observed}{Boolean variable (=1 if moose was observed and 0 otherwise)} 
  \item{voc}{measurement of visual obstruction}
  \item{grpsize}{group size (number of observed moose in each independently sighted group)}
 }
}
\description{
Experimental (test trials) data set used to estimate detection probabilities
for moose in MN
}
\examples{

data(exp.m)
exp.m[1:5,]

}
\references{
Giudice, J H. and Fieberg, J. and Lenarz, M. S.  2012. Spending
Degrees of Freedom in a Poor Economy: A Case Study of Building a
Sightability Model for Moose in Northeastern Minnesota.  Journal of Wildlife
Management 76(1):75-87.
}
\keyword{datasets}
