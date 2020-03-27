\name{vardiff}
\alias{vardiff}
\title{Function to estimate the variance of the difference between two population estimates
}
\description{
Function to estimate the variance of the difference between two population estimates formed using the same sightability model (to correct for detection).}
\usage{
vardiff(sight1, sight2)
}
\arguments{
  \item{sight1}{ Sightability model object for the first population estimate (formed by calling Sight.Est function) 
}
  \item{sight2}{ Sightability model object for the second population estimate (formed by calling Sight.Est function) 
}
}
\details{
  Population estimates constructed using the same sightability model will NOT be independent (they will typically exhibt positive covariance).  
  This function estimates the covariance due to using the same sightability model and subtracts it from the summed variance.
}
\value{
  \item{vardiff}{numeric = var(tau^[1])+var(tau^[2])-2*Cov(tau^[1],tau^[2])}
}
\author{
John Fieberg
}


\examples{

# Example using moose survey data 
  data(obs.m) # observational moose survey data
  data(exp.m) # experimental moose survey data
  data(sampinfo.m) # information on sampling rates
 
# Estimate population size in 2006 and 2007 
 sampinfo <- sampinfo.m[sampinfo.m$year == 2007, ]
 tau.2007 <- Sight.Est(observed ~ voc, odat = obs.m[obs.m$year == 2007, ], 
                         sdat = exp.m, sampinfo.m[sampinfo.m$year == 2007, ], 
                         method = "Wong", logCI = TRUE, alpha = 0.05, Vm.boot = FALSE) 
 tau.2006 <- Sight.Est(observed ~ voc, odat = obs.m[obs.m$year == 2006, ],
                         sdat = exp.m, sampinfo.m[sampinfo.m$year == 2006, ],
                         method = "Wong", logCI = TRUE, alpha = 0.05, Vm.boot = FALSE) 

# naive variance
  tau.2007$est[2]+tau.2006$est[2]

# variance after subtracting positvie covariance
  vardiff(tau.2007, tau.2006)
}
\keyword{methods}