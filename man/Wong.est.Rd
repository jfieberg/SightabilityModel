\name{Wong.est}
\alias{Wong.est}
\title{
Sightability estimate with variance components estimator from Wong (1996)
}
\description{
Estimates population size, with variance estimated using Wong's (1996) estimator.  This function will usualled be called by Sight.Est function (but see details).
}
\usage{
Wong.est(total, srates, nh, Nh, stratum, subunit, covars, beta, varbeta, smat = NULL)
}
\arguments{
  \item{total}{
 Number of animals in each independently sighted group
}
  \item{srates}{
  Vector of plot-level sampling probabilities (same dimension as \code{total}).
}
  \item{nh}{
  Number of sample plots in each stratum  

}
  \item{Nh}{
  Number of population plots in each stratum  

}
  \item{stratum}{
   Stratum identifiers (associated with the independently observed animal groups)
}
  \item{subunit}{
  Plot ID (associated with the independently observed animal groups)
}
  \item{covars}{
  Matrix of sightability covariates (associated with the independently observed animal groups)
}
  \item{beta}{
  Logistic regression parameter estimates (from fitted sightability model)
}
  \item{varbeta}{
  Estimated variance-covariance matrix for the logistic regression parameter estimates (from fitted sightability model)
}
  \item{smat}{
   Estimated variance-covariance matrix for the inflation factors (1/probability of detection).  This is an n.animal x n.animal matrix, and is usually calculated within the Wong.est function. 
   Non-null values can be passed to the function (e.g., if a bootstrap is used to estimate uncertainty due to the estimated detection parameters).
}
}
\details{
  This function is called by Sight.Est, but may also be called directly by the user (e.g., in cases where the original sighability [test trial] data are not available, but 
  the parameters and var/cov matrix from the logistic regression model is available in the literature).
}
\value{
  \item{tau.hat}{Sightability estimate of population size, tau^}
  \item{Vartot}{Estimated variance of tau^}
  \item{VarSamp}{Estimated variance component due to sampling aerial units}
  \item{VarSight}{Estimated variance component due to sighting process (i.e., series of binomial rv for each animal group)}
  \item{Varmod}{Estimated variance component due to estimating detection probabilities using test trial data}
}
\references{
    Rice CG, Jenkins KJ, Chang WY (2009).  Sightability Model for Mountain Goats." The Journal of Wildlife Management, 73(3), 468- 478.
  
  Steinhorst, R. K., and M.D. Samuel. (1989).  Sightability adjustment methods for aerial surveys of wildlife populations.  Biometrics 45:415-425.

  Wong, C. (1996).  Population size estimation using the modified Horvitz-Thompson estimator with estimated sighting probabilities.  Dissertation, Colorado State University, Fort Collins, USA.
}
\author{
John Fieberg
}

\seealso{
\code{\link{Sight.Est}}, \code{\link{SS.est}}
}
\keyword{methods}