\name{covtheta}
\alias{covtheta}
\title{
Estimates var/cov matrix of inflation factors (1/prob detection) using a non-parametric bootstrap.   
}
\description{
Estimates var/cov matrix of inflation factors (1/prob detection) using a non-parametric bootstrap.  Called by function Sight.Est if Vm.boot = TRUE.
}
\usage{
covtheta(total, srates, stratum, subunit, covars, betas, varbetas, nboots)
}
\arguments{
  \item{total}{Number of animals in each independently sighted group
}
  \item{srates}{Plot sampling probability (associated with the independently observed animal groups)
}
  \item{stratum}{Stratum identifiers (associated with the independently observed animal groups)
}
  \item{subunit}{Plot ID (associated with the independently observed animal groups)
}
  \item{covars}{Matrix of sightability covariates (associated with the independently observed animal groups)
}
  \item{betas}{Logistic regression parameter estimates (from fitted sightability model)
}
  \item{varbetas}{Estimated variance-covariance matrix for the logistic regression parameter estimates (from fitted sightability model)
}
  \item{nboots}{Number of bootstrap resamples.
}
}
\keyword{methods}

\value{
  \item{smat }{Estimated variance-covariance matrix for the inflation factors theta = (1/probability of detection).  This is an n.animal x n.animal matrix.
}
}
\author{
John Fieberg
}

\seealso{
 \code{\link{Sight.Est}}
}
