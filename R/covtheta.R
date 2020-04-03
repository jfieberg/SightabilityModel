#' Estimates var/cov matrix of inflation factors (1/prob detection) using a
#' non-parametric bootstrap.
#' 
#' Estimates var/cov matrix of inflation factors (1/prob detection) using a
#' non-parametric bootstrap.  Called by function Sight.Est if Vm.boot = TRUE.
#' 
#' 
#' @param total Number of animals in each independently sighted group
#' @param srates Plot sampling probability (associated with the independently
#' observed animal groups)
#' @param stratum Stratum identifiers (associated with the independently
#' observed animal groups)
#' @param subunit Plot ID (associated with the independently observed animal
#' groups)
#' @param covars Matrix of sightability covariates (associated with the
#' independently observed animal groups)
#' @param betas Logistic regression parameter estimates (from fitted
#' sightability model)
#' @param varbetas Estimated variance-covariance matrix for the logistic
#' regression parameter estimates (from fitted sightability model)
#' @param nboots Number of bootstrap resamples.
#' @return \item{smat }{Estimated variance-covariance matrix for the inflation
#' factors theta = (1/probability of detection).  This is an n.animal x
#' n.animal matrix. }
#' @author John Fieberg
#' @seealso \code{\link{Sight.Est}}
#' @keywords methods
#' @export covtheta
covtheta <-
function(total, srates, stratum, subunit, covars, betas, varbetas, nboots){

  #  Arguments (the first four should have one observation for each individual in the obs. dataset
  #  total = total number of animals in the observed group
  #  srates = sampling rate for the stratum (for each observed animal)
  #  subunit = subunit identifier 
  #  covars = matrix of covariates used in the sightability model
  #  beta = vector of parameter estimates from logistic model
  #  varbeta = var/cov matrix for above parameter estimates
  # nboots = number of bootstrap replicates
  
  # Check on length of vector arguments
    n <- length(total)
    if(length(srates) != n) {stop("Srates vector needs to be same dimension as total vector")}
    if(length(stratum) != n) {stop("Stratum vector needs to be same dimension as total vector")}
    if(length(srates) != n) {stop("Subunit vector needs to be same dimension as total vector")}

  # Make sure strata are numbered 1 to h
    stratum <- as.numeric(factor(stratum))

  # Form xmatrix
    xdata <- as.matrix(cbind(rep(1, n), covars))
    ncovars <- dim(xdata)
    nxdata <- nrow(xdata)
    if(ncovars[1] != n) {stop("Covariate matrix needs to be same length as total vector")}  
      
  # Make sure beta is a matrix of dim ncovar x 1 
    thetas <- matrix(NA, nboots, n)
    for(j in 1:nboots){
      tempbet <- betas[j, ]
      tempvarbeta <- varbetas[j, , ]
      
  #  Estimate theta for each animal = 1/prob(detection)
      thetas[j, ] <- 1+exp(-xdata%*%tempbet-diag(xdata%*%tempvarbeta%*%t(xdata))/2)  
  } 
    
    
  # Use empirical variance for Cov(thetas)
    smat <- cov(thetas) 
    return(smat)
}
