#' Compute the detection probability given a sightability model
#' 
#' @param data Data.frame containing covariates for sightability model
#' @param sight.model Formula with sightability model
#' @param sight.beta Parameter estimates (from fitted sightability model
#' @param sight.beta.cov Estimated variance-covariance matrix for  parameter estimates
#'        from fitted sightability model.
#' @param check.args Should the sightability model arguments be checked for consistency/
#' @return Vector of detection probabilities
#' @template author 
#' @seealso \code{\link{compute.SCF}}
#' @importFrom stats model.matrix
#' @keywords methods
#' @examples
#' sightability.table <- data.frame(VegCoverClass=1:5)
#' sight.beta <- c(4.2138, -1.5847)
#' sight.beta.cov <- matrix(c(0.7821634, -0.2820000,-0.2820000,  0.1114892), nrow=2)
#' sightability.table$detect.prob <- compute.detect.prob( sightability.table, 
#'                                                       ~VegCoverClass, 
#'                                                       sight.beta, 
#'                                                       sight.beta.cov)
#' sightability.table$SCF         <- compute.SCF        ( sightability.table,
#'                                                        ~VegCoverClass, 
#'                                                        sight.beta, 
#'                                                        sight.beta.cov)
#' sightability.table
#' #"Note that the SCF != 1/detect.prob because of correction terms for covariance of beta.terms"
#' 
#' @export compute.detect.prob
#' 
compute.detect.prob <- function(data, sight.model, sight.beta, sight.beta.cov, 
                        check.args=FALSE){
  # Compute the probability of detection from a sightability model
  # See 
  #   Fieberg, J.R. (2012). 
  #   Estimating Population Abundance Using Sightability Models: R SightabilityModel Package. 
  #   Journal of Statistical Software, 51(9), 1-20. 
  #   https://doi.org/10.18637/jss.v051.i09
  #
  # data - data frame containing covariates in the sightability model
  # sight.mode - formula for sightability model
  # sight.beta - estimated coefficients for sightability model
  # sight.beta.cov - variance/covariance matrix for beta coefficients
  # check.args - should the arguments be checked for validity
  
  if(check.args){
    check.sightability.model.args(data, sight.model, sight.beta, sight.beta.cov)
  }
  
  # Get the design matrix
  dmat <- stats::model.matrix(sight.model, data=data)
  #browser()
  # compute the terms with optional adjustments
  part1 <- -dmat%*%sight.beta
  prob <- 1/(1+exp(part1))
  as.vector(prob)
  
}
