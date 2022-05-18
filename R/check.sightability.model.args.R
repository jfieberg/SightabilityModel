#' Check the sightability model arguments for consistency
#' 
#' @param data Data.frame containing covariates for sightability model
#' @param sight.model Formula with sightability model
#' @param sight.beta Parameter estimates (from fitted sightability model
#' @param sight.beta.cov Estimated variance-covariance matrix for  parameter estimates
#'        from fitted sightability model.
#' @return Error condition or invisible
#' @template author 
#' @importFrom stats model.matrix
#' @importFrom formula.tools rhs.vars
#' @keywords methods
#' @examples
#' sightability.table <- data.frame(VegCoverClass=1:5)
#' sight.beta <- c(4.2138, -1.5847)
#' sight.beta.cov <- matrix(c(0.7821634, -0.2820000,-0.2820000,  0.1114892), nrow=2)
#' check.sightability.model.args( sightability.table, 
#'                                ~VegCoverClass, 
#'                                sight.beta, 
#'                                sight.beta.cov)

#' \dontrun{
#' check.sightability.model.args( sightability.table, 
#'                               ~VegCoverClass2, 
#'                               sight.beta,
#'                               sight.beta.cov)
#' check.sightability.model.args( sightability.table, 
#'                                ~VegCoverClass,
#'                                 sight.beta[1],
#'                                sight.beta.cov)
#' }
#' 
#' @export check.sightability.model.args
#' 
check.sightability.model.args <- function(data, sight.model, sight.beta, sight.beta.cov){
  # Check the sightability model arguments for consistency
  #
  # data - data frame containing covariates in the sightability model
  # sight.model- formula for sightability model
  # sight.beta - estimated coefficients for sightability model
  # sight.beta.cov - variance/covariance matrix for beta coefficients

  if(!is.data.frame(data))stop("Data must be a data frame")
  
  
  if(is.null(sight.model))stop("A sightability formula must be specified")
  if(!plyr::is.formula(sight.model))stop("A sightability model must be a valid formula")
  xvars <- formula.tools::rhs.vars(sight.model)
  
  if( !all(xvars %in% names(data)))stop("Sightability model uses variables not in data frame")
  
  # check if any missing values
  if( any(is.na(data[,xvars])))stop("Sightability model does not allow missing values in variables in sightability model")
  
  # make sure beta coefficients are specified and same length as model
  if(is.null(sight.beta))stop("You must specify beta coefficients for sightability model")
  if(!is.numeric(sight.beta))stop("Beta coefficients must be numeric")
  
  # see if variables in data frame match the beta coefficients et 
  design.mat <- stats::model.matrix(sight.model, data)
  
  if(length(sight.beta) != ncol(design.mat))stop("Length of beta coefficients doesn't match sightability model")
  
  # make sure that beta.vcov is right size
  if(is.null(sight.beta.cov))stop("You must specify beta covariance matrix for sightability model")
  if(!is.numeric(sight.beta.cov))stop("Beta covariance must be numeric")
  if(!is.matrix(sight.beta.cov))stop("Beta covariance must be a matrix")
  if(nrow(sight.beta.cov) != length(sight.beta) | ncol(sight.beta.cov) != length(sight.beta))
    stop("Beta covariance must be square and same size as beta coefficients")
  
  invisible()
}

# check.sightability.model.args( sightability.table, ~VegCoverClass2, sight.beta, sight.beta.cov)
# check.sightability.model.args( sightability.table, ~VegCoverClass, sight.beta[1], sight.beta.cov)
