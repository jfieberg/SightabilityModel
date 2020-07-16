#' @param sight.formula A formula that identifies the model used
#'   to estimate sightability. For example \code{observed ~ VegCoverClass} would indicate
#'   that sightability is a function of the \code{VegCoverClass} variable in the survey 
#'   data. The left hand variable is arbitrary. The right hand variables must be present
#'   in the survey.data data frame.
#' @param sight.beta The vector of estimated coefficients for the logistic regression sightability model.
#' @param sight.beta.cov The covariance matrix of \code{sight.beta}
#' @param sight.logCI Should confidence intervals for the sightability adjusted estimates be computed
#'    using a normal-based confidence interval on \code{log(abundance)}
#' @param sight.var.method What method should be used to estimate the variances after adjusting
#'    for sightability.
#'    
