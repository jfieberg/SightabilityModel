#' Sightability Model Estimator
#' 
#' Estimates population abundance by 1) fitting a sightability (logistic
#' regression) model to "test trial" data; 2) applying the fitted model to
#' independent (operational) survey data to correct for detection rates < 1.
#' 
#' Variance estimation methods: method = Wong implements the variance estimator
#' from Wong (1996) and is the recommended approach.  Method = SS implements
#' the variance estimator of Steinhorst and Samuel (1989), with a modification
#' detailed in the Appendix of Samuel et al. (1992).
#' 
#' Estimates of the variance may be biased low when the number of test trials
#' used to estimate model parameters is small (see Wong 1996, Fieberg and
#' Giudice 2008).  A bootstrap can be used to aid the estimation process by
#' specifying Vm.boot = TRUE [note: this method is experimental, and can be
#' time intensive].
#' 
#' Confidence interval construction: often the sampling distribution of tau^ is
#' skewed right.  If logCI = TRUE, the confidence interval for tau^ will be
#' constructed under an assumption that (tau^ - T) has a lognormal distribution,
#' where T is the total number of animals seen.  In this case, the upper and
#' lower limits are constructed as follows [see Wong(1996, p. 64-67)]:
#' 
#' LCL = T + [(tau^-T)/C]*sqrt(1+cv^2), UCL = T+[(tau^-T)*C]*sqrt(1+cv^2),
#' where cv^2 = var(tau^)/(tau^-T)^2 and C = exp[z[alpha/2]*sqrt(ln(1+cv^2))].
#' 
#' @param form a symbolic description of the sightability model to be fit
#' (e.g., "y ~ x1 + x2 + ..."), where y is a binary response variable (= 1 if
#' the animal is seen and 0 otherwise) and x1, x2, ... are a set of predictor
#' variables thought to influence detection
#' @param sdat 'sightability' data frame.  Each row represents an independent
#' sightability trial, and columns contain the response (a binary random
#' variable = 1 if the animal was observed and 0 otherwise) and the covariates
#' used to model detection probabilities.
#' @param odat 'observational survey' data frame containing the following
#' variable names (\emph{stratum, subunit, total}) along with the same
#' covariates used to model detection probabilities (each record corresponds to
#' an independently sighted group of animals).  \emph{stratum} = stratum
#' identifier (will take on a single value for non-stratified surveys);
#' \emph{subunit} = numeric plot unit identifier; \emph{total} = total number
#' of observed animals (for each independently sighted group of animals).
#' @param sampinfo data frame containing sampling information pertaining to the
#' observational survey.  Must include the following variables (\emph{stratum,
#' nh, Nh}).  \emph{stratum} = stratum identifier (must take on the same values
#' as \emph{stratum} variable in observational data set), \emph{nh} = number of
#' sampled units in stratum h, \emph{Nh} = number of population units in
#' stratum h; note (this dataset will contain a single record for
#' non-stratified designs).
#' 
#' @param method method for estimating variance of the abundance estimator.
#' Should be one of ("Wong", "SS").  See details for more information.
#' @param logCI Boolean variable, default (= TRUE), indicates the confidence
#' interval should be constructed under the assumption that (tau^ - T) has a
#' lognormal distribution, where T is the total number of animals observed
#' (see details)
#' @param alpha type I error rate for confidence interval construction
#' @param Vm.boot Boolean variable, when = TRUE indicates a bootstrap should be
#' used to estimate cov(theta[i,j],theta[i',j']), var/cov matrix of the
#' expansion factors (1/detection prob)
#' @param nboot number of bootstrap replicates to use if Vm.boot = TRUE
#' @param bet regression parameters (if the sightability model is not to be fit
#' by Sight.Est).  Make sure the order is consistent with the specification in
#' the "form" argument.
#' @param varbet variance-covariance matrix for beta^ (if the sightability
#' model is not to be fit by Sight.Est).  Make sure the order is consistent
#' with the specification in the "form" argument.
#' @return An object of class \code{sightest}, a list that includes the
#' following elements: \item{sight.model}{the fitted sightability model}
#' \item{est}{abundance estimate [tau.hat] and its estimate of uncertainty
#' [Vartot] as well as variance components due to sampling [Varsamp], detection
#' [VarSight], and model uncertainty [VarMod]} The list also includes the
#' original test trial and operational survey data, sampling information, and
#' information needed to construct a confidence interval for the population
#' estimate.
#' @author John Fieberg, Wildlife Biometrician, Minnesota Department of Natural
#' Resources
#' @references
#' 
#' Fieberg, J.  2012.  Estimating Population Abundance Using Sightability
#' Models: R SightabilityModel Package. Journal of Statistical Software, 51(9),
#' 1-20.  URL https://doi.org/10.18637/jss.v051.i09.
#' 
#' Fieberg, John and Giudice, John. 2008 Variance of Stratified Survey
#' Estimators With Probability of Detection Adjustments. Journal of Wildlife
#' Management 72:837-844.
#' 
#' Samuel, Michael D. and Steinhorst, R. Kirk and Garton, Edward O. and
#' Unsworth, James W. 1992.  Estimation of Wildlife Population Ratios
#' Incorporating Survey Design and Visibility Bias.  Journal of Wildlife
#' Management 56:718-725.
#' 
#' Steinhorst, R. K., and M.D. Samuel. 1989.  Sightability adjustment methods
#' for aerial surveys of wildlife populations.  Biometrics 45:415-425.
#' 
#' Wong, C. 1996.  Population size estimation using the modified
#' Horvitz-Thompson estimator with estimated sighting probabilities.
#' Dissertation, Colorado State University, Fort Collins, USA.
#' @keywords survey models
#' @examples
#' 
#' # Load data frames
#'   data(obs.m) # observational survey data frame
#'   data(exp.m) # experimental survey data frame
#'   data(sampinfo.m) # information on sampling rates (contained in a data frame)
#'  
#' # Estimate population size in 2007 only
#'   sampinfo <- sampinfo.m[sampinfo.m$year == 2007,]
#'   Sight.Est(observed ~ voc, odat = obs.m[obs.m$year == 2007,],
#'     sdat = exp.m, sampinfo, method = "Wong", 
#'     logCI = TRUE, alpha = 0.05, Vm.boot = FALSE) 
#' 
#' 
#' # BELOW CODE IS SOMEWHAT TIME INTENSIVE (fits models using 2 variance estimators to 3 years of data)
#' # Estimate population size for 2004-2007
#' # Compare Wong's and Steinhorst and Samuel variance estimators
#'   tau.Wong <- tau.SS <- matrix(NA,4,3)
#'   count <- 1
#'   for(i in 2004:2007){
#'     sampinfo <- sampinfo.m[sampinfo.m$year == i,]
#' 
#' # Wong's variance estimator 
#'     temp <- Sight.Est(observed ~ voc, odat = obs.m[obs.m$year == i,],
#'        sdat = exp.m, sampinfo, method = "Wong", 
#'        logCI = TRUE, alpha = 0.05, Vm.boot = FALSE) 
#'     tau.Wong[count, ] <- unlist(summary(temp)) 
#'  
#' # Steinhorst and Samuel (with Samuel et al. 1992 modification) 
#'     temp <- Sight.Est(observed ~ voc, odat = obs.m[obs.m$year == i,],  
#'        sdat = exp.m, sampinfo, method = "SS")
#'     tau.SS[count, ] <- unlist(summary(temp)) 
#'     count<-count+1
#'   }  
#'   rownames(tau.Wong) <- rownames(tau.SS) <- 2004:2007
#'   colnames(tau.Wong) <- colnames(tau.SS) <- c("tau.hat","LCL","UCL")
#'   (tau.Wong <- apply(tau.Wong, 1:2, 
#'       FUN=function(x){as.numeric(gsub(",", "", x, fixed = TRUE))}))
#'   (tau.SS <-   (tau.Wong <- apply(tau.Wong, 1:2, 
#'      FUN = function(x){as.numeric(gsub(",", "", x, fixed = TRUE))})))
#' 
#' \dontrun{
#'   require(gplots)
#'   par(mfrow = c(1,1))
#'     plotCI(2004:2007-.1, tau.Wong[,1], ui = tau.Wong[,3], 
#'         li = tau.Wong[,2], type = "l", xlab = "", 
#'         ylab = "Population estimate", xaxt = "n",
#'         xlim=c(2003.8, 2007.2))
#'     plotCI(2004:2007+.1, tau.SS[,1], ui = tau.SS[,3], li = tau.SS[,2], 
#'          type = "b", lty = 2, add = TRUE)
#'     axis(side = 1, at = 2004:2007, labels = 2004:2007)  
#'   }
#'
#' @export Sight.Est
#' @import stats
#' @import utils


Sight.Est <-
function(form, sdat=NULL, odat, sampinfo, method="Wong", logCI=TRUE, alpha=0.05, Vm.boot = FALSE, nboot = 1000, bet = NULL, varbet = NULL){

  # form = model formula (for fitting logistic regression model to sightability dataset)
  # sdat = sightability dataset containing response and sightability covariates
  # odat = dataset containing observed groups, sample unit ids, stratum identifiers, and sampling rates
  # sampinfo = data frame with sampling information (number of sample (nh) and population units (Nh) in each stratum
  # Method = method of variance estimation 
  # logCI = should the confidence interval be constructed under the assumption that tau^ has a lognormal distribution
  # alpha = type I error rate for confidence interval construction
  # Vm.boot = should a bootstrap be used to estimate cov(theta[i,j],theta[i',j']), var/cov matrix of the expansion factors for detection
  # nboot = number of bootstraps if Vm.boot = TRUE
  # beta=regression parameters (can be passed, rather than using Sight.Est to fit the model)
  # varbet= var/cov matrix for regression parameter estimates (can be passed, rather than using Sight.Est to fit the model)

  # Check arguments
    if(is.null(sdat) & (is.null(bet) | is.null(varbet) )){
      stop("Need to specify sightability data set or beta^ and var(beta^)")
    }   
    if(method %in%c("Wong", "SS") != TRUE){
      stop("Method must be Wong or SS")
    }
    if(is.null(bet)){   
      sight.model <- glm(form, family = binomial(), data = sdat)
      beta <- coef(sight.model)
      varbet <- summary(sight.model)$cov.unscaled
    }else{
      sight.model <- NULL
      sight.model$bet <- bet 
      sight.model$varbet <- varbet
      sight.model$note <- "User supplied regression model"
      beta <- bet
      varbet <- varbet
    }  
    if(sum(c("stratum", "subunit", "total")%in%names(odat)) != 3){
      stop("Need to have variables stratum, subunit,  and total in observational survey dataset.  These are CASE-SENSITIVE")
    
    }  
    if(sum(odat$total == 0) > 0){
       print("Dropping records with 0 animals in the observational data set")
       odat <- odat[odat$total > 0, ]
    }   
    if(sum(c("nh","Nh", "stratum")%in%names(sampinfo)) != 3){
      stop("Need to have variables nh, Nh, and stratum in dataset containing sampling information")
    }  

 # Make sure all stratum names are in both observational and sampling information data sets
    s1 <- sort(unique(odat$stratum))
    s2 <- sort(unique(sampinfo$stratum))
    if(length(s1) != length(s2) | sum(s1%in%s2) != length(s2)){
      stop("The same list of stratum must be present in both both observational and sampling information data sets")
    }  
 # create sampling variables
    nh <- sampinfo$nh
    Nh <- sampinfo$Nh
  
 # Make sure strata are numbered 1:h
    sampinfo$stemp <- 1:nrow(sampinfo)
  
 # Want srates to be same legnth as total vector
    sampinfo$samp.rates <- sampinfo$nh/sampinfo$Nh
    odat <- merge(odat, sampinfo, by.x = "stratum", by.y = "stratum")  
  
 #pull off other vairables
    srates <- odat$samp.rates 
    stratum <- odat$stemp
    subunit <- odat$subunit
    total <- odat$total
   
    if(is.null(bet)){ 
      tempnm <- terms(form, data = sdat)
    }else{tempnm <- terms(form, data = odat)}
    tempnm2 <- attr(tempnm, "term.labels")
    if(sum(tempnm2%in%names(odat)) != length(tempnm2)){
       print("The exact same names need to be used for covariates in the sightability and operational survey datasets")
    }
    covars <- odat[, tempnm2]  
 
 #Bootstrap for Cov(theta)? 
    if(Vm.boot == TRUE){
      bets <- matrix(NA, nboot, length(beta))
      varbets <- array(NA,dim=c(nboot, length(beta), length(beta)))
      nsdat <- nrow(sdat)
      for(i in 1:nboot){
        bdat <- sdat[sample(1:nsdat, replace = TRUE), ] # sample observations with replacement (non-parametric bootstrap)
        fit <- glm(form, family = binomial, data = bdat)  # Fit model to bootstrap data
        bets[i,] <- coef(fit)  # pull of coefficients
        varbets[i,,] <- summary(fit)$cov.unscaled  # pull off variance/covariance matrix
      }
      smat <- covtheta(total, srates, stratum, subunit, covars, bets, varbets, nboot)
    }
    else{smat <- NULL}  
    if(method == "Wong"){
      est <- Wong.est(total, srates, nh, Nh, stratum, subunit, covars, beta, varbet, smat)
    }
    if(method == "SS"){
      est <- SS.est(total, srates, nh, Nh, stratum, subunit, covars, beta, varbet, smat)
    }    
    out <- NULL
    out$call <- match.call()
    out$sight.model <- sight.model
    out$est <- est$est
    out$var.method <- est$var.method 
    if(is.null(sdat) != TRUE){out$sdat <- sdat}
      out$odat <- odat
      out$samp <- sampinfo[,c("stratum", "nh", "Nh")]
      out$CI.method <- "lognormal"
      if(logCI == FALSE){  out$CI.method <- "normal"}
      out$alpha <- alpha
      class(out) <- "sightest"
      out
}
