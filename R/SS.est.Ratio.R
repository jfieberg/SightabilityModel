#' Sightability estimate or ratio with variance components estimator from
#' Steinhorst and Samuel (1989) and Samuel et al. (1992).  This is merely a
#' stub and has not been implemented.
#' 
#' Estimates ratio, with variance estimated using Steinhorst and Samuel (1989)
#' and Samuel et al.'s (1992) estimator.  Usually, this function will be
#' called by Sight.Est.Ratio()
#' 
#' 
#' @param numerator,denominator Number of animals for the numerator and
#' denominator of the ratio in each independently sighted group
#' @param srates Plot-level sampling probability
#' @param nh Number of sample plots in each stratum
#' @param Nh Number of population plots in each stratum
#' @param stratum Stratum identifiers (associated with the independently
#' observed animal groups)
#' @param subunit Plot ID (associated with the independently observed animal
#' groups)
#' @param covars Matrix of sightability covariates (associated with the
#' independently observed animal groups)
#' @param beta Logistic regression parameter estimates (from fitted
#' sightability model)
#' @param varbeta Estimated variance-covariance matrix for the logistic
#' regression parameter estimates (from fitted sightability model)
#' @param smat Estimated variance-covariance matrix for the inflation factors
#' (1/probability of detection).  This is an n.animal x n.animal matrix, and is
#' usually calculated within the SS.est.Ratio function.  Non-null values can be
#' passed to the function (e.g., if a bootstrap is used to estimate uncertainty
#' due to the estimated detection parameters).
#' @return \item{ratio.hat}{Sightability estimate of ratio, ratio^}
#' \item{VarRatio}{Estimated variance of ratio^} \item{VarSamp, VarSight,
#' VarMod}{Estimated variance component due to sampling, sightability and model
#' set to NA}
#' @author Carl James Schwarz, cschwarz.stat.sfu.ca@@gmail.com
#' @seealso \code{\link{Sight.Est}}, \code{\link{Wong.est}}
#' @references
#' 
#' Steinhorst, R. K., and M.D. Samuel. 1989.  Sightability adjustment methods
#' for aerial surveys of wildlife populations.  Biometrics 45:415-425.
#' 
#' Wong, C. 1996.  Population size estimation using the modified
#' Horvitz-Thompson estimator with estimated sighting probabilities.
#' Dissertation, Colorado State University, Fort Collins, USA.
#' @keywords methods
#' @export SS.est.Ratio

SS.est.Ratio <-
function(numerator, denominator, srates, nh, Nh, stratum, subunit, covars, beta, varbeta, smat=NULL){

  # estimates the ratio of the numerator/denominator and variance using SS results for the numerator
  # and denominator, and then uses Wang (1996) results to estimate the covariance
  # Questions to cschwarz.stat.sfu.ca@gmail.com
  
  #  Arguments (the first five should have one observation for each individual in the obs. dataset
  # numerator = total number of animals in the observed group for the numerator of the ratio
  # denominator = total number of animals in the observed group for the denominator of the ratio
  # srates = sampling rate for the stratum (for each observed animal)
  # nh = number of sampled units in each stratum
  # Nh = number of population units in each stratum
  # stratum = stratum identifier  
  # Subunit = subunit identifier 
  # covars = matrix of covariates used in the sightability model
  # beta = vector of parameter estimates from logistic model
  # varbeta = var/cov matrix for above parameter estimates
  # smat = variance covariance matrix of estimated inflation factors (theta)
  
  version <- "2020-06-01"
  
  # Check on length of vector arguments
    if(length(numerator) != length(denominator))stop("Numerator and Denominator must be equal length")
    n <- length(numerator)
    if(length(srates) != n) {stop("Srates vector needs to be same dimension as numerator and denomiantor vectors")}
    if(length(stratum) != n) {stop("Stratum vector needs to be same dimension as numerator and denominator vectors")}
    if(length(srates) != n) {stop("Subunit vector needs to be same dimension as numerator and denominator vectors")}

  # Form xmatrix
    xdata <- as.matrix(cbind(rep(1, n), covars))
    ncovars <- dim(xdata)
    if(ncovars[1] != n) {stop("Covariate matrix needs to be same length as numerator and denominator vectors")}  
      
  # Make sure beta is a matrix of dim ncovar x 1 
    beta <- matrix(beta, length(beta), 1)
    if (length(beta) != ncovars[2]){stop("Covars matrix must have same number of columns as the Beta vector")} 
 
  #  Estimate theta for each animal in each group = 1/prob(detection)
    theta <- 1+exp(-xdata%*%beta-diag(xdata%*%varbeta%*%t(xdata))/2)  

  # Xie = prob. detection = 1/ theta  
    xie <- 1/theta
  
  # Point estimate for the numerator and denomiator variables using the Wong.est routines
    numerator.hat   <- SS.est(numerator,   srates, nh, Nh, stratum, subunit, covars, beta, varbeta, smat=NULL)
    denominator.hat <- SS.est(denominator, srates, nh, Nh, stratum, subunit, covars, beta, varbeta, smat=NULL)
    ratio.hat <- as.numeric(numerator.hat$est["tau.hat"] /denominator.hat$est["tau.hat"]) # strip off attributes
    #browser()
 
  #------------------------- Now, compute the variance of the ratio  ----------------------------------
  #  This requires calculation of Cov(theta^[ij], theta^[i'j']) (see Wong 1996)
  #  So, start by calculating variance covariance matrix of theta values = smat (below)
    nxdata <- nrow(xdata)
    sm3 <- matrix(0, nxdata, nxdata) # holder for some of the terms in the expression for smat
   
   if(is.null(smat)){ #Variance/covariance matrix of smat not supplied by user
                    # Use asymptotic estimate based on log-normal assumption (of SS and Wong)  
  # Do as much of the matrix multiplication outside of loop as possible
    xb <- xdata%*%beta  # X'beta
    xbb <- kronecker(xb, t(xb), FUN = "+")  # x1+x2
    xvarbeta <- xdata%*%varbeta%*%t(xdata)  # X Sig X
    for(i in 1:nxdata){
      for(j in i:nxdata){
       xtemp1 <- as.vector(xdata[i, ], mode = "numeric")
       xtemp2 <- as.vector(xdata[j, ], mode = "numeric")
       xtot <- t(xtemp1+xtemp2)
       sm3[i, j] <- sm3[j, i] <- xtot%*%varbeta%*%t(xtot)/2
      }
    }
    smat <- exp(-xbb-sm3)*(exp(xvarbeta)-1)
  }
  
  #------------------------ End Calculation of Cov matrix ------------------------------------------------

  # Get various quantities read
  #  Calculate MKs and pks (MKs = corrected totals by subunit, pks = sampling rate for the subunit)
  #browser()
  # We need to adjust the numerator and denominator variables for sightability
  MKs <- rowsum(numerator*denominator*theta, subunit)
  
  pks <- tapply(srates,subunit, unique)
  nk <- length(pks)
 
  # Turn MKs and pks into matrices with 1 columns 
  MKs <- matrix(MKs, length(MKs), 1)
  pks <- matrix(pks, length(pks), 1)
    
  # get information at the subunit level
  first.subunit <- !duplicated(subunit)
  subunit.info <- data.frame(subunit      =subunit[first.subunit],
                             stratum      =stratum[first.subunit],
                             stratum.srate=srates[first.subunit]) 
  #---------------------- Calculate covariance of numerator and denominator total
  # Refer to equation 2.4.6 of Wong(1996)
  #
  # Term1 = produce of the totals
  E2.4.6.t1 <- as.numeric(numerator.hat$est["tau.hat"] * denominator.hat$est["tau.hat"])
  
  # Term 2 <- component due to sampling
  E2.4.6.t2 <- sum(MKs / pks)
  
  # term 3 - within subunit terms. We loop over the subunits
  # The stratum vector has values from 1:nstrata from the calling function
  E2.4.6.t3 <- 0
  for(i.subunit in 1:nrow(subunit.info)){
     i.subunit.val = subunit.info$subunit[i.subunit]
     i.index <- (1:length(subunit))[subunit==i.subunit.val] # which elements belong to this subunit
     cross.num.denom <- outer(numerator[i.index], denominator[i.index], FUN="*")
     cross.theta     <- outer(theta    [i.index], theta      [i.index], FUN="*")
     temp <- cross.num.denom*(-smat[i.index,i.index]+cross.theta)
     diag(temp) <- 0  # we only want the off diagonal terms
     E2.4.6.t3 <- E2.4.6.t3 + sum(temp)/subunit.info$stratum.srate[i.subunit] 
  }
  
  # term 4 - across subunit terms
  #browser()
  E2.4.6.t4 <- 0
  for(i.subunit in 1:nrow(subunit.info)){
      i.subunit.val = subunit.info$subunit[i.subunit]
      i.index <- (1:length(stratum))[subunit==i.subunit.val] # which elements belong to this stratum
      for(j.subunit in 1:nrow(subunit.info)){
        j.subunit.val = subunit.info$subunit[j.subunit]
        j.index <- (1:length(stratum))[subunit==j.subunit.val] # which elements belong to this stratum
        if(i.subunit != j.subunit){
          cross.num.denom <- outer(numerator[i.index], denominator[j.index], FUN="*")
          cross.theta     <- outer(theta    [i.index], theta      [j.index], FUN="*")
          temp <- cross.num.denom*(-smat[i.index,j.index]+cross.theta)
          pij <- ifelse(subunit.info$stratum[i.subunit]==subunit.info$stratum[j.subunit]
                         # subunits from same stratum are not indepent since sampling is without replacement
                        ,nh[subunit.info$stratum[i.subunit]]/Nh[subunit.info$stratum[i.subunit]]*
                        (nh[subunit.info$stratum[j.subunit]]-1)/(Nh[subunit.info$stratum[j.subunit]]-1)
                         # subunits from different strata are selected independently
                        ,nh[subunit.info$stratum[i.subunit]]/Nh[subunit.info$stratum[i.subunit]]* # different strata independent selection
                         nh[subunit.info$stratum[j.subunit]]/Nh[subunit.info$stratum[j.subunit]])
          E2.4.6.t4 <- E2.4.6.t4 + sum(temp)/pij
        }
      }
  }
    #browser()
  # Now, add then all toegher
  # I think that there is an error in Wong's thesis in 2.4.5 and it should read
    cov.nhat.dhat <- E2.4.6.t1 - (E2.4.6.t2 + E2.4.6.t3 + E2.4.6.t4)
  # rather than
  # cov.nhat.dhat <- E2.4.6.t1 + (E2.4.6.t2 + E2.4.6.t3 + E2.4.6.t4)
    # now add together the variance components
    # see Equation 2.4.4 of Wong(1996)
    VarRatio <- as.numeric(ratio.hat^2*(
                 numerator.hat  $est["VarTot"]/numerator.hat  $est["tau.hat"]^2 +
                 denominator.hat$est["VarTot"]/denominator.hat$est["tau.hat"]^2 -
                 2*cov.nhat.dhat / numerator.hat  $est["tau.hat"] /denominator.hat$est["tau.hat"])
                )

    # we don't separate out the variance components due to sightability and model but see the num/denom variables
    out <- c(ratio.hat=ratio.hat, 
             VarRatio=VarRatio,
             VarSamp = NA,
             VarSight= NA,
             VarMod  = NA)
    out2 <- NULL
    out2$est <- out
    out2$var.method = "SS"
    out2$numerator   <- numerator.hat
    out2$denominator <- denominator.hat
    out2$cov.nhat.dhat.components <- data.frame(E2.4.6.t1, E2.4.6.t2, E2.4.6.t3, E2.4.6.t4,cov.nhat.dhat )
    out2$version <- version
    out2
}
