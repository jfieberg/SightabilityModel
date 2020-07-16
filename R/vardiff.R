#' Function to estimate the variance of the difference between two population
#' estimates
#' 
#' Function to estimate the variance of the difference between two population
#' estimates formed using the same sightability model (to correct for
#' detection).
#' 
#' Population estimates constructed using the same sightability model will NOT
#' be independent (they will typically exhibit positive covariance).  This
#' function estimates the covariance due to using the same sightability model
#' and subtracts it from the summed variance.
#' 
#' @param sight1 Sightability model object for the first population estimate
#' (formed by calling Sight.Est function)
#' @param sight2 Sightability model object for the second population estimate
#' (formed by calling Sight.Est function)
#' @return \item{vardiff}{numeric =
#' var(tau^[1])+var(tau^[2])-2*cov(tau^[1],tau^[2])}
#' @author John Fieberg
#' @keywords methods
#' @examples
#' 
#' 
#' # Example using moose survey data 
#'   data(obs.m) # observational moose survey data
#'   data(exp.m) # experimental moose survey data
#'   data(sampinfo.m) # information on sampling rates
#'  
#' # Estimate population size in 2006 and 2007 
#'  sampinfo <- sampinfo.m[sampinfo.m$year == 2007, ]
#'  tau.2007 <- Sight.Est(observed ~ voc, odat = obs.m[obs.m$year == 2007, ], 
#'                          sdat = exp.m, sampinfo.m[sampinfo.m$year == 2007, ], 
#'                          method = "Wong", logCI = TRUE, alpha = 0.05, Vm.boot = FALSE) 
#'  tau.2006 <- Sight.Est(observed ~ voc, odat = obs.m[obs.m$year == 2006, ],
#'                          sdat = exp.m, sampinfo.m[sampinfo.m$year == 2006, ],
#'                          method = "Wong", logCI = TRUE, alpha = 0.05, Vm.boot = FALSE) 
#' 
#' # naive variance
#'   tau.2007$est[2]+tau.2006$est[2]
#' 
#' # variance after subtracting positvie covariance
#'   vardiff(tau.2007, tau.2006)
#'
#' @export vardiff 
vardiff <-
function(sight1, sight2){
    if(sight1$call$form != sight2$call$form){
      stop("Need same sightability model to calculate the covariance")
    }  
    if(is.null(sight1$sight$note)){
      varbeta <- vcov(sight1$sight)
      beta <- coef(sight1$sight)
    }else{
      varbeta <- sight1$sight$varbet
      beta <- sight1$sight$bet
    }  
    y1 <- sight1$odat$total
    y2 <- sight2$odat$total
    n1 <- nrow(sight1$odat)
    n2 <- nrow(sight2$odat)
    inv.srate1 <- 1/sight1$odat$samp.rates
    inv.srate2 <- 1/sight2$odat$samp.rates
    fo <- sight1$call$form
    class(fo) <- "formula"
    tempnm1 <- terms(fo, data = sight1$odat)
    tempnm2 <- attr(tempnm1, "term.labels")
    covars1 <- sight1$odat[, tempnm2]  
    covars2 <- sight2$odat[, tempnm2]  

  # Do as much of the matrix multiplication outside of loop as possible
    xdat1 <- as.matrix(cbind(rep(1, n1), covars1))
    xdat2 <- as.matrix(cbind(rep(1, n2), covars2))
    xb1 <- xdat1%*%beta  # X'beta
    xb2 <- xdat2%*%beta  # X'beta
    xbb <- kronecker(xb1, t(xb2), FUN = "+")  # X1+x2
    xvarbeta <- xdat1%*%varbeta%*%t(xdat2)  # X Sig X
    smat <- matrix(0, n1, n2) # holder for some of the terms in the expression for smat
    for(i in 1:n1){
      for(j in 1:n2){
        xtemp1 <- as.vector(xdat1[i, ], mode = "numeric")
        xtemp2 <- as.vector(xdat2[j, ], mode = "numeric")
        xtot <- t(xtemp1+xtemp2)
        smat[i, j] <- (xtot%*%varbeta%*%t(xtot))/2 
      }
    }    
  
    smat.cov <- exp(-xbb-smat)*(exp(xvarbeta)-1)
  
    y.p1 <- as.matrix(y1*inv.srate1, n1, 1)
    y.p2 <- as.matrix(y2*inv.srate2, n2, 1)
 
    cov.total <- 2*t(y.p1)%*%smat.cov%*%y.p2
    vardiff <- sight1$est[2]+sight2$est[2]-cov.total
    colnames(vardiff) <- "Var(difference)"
    return(vardiff = vardiff)
}
