# summary method for ratio estimate from sightability model
#' @export 
summary.sightest_ratio <-
function(object,...){
  # ratio of numerator and denominator
    Tot.seen.numerator   <- sum(object$odat$numerator) 
    Tot.seen.denominator <- sum(object$odat$denominator)
  
  # z statistic
    z <- qnorm(1-object$alpha/2)
  
  # Normal assumption confidence interval  
    temp <- rep(object$est[1], 3)+c(0, -z, z)*sqrt(object$est[2])

    names(temp) <- NULL
    temp <- format(round(temp, 3), big.mark=",")
    cat("\n")
    temp2 <- paste("ratio.hat = ", temp[1], ";  ",100*(1-object$alpha), "% CI = (", temp[2], ", ", temp[3], ")")
    print(format(temp2[1], big.mark=","), quote = FALSE)
    out.summary <- list(ratio.hat = temp[1], lcl = temp[2], ucl = temp[3])  
}
