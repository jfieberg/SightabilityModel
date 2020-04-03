# Print method for sightability estimates of a ratio
#' @export

print.sightest_ratio <-
function(x,...){
    cat("Call:\n")
    print(x$call)
    cat("\n------------------- SIGHTABILITY MODEL ---------------------\n")
    print(x$sight.model)
    cat("\n----------------- Population Survey data  ----------------\n")
    cat("\n Stratum Sampling Information\n")
    print(x$samp)
   
    cat("\n Number of animals seen in each stratum for numerator and denominator\n")
    ua.num  <- tapply(x$odat$numerator,   x$odat$stratum, sum)
    ua.denom<- tapply(x$odat$denominator, x$odat$stratum, sum)
    print(format(ua.num,   big.mark=","), quote = FALSE)
    print(format(ua.denom, big.mark=","), quote = FALSE)

        cat("\n-------------- Ratio ESTIMATE (",100*(1-x$alpha),"% CI) ----------------\n")
    z <- qnorm(1-x$alpha/2)
    temp <- rep(x$est[1], 3)+c(0, -z, z)*sqrt(x$est[2])
    names(temp) <- NULL
    temp <- format(round(temp, 3), big.mark=",")
    cat("\n")
    temp2 <- paste("ratio.hat = ", temp[1], ";  ",100*(1-x$alpha), "% CI = (", temp[2], ", ", temp[3], ")")
    print(temp2[1], quote = FALSE)
    cat("\n")
    cat("\n------------------  SE(ratio.hat) --------------------------------\n")
    cat("Variance method: "); print(x$var.method)
    SE <- sqrt(x$est[2])
    names(SE) <- "SE"
    print(format(round(SE, 3),  big.mark=","), quote = FALSE)
}
