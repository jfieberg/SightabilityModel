#' R function that gives the same functionality as the MoosePop program.
#' 
#'  A stratified random sample of blocks in a survey area is conducted.
#'  In each block, groups of moose are observed (usually through an aerial survey).
#'  For each group of moose, the number of moose is recorded along with attributes
#'  such as sex or age.
#'  MoosePopR() assumes that sightability is 100%.
#'  Use the SightabilityPopR() function to adjust for sightability < 100%
#'
#' @template survey.data.input
#' @param conf.level Confidence level used to create confidence intervals.
#' @param survey.lonely.psu How to deal with lonely PSU within strata. See \code{surveyoptions} in the \code{survey} package. 
#' @return  A data frame containing for each stratum and for all strata (identified as stratum id \code{.OVERALL}), the density,
#'    or abundance or ratio estimate along with its estimated standard error and large-sample normal-based confidence
#' @template author 
#' @template references
#' @keywords ~MOOSEPOP ~moose surveys
#' @import formula.tools
#' @import survey
#' @importFrom Matrix bdiag 
#' @examples
#'  
#' ##---- See the vignettes for examples on how to run this analysis.
#' 
####' @import formula.tools
#' @export MoosePopR


MoosePopR <- function(
      survey.data,
      survey.block.area,
      stratum.data
      
      ,density=NULL       # variable to compute density as a formula
      ,abundance=NULL     # variable to compute abundance as a formula
      ,numerator=NULL     # numerator of ratio as a formula
      ,denominator=NULL  # denominator of ratio as a formula
      
      
      ,block.id.var  ="Block.ID" # variable identifying the block   id in the survey data and block area data
      ,block.area.var="Block.Area" # variable identifying the block area

      ,stratum.var       ="Stratum"   # variable identifying the stratum id in the survey.data and stratum.data    
      ,stratum.blocks.var="Stratum.Blocks"  # total number of blocks in the stratum      
      ,stratum.area.var  ="Stratum.Area"    # total area in the stratum      
      
      ,conf.level=0.90       # confidence interval level
      ,survey.lonely.psu="fail"  # default option for lonely psu
){

  options.old <- options()
  options(survey.lonely.psu=survey.lonely.psu)  # how to deal with lonely PSU
  version <- "2020-06-01"
# Error checking
# Be sure that survey.data, survey.block.area, and stratum.data are data.frames
  if( !is.data.frame(survey.data))      stop("survey.data is not a data frame")
  if( !is.data.frame(survey.block.area))stop("survey.block.area is not a data frame")
  if( !is.data.frame(stratum.data))     stop("stratum.data is not a data frame")

# Make sure that important variables are present
  if( is.null(stratum.var) || !is.character(stratum.var) || !length(stratum.var)==1)stop("stratum.var is missing or not a character or not length 1")
  if( is.null(block.id.var)|| !is.character(block.id.var)|| !length(stratum.var)==1)stop("stratum.var is missing or not a character or not length 1")
  if( is.null(stratum.blocks.var)|| !is.character(stratum.blocks.var)|| !length(stratum.blocks.var)==1)
       stop("stratum.blocks.var is missing or not a character or not length 1")
  if( is.null(stratum.area.var)|| !is.character(stratum.area.var)|| !length(stratum.area.var)==1)
       stop("stratum.area.var is missing or not a character or not length 1")
  if( is.null(block.area.var)|| !is.character(block.area.var)|| !length(block.area.var)==1)
       stop("block.area.var is missing or not a character or not length 1")

# Make sure that the stratum.var is present in the survey.data, stratum.data
  if( !stratum.var %in% names(survey.data)) stop("Stratum variable: ", stratum.var," not in survey data")
  if( !stratum.var %in% names(stratum.data))stop("Stratum variable: ", stratum.var," not in stratum data") 
   
# Make sure that stratum number of blocks and stratum area are in the stratum.data df
  if( !stratum.blocks.var %in% names(stratum.data))stop("Stratum variable for number of blocks not in stratum df: ", stratum.blocks.var)
  if( !stratum.area.var   %in% names(stratum.data))stop("Stratum variable for stratum area     not in stratum df: ", stratum.area.var)
# Make sure that these variables are numeric
  if( !is.numeric(stratum.data[,stratum.blocks.var,drop=TRUE]))stop("Stratum number of blocks not numeric")
  if( !is.numeric(stratum.data[,stratum.area.var  ,drop=TRUE]))stop("Stratum area not numeric")
   
# Convert stratum to character in all data frames as needed
  survey.data [, stratum.var] <- as.character(survey.data [, stratum.var, drop=TRUE])
  stratum.data[, stratum.var] <- as.character(stratum.data[, stratum.var, drop=TRUE]) 
# Make sure that stratum variable is not missing anywhere
  if(sum(is.na(survey.data [,stratum.var]))>0)stop("Stratum value missing in survey data")
  if(sum(is.na(stratum.data[,stratum.var]))>0)stop("Stratum value missing in stratum data")
  
# Ensure that stratum var matches across survey data and stratum data
  temp <- setdiff(stratum.data[,stratum.var,drop=TRUE], survey.data[,stratum.var,drop=TRUE])
  if( length(temp)>0)stop("Mis match in stratum ids - check value of : ", temp)
  temp <- setdiff(survey.data[,stratum.var,drop=TRUE], stratum.data[,stratum.var, drop=TRUE])
  if( length(temp)>0)stop("Mis match in stratum ids - check value of : ", temp)

# Make sure that the block.id.var is present in the survey.data, block area data
  if( !block.id.var %in% names(survey.data))      stop("Block.ID variable: ", block.id.var," not in survey data")
  if( !block.id.var %in% names(survey.block.area))stop("Block.ID variable: ", block.id.var," not in survey block area data")  
# Convert block ids to character. Ensure that match between area and survey.data frames
  survey.data      [, block.id.var ] <- as.character( survey.data      [, block.id.var, drop=TRUE])
  survey.block.area[, block.id.var ] <- as.character( survey.block.area[, block.id.var, drop=TRUE])
# Make sure that all survey blocks have an area
  temp <- setdiff(as.vector(survey.data[, block.id.var,drop=TRUE]), as.vector(survey.block.area[, block.id.var,drop=TRUE]))
  if(length(temp)>0)stop("Not all blocks in survey data have an area: ", paste(temp, collapse=" "))
# check for missing values in block id
  if( sum(is.na(survey.data[,block.id.var]))>0)stop("Missing value in block id in survey data")
# check that block areas are numeric and all present
  if( !is.numeric(survey.block.area[,block.area.var,drop=TRUE]))stop("block area must be numeric")
  if(sum(is.na(survey.block.area[,block.area.var,drop=TRUE]))>0)stop("Missing some block areas")

# check the block area to the survey data
# make sure no conflict with variable between the two sources except for block.id.var
  common.names <- intersect(names(survey.block.area), names(survey.data))
  if(length(common.names)>1)stop("Too many common variables in survey.block.area and survey.data: ", paste(common.names, collapse=" "))
  survey.data <- merge(survey.data, survey.block.area, by=block.id.var, all.x=TRUE)
  
# Merge the stratum information to the survey data
  common.names <- intersect(names(stratum.data), names(survey.data))
  if(length(common.names)>1)stop("Too many common variables in stratum info and survey data: ", paste(common.names, collapse=" "))
  survey.data <- merge(survey.data, stratum.data,      by=stratum.var)
  
# Check the density, abundance, numerator, and denominator variables
  Type= NA
  if(!is.null(density))    Type="D"
  if(!is.null(abundance))  Type="A"
  if(!is.null(numerator))  Type="R"
  if(!is.null(denominator))Type="R"
  if(is.na(Type))stop("Must specify one of density, abundance, numerator/denominator")
  
  if(Type=="D" && (!is.null(abundance) | !is.null(numerator) | !is.null(denominator)))
       stop("Only one density, abundance or numerator/denominator pair should be specified")
  if(Type=="A" && (!is.null(density)   | !is.null(numerator) | !is.null(denominator)))
       stop("Only one density, abundance or numerator/denominator pair should be specified")
  if(Type=="R" && (is.null(numerator)  | is.null(denominator)))stop("Both numerator and denominator must be specified")
  if(Type=="R" && (!is.null(density)   | !is.null(abundance)))
       stop("Only one density, abundance or numerator/denominator pair should be specified")
  
  if(Type=="D" && !formula.tools::is.one.sided(density))stop("Density specification must be a right-sided formula")
  if(Type=="A" && !formula.tools::is.one.sided(abundance))stop("Abundance specification must be a right-sided formulat")
  if(Type=="R" && (!formula.tools::is.one.sided(numerator) | !formula.tools::is.one.sided(denominator)))
     stop("Abundance specification must be a right-sided formulat")
     
  # extract the variables from the formula and check that valid
  if(Type=="D"){
    density = formula.tools::rhs.vars(density)
    if(length(density)>1)stop("Only one variable for density formula")
    if(!density %in% names(survey.data))stop("Density variable not in survey data for ",density)
    if(!is.numeric(survey.data[, density]))stop("Density variable in survey data not numeric for ", density)
    if(any(is.na(survey.data[,density])))stop("Missing data not allowed in survey data data values for ", density)
  }
  if(Type=="A"){
    abundance = formula.tools::rhs.vars(abundance)
    if(length(abundance)>1)stop("Only one variable for abundance formula")
    if(!abundance %in% names(survey.data))stop("abundance variable not in survey data for ",abundance)
    if(!is.numeric(survey.data[, abundance]))stop("abundance variable in survey data not numeric for ", abundance)
    if(any(is.na(survey.data[,abundance])))stop("Missing data not allowed in survey data data values for ", abundance)
  }
  if(Type=="R"){
    numerator = formula.tools::rhs.vars(numerator)
    if(length(numerator)>1)stop("Only one variable for numerator formula")
    if(!numerator %in% names(survey.data))stop("numerator variable not in survey data for ",numerator)
    if(!is.numeric(survey.data[, numerator]))stop("numerator variable in survey data not numeric for ", numerator)
    if(any(is.na(survey.data[,numerator])))stop("Missing data not allowed in survey data data values for ", numerator)
    denominator = formula.tools::rhs.vars(denominator)
    if(length(denominator)>1)stop("Only one variable for denominator formula")
    if(!denominator %in% names(survey.data))stop("denominator variable not in survey data for ",denominator)
    if(!is.numeric(survey.data[, denominator]))stop("denominator variable in survey data not numeric for ", denominator)
    if(any(is.na(survey.data[,denominator])))stop("Missing data not allowed in survey data data values for ", denominator)
  }
  
# make sure confidence level is sensible.
  if(conf.level < 0.5 | conf.level > .999)stop("Confidence level not sensible: ", conf.level)
  
# Finally, now it is time to do the actual estimation. 
  Var1 <- c(abundance, density, numerator)
  Var2 <- denominator
  if(is.null(Var2))Var2 <- NA
# Get the block total
  block.totals <- plyr::ddply(survey.data, c(stratum.var,block.id.var), function(x){
                              Var1.total = sum(x[,Var1])
                              Var2.total = ifelse(is.na(Var2),NA, sum(x[,Var2]))
                              data.frame(Var1.total  =Var1.total,
                                         Var2.total  =Var2.total)
  })
  # merge in block area and total stratum area values
  block.totals <- merge(block.totals, survey.block.area, by=block.id.var, all.x=TRUE)
  block.totals <- merge(block.totals, stratum.data,      by=stratum.var,  all.x=TRUE)

  # now for the analysis depending on the type of analysis
  #browser()
  # for both density and abundance we find the ratio estimator using the total block areas as weighting factors
  if(Type %in% c("D","A")){  # density or abundance estimator. A separate ratio estimator is used to combine over strata
    # first the individual strata estimates
    stratum.res <- plyr::ddply(block.totals, stratum.var, function(x){
      if(nrow(x)>1){ # replicate data from this stratum
         moose.design=survey::svydesign(data=x,
            ids=~1,
            fpc=reformulate(stratum.blocks.var))
         ratio.val <- survey::svyratio(numerator=~Var1.total, 
                            denominator=reformulate(block.area.var), 
                            design=moose.design)
         ratio.var.ci <- confint(ratio.val, level=conf.level)
         #browser()
         res.df <- data.frame(
               Var1       = Var1,
               Var2       = block.area.var,
               Var1.obs.total = sum(x$Var1.total),
               Var2.obs.total = sum(x[, block.area.var]),
               estimate   = coef(ratio.val),
               SE         = survey::SE(ratio.val),
               conf.level = conf.level,
               LCL        = ratio.var.ci[1],
               UCL        = ratio.var.ci[2]
            )
      }
      if(nrow(x)==1){  # lonely PSU. Set SE and LCL to NA
        res.df <- data.frame(
               Var1       = Var1,
               Var2       = block.area.var,
               Var1.obs.total = sum(x$Var1.total),
               Var2.obs.total = sum(x[, block.area.var]),
               estimate   = sum(x$Var1.total)/sum(x[, block.area.var]),
               SE         = NA,
               conf.level = conf.level,
               LCL        = NA,
               UCL        = NA
            )
      }
      res.df
    })
    # Now get the overall density using a ratio estimator
    moose.design=survey::svydesign(data=block.totals,
            ids=~1,
            strata=reformulate(stratum.var),
            fpc=reformulate(stratum.blocks.var))
    sep.ratio.val <- survey::svyratio(numerator=~Var1.total, 
                             denominator=reformulate(block.area.var), 
                             design=moose.design,
                             separate=TRUE)
    stratum.EF <- setNames( stratum.data[,stratum.area.var, drop=TRUE], stratum.data[,stratum.var,drop=TRUE]) 
    # be sure that expansion factors match the order in which the ratios are computed
    sep.ratio.val2 <- predict(sep.ratio.val, total=(stratum.EF/sum(stratum.EF))[names(sep.ratio.val$ratios)])
    #browser()
    total.df <- data.frame(
                  Var1.obs.total = sum(block.totals$Var1.total),
                  Var2.obs.total = sum(block.totals[, block.area.var]),
                  estimate       = as.numeric(sep.ratio.val2$total),
                  SE             = as.numeric(sep.ratio.val2$se),
                  conf.level     = conf.level)
    total.df$LCL <- total.df$estimate + qnorm((1-conf.level)/2)    *total.df$SE
    total.df$UCL <- total.df$estimate + qnorm( 1-(1-conf.level)/2) *total.df$SE
    total.df[,stratum.var] <- ".OVERALL"
    stratum.res <- plyr::rbind.fill(stratum.res, total.df)
  }
     
  # abundance estimation is similar except we need to expand by stratum area
  if(Type =="A"){  # abundance estimator. Multiply density estimates by the stratum total area
    stratum.res <- merge(stratum.res, stratum.data[,c(stratum.var, stratum.area.var)], all.x=TRUE, sort=FALSE)
    stratum.res[stratum.res[,stratum.var]==".OVERALL",stratum.area.var] <- sum(stratum.res[, stratum.area.var], na.rm=TRUE)
    stratum.res$estimate <- stratum.res$estimate * stratum.res[, stratum.area.var]
    stratum.res$SE       <- stratum.res$SE       * stratum.res[, stratum.area.var]
    stratum.res$LCL      <- stratum.res$LCL      * stratum.res[, stratum.area.var]
    stratum.res$UCL      <- stratum.res$UCL      * stratum.res[, stratum.area.var]
  }
     
  # ratio estimator. We need to get the population totals for the numerator and denominator variables
  # and then find the ratio and it se using a Taylor series
  if(Type %in% c("R")){  # ratio estimator
    # first the individual strata estimates
    stratum.res <- plyr::ddply(block.totals, stratum.var, function(x){
      if(nrow(x)>1){ # replicate data from this stratum
        moose.design=survey::svydesign(data=x,
            ids=~1,
            fpc=reformulate(stratum.blocks.var))
        ratio.val <- survey::svyratio(numerator    =~Var1.total, 
                            denominator =~Var2.total, 
                            design      =moose.design)
        ratio.var.ci <- confint(ratio.val, level=conf.level)
        #browser()
        res.df <- data.frame(
                   Var1       = Var1,
                   Var2       = Var2,
                   Var1.obs.total = sum(x$Var1.total),
                   Var2.obs.total = sum(x$Var2.total),
                   estimate   = coef(ratio.val),
                   SE         = survey::SE(ratio.val),
                   conf.level = conf.level,
                   LCL        = ratio.var.ci[1],
                   UCL        = ratio.var.ci[2]
                )
      }
      if(nrow(x)==1){ # lonely PSU
                res.df <- data.frame(
                   Var1       = Var1,
                   Var2       = Var2,
                   Var1.obs.total = sum(x$Var1.total),
                   Var2.obs.total = sum(x$Var2.total),
                   estimate   = sum(x$Var1.total)/sum(x$Var2.total),
                   SE         = NA,
                   conf.level = conf.level,
                   LCL        = NA,
                   UCL        = NA
                )

      }
      res.df
    })
    # now we need to estimate the numerator and denominator over the entire study using a separate ratio estimator
    # for both variables
    # we need to account for the covariance among the two variables within each strata and the 
    # separate ratio estimator does not return the full variance covariance matrix.
    # So we estimate the density in each stratum for both variables to get the full vcov matrix
    stratum.cov <- plyr::dlply(block.totals, stratum.var, function(x){
      if(nrow(x)>1){
         moose.design=survey::svydesign(data=x,
               ids=~1,
               fpc=reformulate(stratum.blocks.var))
         ratio.est <- survey::svyratio(numerator=~Var1.total+Var2.total, 
                             denominator = reformulate(block.area.var), 
                             design=moose.design, covmat=TRUE)
      }
      if(nrow(x)==1){
         stop("Sorry, unable to deal with lonely PSU when finding a ratio estimate")
      }
      ratio.est
    })
    # extract the density of var1, var2 in stratum 1, then var1 var2 in stratum 2 etc
    all.est <- as.vector(t(plyr::laply(stratum.cov, function(x){ unlist(coef(x))})))
    all.cov <- Matrix::bdiag(plyr::llply(stratum.cov, function(x){ vcov(x)}))
    temp <- as.data.frame(stratum.data)
    row.names(temp) <- temp[, stratum.var, drop=TRUE]
    all.weight <- rbind(as.vector(t(cbind(temp[ names(stratum.cov)    , stratum.area.var,drop=TRUE], 0))),
                        as.vector(t(cbind(0, temp[ names(stratum.cov)    , stratum.area.var,drop=TRUE]))))
    pop.totals.est <- all.weight %*% all.est
    pop.totals.cov <- all.weight %*% all.cov %*% t(all.weight)
    pop.ratio = pop.totals.est[1] / pop.totals.est[2]
    pop.ratio.se = msm::deltamethod( ~x1/x2, pop.totals.est, pop.totals.cov, ses=TRUE)
    #browser()
    total.df <- data.frame(
                        Var1.obs.total = sum(block.totals$Var1.total),
                        Var2.obs.total = sum(block.totals$Var2.total),
                        estimate   = pop.ratio,
                        SE         = pop.ratio.se,
                        conf.level = conf.level)
    total.df$LCL <- total.df$estimate + qnorm((1-conf.level)/2)    *total.df$SE
    total.df$UCL <- total.df$estimate + qnorm( 1-(1-conf.level)/2) *total.df$SE
    total.df[,stratum.var] <- ".OVERALL"
        
    stratum.res <- plyr::rbind.fill(stratum.res, total.df)
  } # end of ratio estimator
  options(options.old)  # restore options
  stratum.res
}
