#' R function that interfaces with the SightabilityModel package and gives similar functionality
#' as the AerialSurvey program
#' 
#'  A stratified random sample of blocks in a survey area is conducted.
#'  In each block, groups of moose are observed (usually through an aerial survey).
#'  For each group of moose, the number of moose is recorded along with attributes
#'  such as sex or age.
#'  
#'  The SightabilityPopR() function  adjusts for sightability < 100\%.
#'
#' @template survey.data.input
#' @template sightability.input
#' @param conf.level Confidence level used to create confidence intervals.
#' @return  A data frame containing for each stratum and for all strata (identified as stratum id \code{.OVERALL}), the density,
#'    or abundance or ratio estimate along with its estimated standard error and large-sample normal-based confidence interval.
#'    Additional information on the components of variance is also reported.
#' @template author 
#' @template references
#' @keywords ~AerialSurvey ~MoosePop
#' @examples
#'  
#' ##---- See the vignettes for examples on how to run this analysis.
#' 
#' @export SightabilityPopR



SightabilityPopR <- function(
      survey.data,
      survey.block.area,
      stratum.data
      
      ,density=NULL           # variable for estimation as right-sided formulat
      ,abundance=NULL
      ,numerator=NULL, denominator=NULL
      
      ,sight.formula = observed ~ 1    # sightability function
      ,sight.beta = 10
      ,sight.beta.cov = matrix(0, nrow=1, ncol=1)
      ,sight.logCI=TRUE
      ,sight.var.method  = c("Wong","SS")[1]
      
      ,block.id.var  ="Block.ID" # variable identifying the block   id in the survey data and block area data
      ,block.area.var="Block.Area" # variable identifying the block area

      ,stratum.var   ="Stratum"   # variable identifying the stratum id in the survey.data and stratum.data    
      ,stratum.blocks.var="Stratum.Blocks"  # total number of blocks in the stratum      
      ,stratum.area.var  ="Stratum.Area"    # total area in the stratum      
      
      ,conf.level=0.90       # confidence interval level
){

  version <- "2021-01-01"
 
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
  if(length(common.names)>1)warning("Multiple common variables in survey.block.area and survey.data: ", paste(common.names, collapse=", "),
                                    "; Will merge on the ", block.id.var," variable only.", immediate.=TRUE)
  survey.data <- merge(survey.data, survey.block.area[,c(block.id.var, block.area.var)], by=block.id.var, all.x=TRUE)
  
# Merge the stratum information to the survey data
  common.names <- intersect(names(stratum.data), names(survey.data))
  if(length(common.names)>1)warning("Multiple common variables in stratum info and survey data: ", paste(common.names, collapse=", "),
                                    "; Will merge on the ", stratum.var," variable only.", immediate.=TRUE)
  survey.data <- merge(survey.data, stratum.data[,c(stratum.var, stratum.blocks.var, stratum.area.var)],      by=stratum.var)
  
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
    density = rhs.vars(density)
    if(length(density)>1)stop("Only one variable for density formula")
    if(!density %in% names(survey.data))stop("Density variable not in survey data for ",density)
    if(!is.numeric(survey.data[, density]))stop("Density variable in survey data not numeric for ", density)
    if(any(is.na(survey.data[,density])))stop("Missing data not allowed in survey data data values for ", density)
  }
  if(Type=="A"){
    abundance = rhs.vars(abundance)
    if(length(abundance)>1)stop("Only one variable for abundance formula")
    if(!abundance %in% names(survey.data))stop("abundance variable not in survey data for ",abundance)
    if(!is.numeric(survey.data[, abundance]))stop("abundance variable in survey data not numeric for ", abundance)
    if(any(is.na(survey.data[,abundance])))stop("Missing data not allowed in survey data data values for ", abundance)
  }
  if(Type=="R"){
    numerator = rhs.vars(numerator)
    if(length(numerator)>1)stop("Only one variable for numerator formula")
    if(!numerator %in% names(survey.data))stop("numerator variable not in survey data for ",numerator)
    if(!is.numeric(survey.data[, numerator]))stop("numerator variable in survey data not numeric for ", numerator)
    if(any(is.na(survey.data[,numerator])))stop("Missing data not allowed in survey data data values for ", numerator)
    denominator = rhs.vars(denominator)
    if(length(denominator)>1)stop("Only one variable for denominator formula")
    if(!denominator %in% names(survey.data))stop("denominator variable not in survey data for ",denominator)
    if(!is.numeric(survey.data[, denominator]))stop("denominator variable in survey data not numeric for ", denominator)
    if(any(is.na(survey.data[,denominator])))stop("Missing data not allowed in survey data data values for ", denominator)
  }

# make sure confidence level is sensible.
  if(conf.level < 0.5 | conf.level > .999)stop("Confidence level not sensible: ", conf.level)

# make sure that variables in sightability model are present on the survey.data dataframe
  #browser()
  if(is.null(sight.formula))stop("A sightability formula must be specified")
  if(!plyr::is.formula(sight.formula))stop("A sightability model must be a valid formula")
  xvars <- formula.tools::rhs.vars(sight.formula)
  if( !all(xvars %in% names(survey.data)))stop("Sightability model uses variables not in survey data")
  # check if any missing values
  if( any(is.na(survey.data[,xvars])))stop("Sightability model does not allow missing values in variables in sightability model")

# make sure beta coefficients are specified and same length as model
  if(is.null(sight.beta))stop("You must specify beta coefficients for sightability model")
  if(!is.numeric(sight.beta))stop("Beta coefficients must be numeric")
  if(length(sight.beta) != (1+length(xvars)))stop("Length of beta coefficients doesn't match sightability model")

# make sure that beta.vcov is right size
  if(is.null(sight.beta.cov))stop("You must specify beta covariance matrix for sightability model")
  if(!is.numeric(sight.beta.cov))stop("Beta covariance must be numeric")
  if(!is.matrix(sight.beta.cov))stop("Beta covariance must be a matrix")
  if(nrow(sight.beta.cov) != length(sight.beta) | ncol(sight.beta.cov) != length(sight.beta))
      stop("Beta covariance must be square and same size as beta coefficients")
 
# check that sight.logCI is logical
  if(!is.logical(sight.logCI) & !is.na(sight.logCI))stop("Illegal value for sight.logCI")
  
# make sure that sight.var method is one of Wong or SS
  if(!sight.var.method[1] %in% c("Wong","SS"))stop("Sightability variance method must be Wong or SS")

  
# Finally, now it is time to do the actual estimation of population abundance/densit corrected
# for sightability
  Var1 <- c(abundance, density, numerator)
  Var2 <- denominator
  if(is.null(Var2))Var2 <- NA

  # we must add variables "stratum", "subunit", to the data frame for passing to Sight.est and Sight.est.Ratio functions
  survey.data$stratum <- survey.data[, stratum.var ]
  survey.data$subunit <- survey.data[, block.id.var]
     
  # we must update the stratum information table
  # Number of observed blocks by stratum
  measured.blocks <- plyr::ddply(survey.data, stratum.var, function(x){
    nh=length(unique(x[,block.id.var]))
     data.frame(nh=nh)
  })
  stratum.data <- merge(stratum.data, measured.blocks, by=stratum.var)
  stratum.data$stratum <- stratum.data[,stratum.var]
  stratum.data$Nh      <- stratum.data[,stratum.blocks.var]
  #browser() 
  
  # for both density and abundance we find the total abundance using the SightabilityModel
  if(Type %in% c("D","A")){  # density or abundance estimator. 
    # we first get the (corrected) abundance estimator for each stratum
    # add the total variable to the data frame
    survey.data$total   <- survey.data[, Var1  ]
    stratum.res <- plyr::ddply(survey.data, stratum.var, function(x, stratum.data){
      # get the stratum data for this stratum
      stratum.data <- stratum.data[ stratum.data$stratum==x$stratum[1],,drop=FALSE]
      # drop observations with 0 animals measured in x
      x <- x[ x$total>0,]
      est.total <- SightabilityModel::Sight.Est(
                       form    = sight.formula, #sightability functional form
                       odat    = x,                # observed data
                       sampinfo= stratum.data,  # stratum information
                       bet     = sight.beta,
                       varbet  = sight.beta.cov,
                       logCI   = sight.logCI,
                       method  = sight.var.method)
      #browser()
      res.df <- data.frame(
                   Var1       = Var1,
                   Var2       = NA,
                   Var1.obs.total = sum(x[,Var1]),
                   Var2.obs.total = NA,
                   estimate   = est.total$est["tau.hat"],
                   SE         = sqrt(est.total$est["VarTot"]),
                   conf.level = conf.level
                 ) 
      res.df$LCL  = res.df$estimate + qnorm((1-conf.level)/2)    *res.df$SE
      res.df$UCL  = res.df$estimate + qnorm( 1-(1-conf.level)/2) *res.df$SE
 
      #browser()
      res.df <- cbind(res.df, data.frame(as.list(est.total$est)))
    }, stratum.data=stratum.data)
    #browser()
    # Now get the overall abundance
    survey.data <- survey.data[ survey.data$total>0,] # drop records with 0 animals observed
    est.total <- SightabilityModel::Sight.Est(
                       form    = sight.formula, #sightability functional form
                       odat    = survey.data,  # observed data
                       sampinfo= stratum.data,  # stratum information
                       bet     = sight.beta,
                       varbet  = sight.beta.cov,
                       logCI   = sight.logCI,
                       method  = sight.var.method)
    #browser()
    total.df <- data.frame(
           Var1.obs.total = sum(survey.data[,Var1]),
           Var2.obs.total = NA,
           estimate=est.total$est["tau.hat"],
           SE         =sqrt(est.total$est["VarTot"]),
           conf.level=conf.level)
    total.df$LCL <- total.df$estimate + qnorm((1-conf.level)/2)    *total.df$SE
    total.df$UCL <- total.df$estimate + qnorm( 1-(1-conf.level)/2) *total.df$SE
    total.df[,stratum.var] <- ".OVERALL"
    total.df <- cbind(total.df, data.frame(as.list(est.total$est)))
    stratum.res <- plyr::rbind.fill(stratum.res, total.df)
  }
     
  # abundance estimation is similar except we need to expand by stratum area
  if(Type =="D"){  # density estimation. Divide abundance estimates by the stratum total area
    #browser()
    stratum.res <- merge(stratum.res, stratum.data[,c(stratum.var, stratum.area.var)], all.x=TRUE, sort=FALSE)
    stratum.res[stratum.res[,stratum.var]==".OVERALL",stratum.area.var] <- sum(stratum.res[, stratum.area.var], na.rm=TRUE)
    stratum.res$estimate <- stratum.res$estimate / stratum.res[, stratum.area.var]
    stratum.res$SE       <- stratum.res$SE       / stratum.res[, stratum.area.var]
    stratum.res$LCL      <- stratum.res$LCL      / stratum.res[, stratum.area.var]
    stratum.res$UCL      <- stratum.res$UCL      / stratum.res[, stratum.area.var]
  }
     
  # ratio estimator. 
  # and then find the ratio and it se using a Taylor series
  if(Type %in% c("R")){  # ratio estimator
    # get the stratum specific estimates
    # add the numerator and denominator variables to the survey data
    survey.data$numerator   <- survey.data[,numerator]
    survey.data$denominator <- survey.data[,denominator]
    stratum.res <- plyr::ddply(survey.data, stratum.var, function(x, stratum.data){
      # get the stratum data for this stratum
      stratum.data <- stratum.data[ stratum.data$stratum==x$stratum[1],,drop=FALSE]
      est.ratio <- Sight.Est.Ratio(
                       form    = sight.formula, #sightability functional form
                       odat    = x,                # observed data
                       sampinfo= stratum.data,  # stratum information
                       bet     = sight.beta,
                       varbet  = sight.beta.cov,
                       logCI   = sight.logCI,
                       method  = sight.var.method)
      #browser()
      res.df <- data.frame(
                   Var1           = numerator,
                   Var2           = denominator,
                   Var1.obs.total = sum(x[,numerator]),
                   Var2.obs.total = sum(x[,denominator]),
                   estimate   = est.ratio$est["ratio.hat"],
                   SE         = sqrt(est.ratio$est["VarRatio"]),
                   conf.level = conf.level
                 ) 
      res.df$LCL  = res.df$estimate + qnorm((1-conf.level)/2)    *res.df$SE
      res.df$UCL  = res.df$estimate + qnorm( 1-(1-conf.level)/2) *res.df$SE
 
      #browser()
      res.df <- cbind(res.df, data.frame(as.list(est.ratio$est)))
    }, stratum.data=stratum.data)
    #browser()
    # Now get the overall ratio
    est.ratio <- Sight.Est.Ratio(
                       form    = sight.formula, #sightability functional form
                       odat    = survey.data,  # observed data
                       sampinfo=stratum.data,  # stratum information
                       bet     = sight.beta,
                       varbet  = sight.beta.cov,
                       logCI   = sight.logCI,
                       method  = sight.var.method)
    #browser()
    total.df <- data.frame(
          Var1.obs.total = sum(survey.data[,numerator]),
          Var2.obs.total = sum(survey.data[,denominator]),
          estimate=est.ratio$est["ratio.hat"],
          SE         =sqrt(est.ratio$est["VarRatio"]),
          conf.level =conf.level)
    total.df$LCL <- total.df$estimate + qnorm((1-conf.level)/2)    *total.df$SE
    total.df$UCL <- total.df$estimate + qnorm( 1-(1-conf.level)/2) *total.df$SE
    total.df[,stratum.var] <- ".OVERALL"
    total.df <- cbind(total.df, data.frame(as.list(est.ratio$est)))
    stratum.res <- plyr::rbind.fill(stratum.res, total.df)        
  } # end of ratio estimator
  stratum.res
}

