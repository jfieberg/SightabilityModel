#' Generate a bootstrap replicate of data for call to MoosePopR_DomStrat()
#' 
#' This function takes the data from a classical/domain stratification and generates
#' a bootstrap replicate suitable for analysis using MoosePopR_DomStrat().
#' A sightability model is allowed which "adjusts" the input data for sightability.
#' This can also be used for SightabilityPopR() models by forcing block areas to 1 and
#' the total block area in stratum to the number of blocks to mimic a mean-per-unit estimator.
#' See the vignette for examples of usage.

#' @template stratum.data.input
#' @template domain.var
#' @template block.id.var
#' @template block.area.var
#' @template stratum.var

#' @template selected.unit.data.input
#' 
#' @template waypoint.data.input
#' @template den.abund.num.denom.var
#' 
#' @param sight.model A formula that identifies the model used
#'   to estimate sightability. For example \code{observed ~ VegCoverClass} would indicate
#'   that sightability is a function of the \code{VegCoverClass} variable in the survey 
#'   data. The left hand variable is arbitrary. The right hand variables must be present
#'   in the survey.data data frame.
#' @param sight.beta The vector of estimated coefficients for the logistic regression sightability model.
#' @param sight.beta.cov The covariance matrix of \code{sight.beta}

#' @template conf.level
#' @param survey.lonely.psu How to deal with lonely PSU within strata. See \code{surveyoptions} in the \code{survey} package.
#'
#' @param check.args Should arguments be checked. Turn off for extensive bootstrapping to save time.
#'   
#' @return  A list containing the input data (\code{input.data}), 
#' the bootstrap replicate (\code{boot.data}), and a data frame (\code{boot.res}) with the estimated density,
#'    or abundance or ratio along with its estimated standard error and large-sample normal-based confidence interval.
#'    The density/abundance/ratio over all strata is also given on the last line of the data.frame.
#' @template author 
#' @template references
#' @keywords ~MOOSEPOP ~moose surveys
#' @import plyr
#' @importFrom formula.tools rhs.vars
#' @importFrom mvtnorm rmvnorm
#' @examples
#'  
#' ##---- See the vignettes for examples on how to use this function
#' 
#' @export MoosePopR_DomStrat_bootrep




MoosePopR_DomStrat_bootrep <- function(
         stratum.data
        ,selected.unit.data
        ,waypoint.data
        
        ,density=NULL, abundance=NULL, numerator=NULL, denominator=NULL
        
        ,sight.model=NULL, sight.beta=NULL, sight.beta.cov=NULL
        
      ,stratum.var             = "Stratum"   # variable identifying the stratum id in the survey.data and stratum.data
      ,domain.var              = "Domain"    # variable identifying the domain id in the various data frames
      ,stratum.total.blocks.var = "Total.Blocks"  # total number of blocks in the stratum      
      ,stratum.total.area.var  = "Total.Area"    # total area in the stratum      

      ,block.id.var  ="Block.ID" # variable identifying the block   id in the selected.unit.data and waypoint.data frames
      ,block.area.var="Block.Area" # variable identifying the block area measuredin the selected.unit.data frame

      ,conf.level              = 0.90       # confidence interval level
      ,survey.lonely.psu       = "fail",
       check.args=TRUE)
{
#  https://stackoverflow.com/questions/15891459/get-a-list-of-all-function-parameters-from-inside-the-function
   input.data <- as.list(environment(), all=TRUE)

  # check the arguments?
if(check.args){
  # Check that input data are all data frames
  if( !is.data.frame(stratum.data ))     stop("stratum.data is not a data frame")
  if( !is.data.frame(selected.unit.data))stop("selected.unit.data is not a data frame")
  if( !is.data.frame(waypoint.data))     stop("waypoint.data is not a data frame")

# Check the stratum/domain variables  
  if( is.null(stratum.var) || !is.character(stratum.var) || !length(stratum.var)==1)stop("stratum.var is missing or not a character or not length 1")
  if( is.null(domain.var)  || !is.character(domain.var)  || !length(domain.var) ==1)stop("domain.var  is missing or not a character or not length 1")
  
# Make sure that the stratum.var and domain.var is present in all of the data frames
  if( !stratum.var %in% names(stratum.data       ))stop("Stratum variable: ", stratum.var," not in stratum data") 
  if( !stratum.var %in% names(selected.unit.data ))stop("Stratum variable: ", stratum.var," not in selected unit data")
  if( !stratum.var %in% names(waypoint.data      ))stop("Stratum variable: ", stratum.var," not in waypoint data")
  if( !domain.var  %in% names(stratum.data       ))stop("Domain variable:  ", domain.var," not in stratum data") 
  if( !domain.var  %in% names(selected.unit.data ))stop("Domain variable:  ", domain.var," not in selected unit data")
  if( !domain.var  %in% names(waypoint.data      ))stop("Domain variable:  ", domain.var," not in waypoint data")

#-----------------------     
# Check the stratum.data
  if( is.null(stratum.total.blocks.var)|| !is.character(stratum.total.blocks.var)|| !length(stratum.total.blocks.var)==1)
       stop("stratum.blocks.var is missing or not a character or not length 1")
  if( is.null(stratum.total.area.var)|| !is.character(stratum.total.area.var)|| !length(stratum.total.area.var)==1)
       stop("stratum.area.var is missing or not a character or not length 1")
  
# Make sure that stratum number of blocks and stratum area are in the stratum.data df
  if( !stratum.total.blocks.var %in% names(stratum.data))stop("Stratum variable for number of blocks not in stratum df: ", stratum.total.blocks.var)
  if( !stratum.total.area.var   %in% names(stratum.data))stop("Stratum variable for stratum area     not in stratum df: ", stratum.total.area.var)

  #browser()
# Check that total number of blocks and total.area are numbers
  if( !is.numeric(stratum.data[,stratum.total.blocks.var,drop=TRUE]))stop("Stratum number of blocks not numeric")
  if( !is.numeric(stratum.data[,stratum.total.area.var  ,drop=TRUE]))stop("Stratum area not numeric")

# check if missing values>
  if( sum(is.na(stratum.data[,stratum.total.blocks.var]))>0)stop("Stratum humber of blocks cannot be missing ")
  if( sum(is.na(stratum.data[,stratum.total.area.var  ]))>0)stop("Stratum total area cannot be missing ")
  
#------------------  
# check the selected units data frame. This is a list of the sampling units selected in each stratum/domain combination
# along with the area measured.
  if( is.null(block.id.var)  || !is.character(block.id.var)  || !length(block.id.var)  ==1)
       stop("block.id.var is missing or not a character or not length 1")
  if( is.null(block.area.var)|| !is.character(block.area.var)|| !length(block.area.var)==1)
       stop("block.area.var is missing or not a character or not length 1")
  if( !block.id.var   %in% names(selected.unit.data))stop("Block.ID variable:   ", block.id.var,  " not in selected.unit.data")  
  if( !block.area.var %in% names(selected.unit.data))stop("Block.Area variable: ", block.area.var," not in selected.unit.data") 
  
  if( !is.numeric(selected.unit.data[,block.area.var,drop=TRUE]))stop("block area must be numeric")

# any missing values
  if( sum(is.na(selected.unit.data[,block.area.var]))>0)stop("Block area cannot be missing in selected.unit.data ")

#----------------    
# check the waypoint data
  if( !block.id.var  %in% names(waypoint.data))stop("Block.ID variable : ", block.id.var, "not on waypoint.data frame")
  
#----------------------------  
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
  
  if(Type=="D" && !formula.tools::is.one.sided(density  ))stop("Density specification must be a right-sided formula")
  if(Type=="A" && !formula.tools::is.one.sided(abundance))stop("Abundance specification must be a right-sided formulat")
  if(Type=="R" && (!formula.tools::is.one.sided(numerator) | !formula.tools::is.one.sided(denominator)))
     stop("Abundance specification must be a right-sided formulat")
     
  # extract the variables from the formula and check that valid
  if(Type=="D"){
    density.var = formula.tools::rhs.vars(density)
    if(length(density.var)>1)stop("Only one variable for density formula")
    if(!density.var %in% names(waypoint.data))stop("Density variable not in waypoint data for ",density.var)
    if(!is.numeric(waypoint.data[, density.var,drop=TRUE]))stop("Density variable in waypoint.data not numeric for ", density.var)
    if(any(is.na(waypoint.data[,density.var,drop=TRUE])))stop("Missing data not allowed in waypoint.data data values for ", density.var)
  }
  if(Type=="A"){
    abundance.var = formula.tools::rhs.vars(abundance)
    if(length(abundance.var)>1)stop("Only one variable for abundance formula")
    if(!abundance.var %in% names(waypoint.data))stop("abundance variable not in waypoint.data for ",abundance.var)
    if(!is.numeric(waypoint.data[, abundance.var,drop=TRUE]))stop("abundance variable in waypoint.data not numeric for ", abundance.var)
    if(any(is.na(waypoint.data[,abundance.var,drop=TRUE])))stop("Missing data not allowed in waypoint.data data values for ", abundance.var)
  }
  if(Type=="R"){
    numerator.var = formula.tools::rhs.vars(numerator)
    if(length(numerator.var)>1)stop("Only one variable for numerator formula")
    if(!numerator.var %in% names(waypoint.data))stop("numerator variable not in waypoint.data for ",numerator.var)
    if(!is.numeric(waypoint.data[, numerator.var,drop=TRUE]))stop("numerator variable in waypoint.data not numeric for ", numerator.var)
    if(any(is.na(waypoint.data[,numerator.var,drop=TRUE])))stop("Missing data not allowed in waypoint.data data values for ", numerator.var)
    denominator.var = formula.tools::rhs.vars(denominator)
    if(length(denominator.var)>1)stop("Only one variable for denominator formula")
    if(!denominator.var %in% names(waypoint.data))stop("denominator variable not in waypoint.data for ",denominator.var)
    if(!is.numeric(waypoint.data[, denominator.var,drop=TRUE]))stop("denominator variable in waypoint.data not numeric for ", denominator.var)
    if(any(is.na(waypoint.data[,denominator.var,drop=TRUE])))stop("Missing data not allowed in waypoint.data data values for ", denominator.var)
  }
  
# make sure confidence level is sensible.
  if(conf.level < 0.5 | conf.level > .999)stop("Confidence level not sensible: ", conf.level)
  
# check the sightability model if specifiec
  if(!is.null(sight.model)){
    check.sightability.model.args(waypoint.data, sight.model=sight.model, sight.beta=sight.beta, sight.beta.cov=sight.beta.cov)
  }
} # end of argument checking
  
   # Which strata are 100% sampling, typically only 1 unit in stratum total table
   census.strata <- stratum.data[,stratum.var,drop=TRUE][ stratum.data[,stratum.total.blocks.var,drop=TRUE]==1 ]
   
   # We will bootstrap the S1/S2 data separately
   # First determine which sample units have s1 only and which have s1/s2
   ..Domain <- NULL # avoid R CMD check errors
   selected.unit.data$..Domain <- selected.unit.data[,domain.var,drop=TRUE]
   SU.Selected.ndomains <- plyr::ddply(selected.unit.data, c(stratum.var,block.id.var), plyr::summarize,
                                       n.domain= length(unique(..Domain)))
   SU.boot <- plyr::ddply(SU.Selected.ndomains, c(stratum.var,"n.domain"), function(x){
        boot.sample <- sample(1:nrow(x), size=nrow(x), replace=TRUE)  
        x <- x[boot.sample,]
        x$SamplingUnitID.new <- paste0(1:nrow(x),"..",x[,block.id.var,drop=TRUE], "..", x$n.domain[1])
        x$SamplingUnitID.new <- paste0(x$SamplingUnitID.new, 
                                       ifelse((x[,stratum.var,drop=TRUE] %in% census.strata),"..census",""))
        x
   })
   
   # Generate a new set of selected.units and waypoints for each element of the bootstrap sample
   selected.unit.boot <- plyr::adply(SU.boot, 1, function(x){
       SU.select <- selected.unit.data[ selected.unit.data[, block.id.var,drop=TRUE] == x[,block.id.var,drop=TRUE],]
       SU.select$SamplingUnitID.old <- SU.select[, block.id.var,drop=TRUE]
       SU.select$SamplingUnitID     <- x$SamplingUnitID.new
       #browser()
       SU.select
   })

   waypoint.boot <- plyr::adply(SU.boot,1, function(x){
       WP.select <- waypoint.data[ waypoint.data[,block.id.var,drop=TRUE] == x$SamplingUnitID,]
       WP.select$SamplingUnitID.old <- WP.select[, block.id.var,drop=TRUE]
       WP.select$SamplingUnitID <- x$SamplingUnitID.new
       WP.select
   })
   #browser()
   # adjust for sightability if the sighability model is specified
   # assume a multivariate normal around the current beta with the specified covariance matrix
   # then generate the density, abundance, numerator, denominator variables using the SCF
   sight.beta.boot = sight.beta
   if(!is.null(sight.model)){
      sight.beta.boot <- as.vector(mvtnorm::rmvnorm(1, mean=sight.beta,
                                             sigma=sight.beta.cov))
      #browser()
      waypoint.boot$SCF <- compute.SCF(waypoint.boot, 
                                          sight.model    = sight.model,
                                          sight.beta     = sight.beta.boot,
                                          sight.beta.cov = sight.beta.cov)
      if(!is.null(density)){
         var <- formula.tools::rhs.vars(density)
         waypoint.boot[,var] <- waypoint.boot[,var] * waypoint.boot$SCF
      }
      if(!is.null(abundance)){
         var <- formula.tools::rhs.vars(abundance)
         waypoint.boot[,var] <- waypoint.boot[,var] * waypoint.boot$SCF
      }
      if(!is.null(numerator)){
         var <- formula.tools::rhs.vars(numerator)
         waypoint.boot[,var] <- waypoint.boot[,var] * waypoint.boot$SCF
      }
      if(!is.null(denominator)){
         var <- formula.tools::rhs.vars(denominator)
         waypoint.boot[,var] <- waypoint.boot[,var] * waypoint.boot$SCF
      }
   }
   
   boot.data <- list(
      stratum.boot      =stratum.data,
      selected.unit.boot=selected.unit.boot,
      waypoint.boot     =waypoint.boot,
      sight.beta.boot   =sight.beta.boot
   )
   
   res <- MoosePopR_DomStrat(
                    stratum.data         =stratum.data
                   ,selected.unit.data   =selected.unit.boot # units that were selected
                   ,waypoint.data        =waypoint.boot      # waypoint data 
                    
                   ,density              =density, 
                    abundance            =abundance, 
                    numerator            =numerator, 
                    denominator          =denominator,
                    
                    stratum.var          =stratum.var,
                    domain.var           =domain.var
                   ,stratum.total.blocks.var= stratum.total.blocks.var  # total number of blocks in the stratum      
                   ,stratum.total.area.var  =stratum.total.area.var # total area of blocks in the stratum
                   
                   ,block.id.var            ="SamplingUnitID"       # newly created sampling units numbers
                   ,block.area.var          =block.area.var
                   ,conf.level              =conf.level
                   ,survey.lonely.psu       =survey.lonely.psu
 
)

   list(input.data=input.data,
        boot.data =boot.data,
        boot.res  =res
   )
}


