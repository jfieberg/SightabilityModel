#' Classical and Domain Stratification using MoosePopR()
#' 
#' This function allows for classical or domain stratification when using MoosePopR().
#' Caution **SE are NOT adjusted for measurements on multiple domains on the same
#' sampling unit. Bootstrapping may be required**. Consult the vignette for more details.
#'  
#'  MoosePopR_DomStrat() assumes that sightability is 100\%.
#'  Use the SightabilityPopR_DomStrat() function to adjust for sightability < 100\%.
#'

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
#' @template conf.level
#' @param survey.lonely.psu How to deal with lonely PSU within strata. See \code{surveyoptions} in the \code{survey} package.
#'  
#' @return  A data frame containing for each stratum and for all combinations of strata and domains
#'  (identified as stratum id \code{.OVERALL}), the density,
#'    or abundance or ratio estimate along with its estimated standard error and large-sample normal-based confidence interval.
#' @template author 
#' @template references
#' @keywords ~MOOSEPOP ~moose surveys
#' @import formula.tools
#' @examples
#'  
#' ##---- See the vignettes for examples on how to run this analysis.
#' 
#' @export MoosePopR_DomStrat

MoosePopR_DomStrat <- function(
       stratum.data       # information on strata totals such as number of blocks and total area 
      ,selected.unit.data # units that were selected
      ,waypoint.data      # waypoint data 
      
      ,density=NULL
      ,abundance=NULL
      ,numerator=NULL
      ,denominator=NULL
                           
      ,stratum.var             = "Stratum"   # variable identifying the stratum id in the survey.data and stratum.data
      ,domain.var              = "Domain"    # variable identifying the domain id in the various data frames
      ,stratum.total.blocks.var= "Total.Blocks"  # total number of blocks in the stratum      
      ,stratum.total.area.var  = "Total.Area"    # total area in the stratum      

      ,block.id.var  ="Block.ID" # variable identifying the block   id in the selected.unit.data and waypoint.data frames
      ,block.area.var="Block.Area" # variable identifying the block area measured in the selected.unit.data frame

      ,conf.level              = 0.90       # confidence interval level
      ,survey.lonely.psu       = "fail")
  {

  version <- "2022-06-01"

# Error checking
 
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

  if(any(duplicated(selected.unit.data[,block.id.var,drop=TRUE])))warning("SE not adjusted for measurements on multiple domains on same sampling unit")

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

   #browser()
   # we create a separate stratum for each combination of stratum and domain
   # get the data ready for the call to MoosePopR
   stratum.data$Stratum.Domain <- paste0(stratum.data[,stratum.var,drop=TRUE], "..", stratum.data[,domain.var,drop=TRUE])
   stratum.data$Stratum.Blocks <- stratum.data[,stratum.total.blocks.var, drop=TRUE]
   stratum.data$Stratum.Area   <- stratum.data[,stratum.total.area.var,   drop=TRUE  ]
   
   # get the selected.units.data.frame ready for call to MoosePop
   survey.block.area <- selected.unit.data
   survey.block.area$Stratum.Domain <- paste0(survey.block.area[,stratum.var,drop=TRUE], "..", survey.block.area [,domain.var,drop=TRUE])
   survey.block.area$Block.ID       <- paste0(survey.block.area$Stratum.Domain, "..", survey.block.area[, block.id.var,drop=TRUE])
   survey.block.area$Block.Area     <- survey.block.area[,block.area.var,drop=TRUE]
   

   survey.data <- waypoint.data
   survey.data$Stratum.Domain <- paste0(survey.data[,stratum.var,drop=TRUE], "..", survey.data[, domain.var,drop=TRUE])
   survey.data$Block.ID <-       paste0(survey.data$Stratum.Domain, "..", survey.data[,block.id.var,drop=TRUE])
    
   #browser()
   # check that stratum.domain matches across the stratum.data and selected.unit.data
   t1 <- setdiff(stratum.data$Stratum.Domain, survey.block.area$Stratum.Domain)
   if(length(t1)>0)stop("Stratum.Domain in stratum data, but not in selected.unit.data ", t1)
   t1 <- setdiff(stratum.data$Stratum.Domain, survey.data$Stratum.Domain)
   if(length(t1)>0)stop("Stratum.Domain in stratum data, but in waypoint.data ", t1)
  
   t1 <- setdiff(survey.block.area$Stratum.Domain, stratum.data$Stratum.Domain)
   if(length(t1)>0)stop("Stratum.Domain in selected.unit.data not in stratum.data ", t1)
   t1 <- setdiff(survey.block.area$Stratum.Domain, survey.data$Stratum.Domain)
   if(length(t1)>0)stop("Stratum.Domain in selected.unit.data not in waypoint.data ", t1)
   
   t1 <- setdiff(survey.data$Stratum.Domain, stratum.data$Stratum.Domain)
   if(length(t1)>0)stop("Stratum.Domain in waypoint data not in stratum.data ", t1)
   t1 <- setdiff(survey.data$Stratum.Domain, survey.block.area$Stratum.Domain)
   if(length(t1)>0)stop("Stratum.Domain in waypoint data not in selected.unit.data ", t1)
   
   # check that block id match back and forth
   t1 <- setdiff(survey.data$Block.ID, survey.block.area$Block.ID)
   if(length(t1)>0)stop("Some units in waypoint.data but not in selected.unit.data", t1)
   t1 <- setdiff(survey.block.area$Block.ID, survey.data$Block.ID)
   if(length(t1)>0)stop("Some units in selected.unit.data but not in waypoint.data", t1, 
                        " Did you forget to add dummy waypoints for blocks with 0 observed animals?")
 
   # pass the minimum variables needed to MoosePopR  
   stratum.data     <- stratum.data[,c("Stratum.Domain","Stratum.Blocks","Stratum.Area")]
   survey.block.area<- survey.block.area[, c("Block.ID","Block.Area")]#, "Stratum.Domain",
   survey.data      <- survey.data[, !names(survey.data) %in% c(stratum.var, domain.var, block.id.var)]

   #browser()
   res <- MoosePopR(survey.data, 
                    survey.block.area    =survey.block.area, 
                    stratum.data         =stratum.data,
                    density              =density, 
                    abundance            =abundance, 
                    numerator            =numerator, 
                    denominator          =denominator,
                    
                    stratum.var          ="Stratum.Domain",
                    
                    survey.lonely.psu=survey.lonely.psu,
                    conf.level=conf.level
          )
   res
}





 
