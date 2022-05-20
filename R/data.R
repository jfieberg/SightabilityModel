

#' Experimental (test trials) data set used to estimate detection probabilities
#' for moose in MN
#' 
#' Experimental (test trials) data set used to estimate detection probabilities
#' for moose in MN
#' 
#' 
#' @name exp.m
#' @docType data
#' @format A data frame with 124 observations on the following 4 variables.
#' \describe{ 
#'   \item{year}{year of the experimental survey (test trial)}
#'   \item{observed}{Boolean variable (=1 if moose was observed and 0 otherwise)} 
#'   \item{voc}{measurement of visual obstruction}
#'   \item{grpsize}{group size (number of observed moose in each independently sighted group)}
#'  }
#' @references Giudice, J H. and Fieberg, J. and Lenarz, M. S.  2012. Spending
#' Degrees of Freedom in a Poor Economy: A Case Study of Building a
#' Sightability Model for Moose in Northeastern Minnesota.  Journal of Wildlife
#' Management 76(1):75-87.
#' @keywords datasets
#' @examples
#' 
#' data(exp.m)
#' exp.m[1:5,]
#' 
NULL





#' Mountain Goat Sightability Model Information
#' 
#' Model averaged regression parameters and unconditional variance-covariance
#' matrix for mountain goat sightability model (Rice et al. 2009)
#' 
#' 
#' @name g.fit
#' @docType data
#' @format The format is: beta.g = list of regression parameters (intercept and
#' parameters associated with GroupSize, Terrain, and X.VegCover) varbeta.g =
#' variance-covariance matrix (associated with beta.g)
#' @references Rice C.G., Jenkins K.J., Chang W.Y. (2009). A Sightability Model
#' for Mountain Goats. The Journal of Wildlife Management, 73(3), 468-478.
#' @keywords datasets
#' @examples
#' 
#' data(g.fit)
#' 
NULL





#' Mountain Goat Survey Data from Olympic National park
#' 
#' Mountain Goat Survey Data from Olympic National park collected in 2004
#' 
#' 
#' @name gdat
#' @docType data
#' @format A data frame with 113 observations on the following 9 variables.
#' \describe{ 
#'   \item{GroupSize}{number of animals observed in each
#'         independently sighted group [cluster size]} 
#'   \item{Terrain}{measure
#'         of terrain obstruction} \item{pct.VegCover}{measure of vegetative
#'         obstruction} \item{stratum}{stratum identifier}
#'   \item{total}{number of animals observed in each independently
#'         sighted group [same as GroupSize]} 
#'   \item{subunit}{a numeric vector, Plot ID} 
#' }
#' @references Jenkins, K. J., Happe, P.J., Beirne, K.F, Hoffman, R.A.,
#' Griffin, P.C., Baccus, W. T., and J. Fieberg.  In press.  Recent population
#' trends in mountain goats in the Olympic mountains.  Northwest Science.
#' @source Patti Happe (Patti_Happe@nps.gov)
#' @keywords datasets
#' @examples
#' 
#' data(gdat)
#' 
NULL





#' MN moose survey data
#' 
#' Operational survey data for moose in MN (during years 2004-2007).  Each
#' record corresponds to an independently sighted group of moose, with variables
#' that capture individual covariates (used in the detection model) as well as
#' plot-level information (stratum identifier, sampling probability, etc).
#' 
#' 
#' @name obs.m
#' @docType data
#' @format A data frame with 805 observations on the following 11 variables.
#' \describe{ 
#'   \item{year}{year of survey}
#'   \item{stratum}{stratum identifier} \item{subunit}{sample
#'        plot ID} 
#'   \item{total}{number of moose observed}
#'   \item{cows}{number of cows observed} \item{calves}{number of
#'        calves observed} 
#'   \item{bulls}{number of bulls observed}
#'   \item{unclass}{number of unclassified animals observed (could not
#'         identify sex/age class)} 
#'   \item{voc}{ measurement of visual
#'        obstruction} 
#'   \item{grpsize}{group size (cluter size)} 
#' }
#' @references Giudice, J H. and Fieberg, J. and Lenarz, M. S.  2012. Spending
#' Degrees of Freedom in a Poor Economy: A Case Study of Building a
#' Sightability Model for Moose in Northeastern Minnesota.  Journal of Wildlife
#' Management 76(1):75-87.
#' @keywords datasets
#' @examples
#' 
#' data(obs.m) 
#' obs.m[1:5, ]
#' 
NULL





#' Data set containing sampling information for observation survey of moose in
#' MN
#' 
#' Data set containing sampling information from a survey of moose in MN
#' (during years 2004-2007)
#' 
#' 
#' @name sampinfo.m
#' @docType data
#' @format A data frame with 12 observations on the following 5 variables.
#' \describe{ 
#'   \item{year}{year of survey}
#'   \item{stratum}{stratum identifier} 
#'   \item{Nh}{number of
#'        population units in stratum h} 
#'   \item{nh}{number of sample units in
#'        stratum h} 
#' }
#' @references Giudice, J H. and Fieberg, J. and Lenarz, M. S.  2012. Spending
#' Degrees of Freedom in a Poor Economy: A Case Study of Building a
#' Sightability Model for Moose in Northeastern Minnesota.  Journal of Wildlife
#' Management 76(1):75-87.
#' @keywords datasets
#' @examples
#' 
#' data(sampinfo.m)
#' sampinfo.m
#' 
NULL





#' Wildlife Sightability Modeling
#' 
#' Uses logistic regression to model the probability of detection as a function
#' of covariates.  This model is then used with observational survey data to
#' estimate population size, while accounting for uncertain detection.  See
#' Steinhorst and Samuel (1989).
#' 
#' \tabular{ll}{ Package: \tab SightabilityModel\cr Type: \tab Package\cr
#' Version: \tab 1.3\cr Date: \tab 2014-10-03\cr License: \tab GPL-2\cr
#' LazyLoad: \tab yes\cr }
#' 
#' @name SightabilityModel-package
#' @aliases SightabilityModel-package SightabilityModel
#' @docType package
#' @author John Fieberg
#' 
#' Maintainer: John Fieberg <jfieberg@@umn.edu>
#' @references Fieberg, J.  2012.  Estimating Population Abundance Using
#' Sightability Models: R SightabilityModel Package. Journal of Statistical
#' Software, 51(9), 1-20.  URL https://doi.org/10.18637/jss.v051.i09 
#' 
#' Steinhorst, Kirk R. and Samuel, Michael D. 1989. Sightability Adjustment
#' Methods for Aerial Surveys of Wildlife Populations. Biometrics 45:415--425.
#' @keywords package
NULL



