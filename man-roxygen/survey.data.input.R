#' @param survey.data A data frame containing counts of moose in each group along with a variable identifying 
#'    the stratum (see stratum.var) and block (see block.id.var)
#' @param survey.block.area A data frame containing for each block, the block id (see block.id.var), the
#'    area of the block (see block.area.var). The data frame can contain information for other blocks that
#'    were not surveyed (e.g. for the entire population of blocks) and additional block information will be ignored.
#' @param stratum.data A data frame containing for each stratum, the stratum id (see stratum.var), the total
#'    number of blocks in the stratum (see stratum.blocks.var) and the total area of the stratum (see stratum.area.var)
#' @param density,abundance,numerator,denominator Right-handed formula identifying the variable(s) in the
#'    survey.data data frame for which the density, abundance, or ratio (numerator/denominator) are to be estimated.
#' @param block.id.var Name of the variable in the survey.data data frame and survey.block.area data frame 
#'    that identifies the block.id 
#'    that links the block between the survey data and the block information.
#' @param block.area.var Name of the variable in the survey.block.area data frame that contains the area of the blocks.
#' @param stratum.var Name of the variable in the survey.data data frame and thee stratum.data data frame that identifies
#'    the stratum.
#' @param stratum.blocks.var Name of the variable in the stratum.data data frame that contains the total number of
#'    blocks in the stratum.
#' @param stratum.area.var  Name of the variable in the stratum.data data.frame that contains the total stratum area.
#'
