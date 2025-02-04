#' Calculate percent not null
#' @description This function calculates the percentage of columns (loci) that are NOT NA values i.e., that have valid SNP data.
#' @param x Takes data (usually a vector) and calculated the percentage of observations that are not NA values. This is used to calculate the percentage of loci with valid SNP data.   
#' @export
calcPercentNotNull <- function(x) 100 - mean(is.na(x)*100)
