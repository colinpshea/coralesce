#' Calculate percent not null
#' @description This function calculates the percentage of columns (loci) that are not `NA` or `NULL` values i.e., that have valid SNP data.
#' @param x Any data frame, matrix, vector.
#' @returns This function takes data (a vector, data frame, or matrix) and calculates the percentage of observations that are not `NA` or `NULL` values. This is used to calculate the percentage of loci with valid SNP data.  
#' @export
calcPercentNotNull <- function(x) 100 - mean(is.na(x)*100)
