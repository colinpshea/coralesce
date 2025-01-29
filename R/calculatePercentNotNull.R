#' Calculate percent not null
#' @description This function calculated the percentage of columns (loci) that are NOT NA values i.e., that have valid SNP data. 
#' @export
calcPercentNotNull <- function(x) 100 - mean(is.na(x)*100)
