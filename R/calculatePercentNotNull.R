#' Calculate percent not null
#'
#' @description Calculates the percentage of values in `x` that are not `NA`
#'   (i.e., that carry valid SNP data). Used to summarise the proportion of loci
#'   with usable data for a colony or a pairwise comparison.
#' @param x A vector, matrix, or data frame.
#' @returns A single numeric value: the percentage of non-`NA` values in `x`.
#' @export
calcPercentNotNull <- function(x) 100 * mean(!is.na(x))
