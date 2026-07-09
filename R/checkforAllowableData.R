#' Flag values that are valid IUPAC allele pairs
#'
#' @description Tests, element by element, whether the supplied values are
#'   allowable allele pairs as listed in the `IUPAC` reference table. Used to
#'   identify which columns of an input data set contain valid SNP data.
#' @param x A vector of values to test (typically one column of the input data).
#' @returns A logical vector the same length as `x`, `TRUE` where the value is a
#'   valid allele pair in `IUPAC$Allelepairs`.
#' @export
checkforAllowableData <- function(x) x %in% IUPAC$Allelepairs
