#' Remove invariant loci
#'
#' @description Drops loci (columns) that carry only a single distinct allele
#'   across colonies, as required before kinship calculations. All-`NA` loci are
#'   retained (they contribute nothing to kinship because such comparisons are
#'   never fully scorable) to preserve the original behaviour.
#' @param dataset A data frame with `Coral_ID` in column 1 and single-letter
#'   allele codes in the remaining columns.
#' @returns The input data frame with invariant loci removed.
#' @export
omitInvariantLoci <- function(dataset) {
  lociCols <- names(dataset)[-1]
  nDistinct <- vapply(
    dataset[lociCols],
    function(col) length(unique(col[!is.na(col)])),
    integer(1)
  )
  keep <- lociCols[nDistinct != 1]                     # drop exactly-invariant loci
  dataset[, c(names(dataset)[1], keep), drop = FALSE]
}
