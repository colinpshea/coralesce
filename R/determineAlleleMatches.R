#' Assess allele matches across all loci (self and other comparisons)
#'
#' @description Applies [classifyAllelePairs()] to every locus
#'   (`2:ncol(dataset)`) and stacks the results. Includes both self- and
#'   other-comparisons, as required for genet assignment via [groupByGenets()].
#' @param dataset A data frame with `Coral_ID` in column 1 and single-letter
#'   `IUPAC` allele data in the remaining columns.
#' @returns A data frame with columns `coral1`, `coral2`, `allele1`, `allele2`,
#'   `locus`, and `match`, covering all pairwise comparisons (self and other) at
#'   every locus.
#' @export
determineAllAlleleMatches <- function(dataset) {
  if (ncol(dataset) < 2) {
    stop("`dataset` must have Coral_ID plus at least one locus column.")
  }
  do.call(rbind, lapply(2:ncol(dataset), function(j) {
    classifyAllelePairs(dataset, locus = j, include_self = TRUE)
  }))
}
