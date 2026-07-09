#' Assess allele matches across all loci (other comparisons only)
#'
#' @description Applies [classifyAllelePairsOthers()] to every locus
#'   (`2:ncol(dataset)`) and stacks the results. Excludes self-comparisons, as
#'   required for kinship via [kinshipCalcsNoInvar()].
#' @param dataset A data frame with `Coral_ID` in column 1 and single-letter
#'   `IUPAC` allele data in the remaining columns.
#' @returns A data frame with columns `coral1`, `coral2`, `allele1`, `allele2`,
#'   `locus`, and `match` for all comparisons of colonies with other colonies.
#' @export
determineAllAlleleMatchesOthers <- function(dataset) {
  if (ncol(dataset) < 2) {
    stop("`dataset` must have Coral_ID plus at least one locus column.")
  }
  do.call(rbind, lapply(2:ncol(dataset), function(j) {
    classifyAllelePairsOthers(dataset, locus = j)
  }))
}
