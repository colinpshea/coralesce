#' Classify allele pairs at a single locus (excluding self-comparisons)
#'
#' @description Convenience wrapper around [classifyAllelePairs()] that returns
#'   only comparisons of each colony with other colonies (no self-comparisons).
#'   Used by the kinship pipeline, which must exclude within-colony comparisons.
#' @inheritParams classifyAllelePairs
#' @returns A data frame with columns `coral1`, `coral2`, `allele1`, `allele2`,
#'   `locus`, and `match`; self-comparisons are omitted.
#' @export
classifyAllelePairsOthers <- function(dataset, locus) {
  classifyAllelePairs(dataset, locus, include_self = FALSE)
}
