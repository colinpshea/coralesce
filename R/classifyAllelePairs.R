#' Classify allele pairs at a single locus
#'
#' @description Builds all pairwise combinations of colonies at one locus and
#'   flags, for each pair, whether the two colonies share the same allele. By
#'   default comparisons of each colony with itself are included (needed for
#'   genet assignment); set `include_self = FALSE` to return only comparisons
#'   with other colonies (used for kinship). Operates on one locus (column) at
#'   a time.
#' @param dataset A data frame with `Coral_ID` in column 1 and single-letter
#'   `IUPAC` allele data for one locus per remaining column.
#' @param locus The column index of the locus to classify.
#' @param include_self Logical; include self-comparisons (default `TRUE`).
#' @returns A data frame with columns `coral1`, `coral2`, `allele1`, `allele2`,
#'   `locus`, and `match` (`TRUE`/`FALSE`, or `NA` when either allele is
#'   missing).
#' @export
classifyAllelePairs <- function(dataset, locus, include_self = TRUE) {
  locus_name <- names(dataset)[locus]
  ids        <- as.character(dataset[[1]])
  alleles    <- dataset[[locus]]

  # Cross comparisons (unordered pairs of distinct colonies).
  if (length(ids) >= 2) {
    cp    <- t(utils::combn(ids, 2))
    pairs <- data.frame(coral1 = cp[, 1], coral2 = cp[, 2],
                        stringsAsFactors = FALSE)
  } else {
    pairs <- data.frame(coral1 = character(0), coral2 = character(0),
                        stringsAsFactors = FALSE)
  }

  # Self comparisons.
  if (isTRUE(include_self)) {
    pairs <- rbind(pairs,
                   data.frame(coral1 = ids, coral2 = ids,
                              stringsAsFactors = FALSE))
  }

  a <- alleles[match(pairs$coral1, ids)]
  b <- alleles[match(pairs$coral2, ids)]

  data.frame(
    coral1  = pairs$coral1,
    coral2  = pairs$coral2,
    allele1 = a,
    allele2 = b,
    locus   = locus_name,
    match   = a == b,          # NA when either allele is missing
    stringsAsFactors = FALSE
  )
}
