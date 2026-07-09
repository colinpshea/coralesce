#' Assign colonies to genets from pairwise allele matches
#'
#' @description For every pairwise comparison of colonies, computes the percent
#'   of loci that match and the percent of loci with scorable data, applies the
#'   user thresholds to decide which pairs are clones with adequate data, and
#'   assigns colonies to genets via [returnGenetIdentity()]. Colonies whose
#'   self-comparison falls below the data threshold are flagged
#'   `AdequateData = No` and carried through with `genet = NA`.
#' @param CoralAlleleData Wide single-letter allele data (`Coral_ID` plus one
#'   column per locus); used to attach each colony's `pctNull`.
#' @param AlleleMatchResults Pairwise allele-match results from
#'   [determineAllAlleleMatches()].
#' @param PctMatchThreshold Minimum percent allele match for a pair to be called
#'   a clone. Required.
#' @param PctNotNullThreshold Minimum percent of scorable loci required to trust
#'   a comparison. Required.
#' @param getPairwiseAlleleMatches Logical; if `TRUE`, also return the full
#'   pairwise comparison table. Default `FALSE`.
#' @returns A list with `genetAssignment` (columns `Coral_ID`, `genet`,
#'   `pctNull`, `AdequateData`) and `pairwiseAlleleMatches` (the pairwise table,
#'   or `NULL`).
#' @importFrom dplyr select filter mutate if_else left_join arrange rename bind_rows
#' @importFrom tidyr pivot_wider
#' @importFrom magrittr %>%
#' @export
groupByGenets <- function(CoralAlleleData, AlleleMatchResults,
                          PctMatchThreshold = NULL, PctNotNullThreshold = NULL,
                          getPairwiseAlleleMatches = FALSE) {

  if (is.null(PctMatchThreshold) || is.null(PctNotNullThreshold)) {
    stop("Both PctMatchThreshold and PctNotNullThreshold must be supplied.",
         call. = FALSE)
  }

  # Percent of loci with no data, per colony.
  CoralAlleleData$pctNull <-
    100 - apply(subset(CoralAlleleData, select = -Coral_ID), 1, calcPercentNotNull)
  CoralAlleleData <- CoralAlleleData %>% select(Coral_ID, pctNull)

  # One row per colony pair; locus matches spread across columns.
  temp <- AlleleMatchResults %>%
    select(coral1, coral2, locus, match) %>%
    arrange(coral1, coral2, locus) %>%
    pivot_wider(names_from = locus, values_from = match)

  locusCols       <- setdiff(names(temp), c("coral1", "coral2"))
  temp$pctMatch   <- rowMeans(as.matrix(temp[, locusCols]), na.rm = TRUE) * 100
  temp$pctNotNull <- apply(temp[, locusCols], 1, calcPercentNotNull)

  temp <- temp %>%
    mutate(PartOfGenet = if_else(pctMatch >= PctMatchThreshold, "Yes", "No")) %>%
    select(coral1, coral2, pctMatch, pctNotNull, PartOfGenet)

  # Keep clone pairs; drop cross-pairs with too little data; flag self-pairs
  # that lack adequate data.
  clones <- temp %>%
    filter(PartOfGenet == "Yes") %>%
    mutate(flag = if_else(coral1 != coral2 & pctNotNull < PctNotNullThreshold,
                          "drop", "keep")) %>%
    filter(flag == "keep") %>%
    mutate(AdequateData = if_else(coral1 == coral2 & pctNotNull < PctNotNullThreshold,
                                  "No", "Yes")) %>%
    select(coral1, coral2, pctMatch, pctNotNull, PartOfGenet, AdequateData)

  adequateYes <- clones %>% filter(AdequateData == "Yes")

  adequateNo <- clones %>%
    filter(AdequateData == "No") %>%
    mutate(genet = NA_integer_, pctNull = 100 - pctNotNull) %>%
    select(coral1, genet, pctNull, AdequateData) %>%
    rename(Coral_ID = coral1)

  genetAssignment <- returnGenetIdentity(adequateYes) %>%
    mutate(AdequateData = "Yes") %>%
    arrange(genet) %>%
    left_join(CoralAlleleData, by = "Coral_ID") %>%
    select(Coral_ID, genet, pctNull, AdequateData) %>%
    bind_rows(adequateNo)

  pairwise <- if (isTRUE(getPairwiseAlleleMatches)) temp else NULL
  list(genetAssignment = genetAssignment, pairwiseAlleleMatches = pairwise)
}
