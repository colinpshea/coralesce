#' Assign colonies to genets from a genotype data frame
#'
#' @description Data-frame entry point to the genet-assignment pipeline: takes a
#'   genotype data frame directly (no `Data`/`Results` folders) and returns genet
#'   assignments. This is the in-memory counterpart to [runGenets()], useful when
#'   the genotypes are already in R -- in particular, its output is the
#'   `genetAssignment` that [collapseToGenets()] needs, so the whole
#'   clone-corrected workflow can run without files.
#'
#'   Internally it runs the same steps as [runGenets()]: translate paired alleles
#'   to IUPAC codes, set aside all-`NA` colonies, compare all colony pairs, and
#'   assign genets ([groupByGenets()]). Results are keyed by `Coral_ID`.
#' @param data A genotype data frame with a `Coral_ID` column and paired-allele
#'   SNP columns (e.g., `"A:G"`), as produced by [readGeneticData()]. A
#'   `MatchMaker_Index` column, if present, is carried through and used to order
#'   the output. Missing values may be `NA` or `"?"`.
#' @param PctMatchThreshold Minimum percent allele match across loci for two
#'   colonies to be called a clone. Required.
#' @param PctNotNullThreshold Minimum percent of loci with valid data required to
#'   trust a comparison. Required.
#' @param speciesCode Optional short code prefixed to each genet label (e.g.,
#'   `"ACER"` gives `ACER_00001`). If `NULL` (default) genet labels are the
#'   zero-padded number alone (`00001`). Colonies without a genet (inadequate or
#'   all-`NA` data) carry a label ending in `NA`.
#' @param getPairwiseAlleleMatches Logical; if `TRUE`, also return the full
#'   pairwise comparison table. Default `FALSE`.
#' @returns If `getPairwiseAlleleMatches = FALSE`, a genet-assignment data frame
#'   with columns `Coral_ID`, (`MatchMaker_Index` if supplied,) `genet`,
#'   `pctNull`, `AdequateData`. If `TRUE`, a list of that data frame plus the
#'   pairwise table.
#' @seealso [runGenets()] for the folder-based workflow; [collapseToGenets()] to
#'   reduce to one colony per genet; [computeKinship()] for the kinship pipeline.
#' @examples
#' \dontrun{
#' raw <- readGeneticData(file.path(getwd(), "Data", "myfile.csv"))[[1]]
#' ga  <- computeGenets(raw, PctMatchThreshold = 90, PctNotNullThreshold = 50)
#'
#' # fully folder-free clone-corrected diversity:
#' computeKinship(collapseToGenets(raw, ga))
#' }
#' @importFrom dplyr bind_rows left_join arrange mutate select
#' @importFrom stringr str_pad
#' @importFrom magrittr %>%
#' @export
computeGenets <- function(data, PctMatchThreshold = NULL, PctNotNullThreshold = NULL,
                          speciesCode = NULL, getPairwiseAlleleMatches = FALSE) {

  if (!("Coral_ID" %in% names(data))) {
    stop("`data` must contain a `Coral_ID` column.", call. = FALSE)
  }
  find_dups(data)
  if (is.null(PctMatchThreshold) || is.null(PctNotNullThreshold)) {
    stop("Both PctMatchThreshold and PctNotNullThreshold must be supplied.",
         call. = FALSE)
  }

  index <- if ("MatchMaker_Index" %in% names(data)) {
    data[, c("Coral_ID", "MatchMaker_Index")]
  } else {
    NULL
  }

  # Record missing values as "?" so partially-scored loci survive translation.
  alleleCols <- setdiff(names(data), c("Coral_ID", "MatchMaker_Index"))
  if (length(alleleCols) > 0) {
    sub <- data[alleleCols]
    sub[is.na(sub)] <- "?"
    data[alleleCols] <- sub
  }

  hasLoci <- length(alleleCols) > 0 &&
    any(vapply(data[alleleCols],
               function(col) all(col %in% IUPAC$Allelepairs), logical(1)))
  if (!hasLoci) {
    stop("No valid allele-pair loci found in `data`. Genotypes must be paired ",
         "alleles such as \"A:G\", not single-letter codes.", call. = FALSE)
  }

  split    <- isolateAllNAColonies(convertBasePairstoCodes(initdata = data))
  withData <- split[[1]]
  noData   <- split[[2]]

  matches <- determineAllAlleleMatches(dataset = withData)
  gg <- groupByGenets(
    CoralAlleleData          = withData,
    AlleleMatchResults       = matches,
    PctMatchThreshold        = PctMatchThreshold,
    PctNotNullThreshold      = PctNotNullThreshold,
    getPairwiseAlleleMatches = getPairwiseAlleleMatches
  )

  label <- function(g) {
    padded <- str_pad(g, 5, "left", "0")
    if (!is.null(speciesCode)) paste0(speciesCode, "_", padded) else padded
  }

  assignment <- gg$genetAssignment %>%
    bind_rows(noData) %>%
    mutate(genet = label(genet))

  if (!is.null(index)) {
    assignment <- assignment %>%
      left_join(index, by = "Coral_ID") %>%
      arrange(MatchMaker_Index) %>%
      select(Coral_ID, MatchMaker_Index, genet, pctNull, AdequateData)
  } else {
    assignment <- assignment %>%
      select(Coral_ID, genet, pctNull, AdequateData)
  }

  if (isTRUE(getPairwiseAlleleMatches)) {
    list(genetAssignment = assignment, pairwiseAlleleMatches = gg$pairwiseAlleleMatches)
  } else {
    assignment
  }
}
