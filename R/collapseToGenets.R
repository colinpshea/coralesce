#' Collapse colonies to one representative per genet
#'
#' @description Reduces a genotype data frame to a single representative colony
#'   per genet, using a set of genet assignments from [runGenets()]. This is a
#'   pre-filter for clone-corrected analyses: feeding the result to the kinship
#'   pipeline yields genet-level (clone-corrected) gene diversity rather than the
#'   eligible-pool diversity computed over all colonies (see [runKinship()]).
#'
#'   The function only filters rows; it does not alter genotype values, so it can
#'   be applied to either raw paired-allele data (e.g., from [readGeneticData()])
#'   or single-letter coded data (from [convertBasePairstoCodes()]). It matches
#'   rows on `Coral_ID`, so identifiers may contain any characters.
#' @param dataset A data frame with a `Coral_ID` column (plus any genotype
#'   columns). Rows are selected from this frame.
#' @param genetAssignment A genet-assignment data frame with at least `Coral_ID`
#'   and `genet` columns (and, ideally, `pctNull`), as returned per file by
#'   [runGenets()].
#' @param representative How to choose the representative colony within a genet:
#'   `"most_data"` (default) keeps the colony with the most scored loci (lowest
#'   `pctNull`); `"first"` keeps the first colony listed. `"most_data"` falls
#'   back to `"first"` if no `pctNull` column is present.
#' @param drop_unassigned Logical. Colonies without a real genet (missing genet,
#'   or a placeholder genet label ending in `NA` — i.e., inadequate-data or
#'   all-`NA` colonies) are not clonal duplicates of one another. If `FALSE`
#'   (default) each is retained as its own unit; if `TRUE` they are dropped.
#' @returns A subset of `dataset` containing one row per genet plus (unless
#'   dropped) each unassigned colony. Column structure and order are unchanged.
#' @seealso [runGenets()] to produce `genetAssignment`; [runKinship()] and
#'   [kinshipCalcsNoInvar()] for the kinship pipeline this feeds.
#' @examples
#' \dontrun{
#' # Clone-corrected (genet-level) kinship and gene diversity:
#' g   <- runGenets(PctMatchThreshold = 90, PctNotNullThreshold = 50)
#' ga  <- g$genetAssignments[[1]]
#' raw <- readGeneticData(file.path(getwd(), "Data", "myfile.csv"))[[1]]
#'
#' reduced <- collapseToGenets(raw, ga)                 # one colony per genet
#'
#' coded   <- convertBasePairstoCodes(reduced)
#' withDat <- isolateAllNAColonies(coded)[[1]]
#' noInv   <- omitInvariantLoci(withDat)
#' k       <- kinshipCalcsNoInvar(noInv, subset = FALSE)
#' k$PopAvgMKGD                                          # clone-corrected diversity
#' }
#' @export
collapseToGenets <- function(dataset, genetAssignment,
                             representative = c("most_data", "first"),
                             drop_unassigned = FALSE) {
  representative <- match.arg(representative)

  if (!("Coral_ID" %in% names(dataset))) {
    stop("`dataset` must contain a `Coral_ID` column.", call. = FALSE)
  }
  if (!all(c("Coral_ID", "genet") %in% names(genetAssignment))) {
    stop("`genetAssignment` must contain `Coral_ID` and `genet` columns.",
         call. = FALSE)
  }

  ga <- genetAssignment
  ga$genet <- as.character(ga$genet)

  # A colony is "unassigned" if it has no real genet: NA, or a placeholder
  # label ending in NA (inadequate-data and all-NA colonies).
  unassigned <- is.na(ga$genet) | grepl("NA$", ga$genet)
  assigned   <- ga[!unassigned, , drop = FALSE]

  # One representative Coral_ID per genet.
  pick_rep <- function(df) {
    if (representative == "most_data" && "pctNull" %in% names(df)) {
      df <- df[order(df$pctNull), , drop = FALSE]   # lowest pctNull = most data
    }
    df$Coral_ID[1]
  }
  reps <- unlist(lapply(split(assigned, assigned$genet), pick_rep),
                 use.names = FALSE)

  keep_ids <- reps
  if (!drop_unassigned) {
    keep_ids <- c(keep_ids, ga$Coral_ID[unassigned])
  }

  dataset[dataset$Coral_ID %in% keep_ids, , drop = FALSE]
}
