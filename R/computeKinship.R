#' Compute kinship and gene diversity from a genotype data frame
#'
#' @description Data-frame entry point to the kinship pipeline: takes a genotype
#'   data frame directly (no `Data`/`Results` folders) and returns
#'   individual- and population-level mean kinship and gene diversity. This is
#'   the in-memory counterpart to [runKinship()], useful when the genotypes are
#'   already in R -- for example after collapsing to one colony per genet with
#'   [collapseToGenets()] for a clone-corrected estimate.
#'
#'   Internally it runs the same steps as [runKinship()]: translate paired
#'   alleles to IUPAC codes ([convertBasePairstoCodes()]), set aside all-`NA`
#'   colonies ([isolateAllNAColonies()]), optionally drop invariant loci, and
#'   compute kinship ([kinshipCalcsNoInvar()]). Results are keyed by `Coral_ID`.
#' @param data A genotype data frame with a `Coral_ID` column and paired-allele
#'   SNP columns (e.g., `"A:G"`), as produced by [readGeneticData()] or
#'   [collapseToGenets()]. A `MatchMaker_Index` column, if present, is ignored.
#'   Missing values may be `NA` or `"?"`.
#' @param subset Logical; if `TRUE`, also compute the `targetN`-colony summary.
#'   Requires `targetN`. Default `FALSE`.
#' @param targetN Number of least-related colonies to retain when
#'   `subset = TRUE`. Must be >= 2.
#' @returns A list with `PopAvgMKGD`, `MK_init`, and `MK_final` (see
#'   [kinshipCalcsNoInvar()]).
#' @seealso [runKinship()] for the folder-based workflow; [collapseToGenets()]
#'   to reduce to one colony per genet before computing a clone-corrected result.
#' @examples
#' \dontrun{
#' raw <- readGeneticData(file.path(getwd(), "Data", "myfile.csv"))[[1]]
#'
#' # eligible-pool kinship (all colonies, ramets included):
#' computeKinship(raw)
#'
#' # clone-corrected kinship (one colony per genet):
#' ga <- runGenets(PctMatchThreshold = 90, PctNotNullThreshold = 50)$genetAssignments[[1]]
#' computeKinship(collapseToGenets(raw, ga))
#' }
#' @export
computeKinship <- function(data, subset = FALSE, targetN = NULL) {
  if (!("Coral_ID" %in% names(data))) {
    stop("`data` must contain a `Coral_ID` column.", call. = FALSE)
  }

  # Record missing values as "?" so partially-scored loci survive the IUPAC
  # translation (matches readGeneticData()'s handling). Keys are left alone.
  alleleCols <- setdiff(names(data), c("Coral_ID", "MatchMaker_Index"))
  if (length(alleleCols) > 0) {
    sub <- data[alleleCols]
    sub[is.na(sub)] <- "?"
    data[alleleCols] <- sub
  }

  # Fail helpfully up front if nothing looks like paired-allele data, rather
  # than letting convertBasePairstoCodes() error cryptically on an empty pivot.
  hasLoci <- length(alleleCols) > 0 &&
    any(vapply(data[alleleCols],
               function(col) all(col %in% IUPAC$Allelepairs), logical(1)))
  if (!hasLoci) {
    stop("No valid allele-pair loci found in `data`. Genotypes must be paired ",
         "alleles such as \"A:G\", not single-letter codes.", call. = FALSE)
  }

  coded    <- convertBasePairstoCodes(data)
  withData <- isolateAllNAColonies(coded)[[1]]
  kinshipCalcsNoInvar(withData, subset = subset, targetN = targetN)
}
