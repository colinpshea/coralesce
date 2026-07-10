#' Population- and individual-level kinship and gene diversity
#'
#' @description Calculates individual- and population-level mean kinship and gene
#'   diversity (`1 - kinship`) across all colonies and loci, excluding invariant
#'   loci and within-colony comparisons. Optionally restricts the population
#'   summary to the `targetN` least related colonies by iteratively removing the
#'   colony with the highest mean kinship until `targetN` remain.
#' @param dataset Wide single-letter allele data (`Coral_ID` plus one column per
#'   locus) from [convertBasePairstoCodes()].
#' @param subset Logical; if `TRUE`, also compute the `targetN`-colony summary.
#'   Requires `targetN`. Default `FALSE`.
#' @param targetN Number of least-related colonies to retain when
#'   `subset = TRUE`. Must be >= 2.
#' @returns A list with `PopAvgMKGD` (population mean kinship and gene diversity),
#'   `MK_init` (individual mean kinship for all colonies), and `MK_final`
#'   (individual mean kinship for the retained colonies, or `NULL` when
#'   `subset = FALSE`).
#' @seealso [runKinship()] for notes on how clonal replicates (ramets) are
#'   handled and what `PopAvgGD` represents.
#' @importFrom dplyr select group_by summarise arrange desc slice pull filter
#' @importFrom tidyr pivot_longer
#' @importFrom magrittr %>%
#' @export
kinshipCalcsNoInvar <- function(dataset, subset = FALSE, targetN = NULL) {

  if (!subset && !is.null(targetN)) {
    stop("subset = FALSE but a targetN value was supplied. Set subset = TRUE ",
         "or leave targetN = NULL.", call. = FALSE)
  }
  if (subset && is.null(targetN)) {
    stop("subset = TRUE requires a targetN value.", call. = FALSE)
  }
  if (subset && targetN < 2) {
    stop("When subset = TRUE, targetN must be >= 2.", call. = FALSE)
  }

  dat1 <- omitInvariantLoci(dataset)
  kin  <- kinshipCalcs(determineAllAlleleMatchesOthers(dat1))

  MK_init <- indMeanKinship(kin)
  PopAvgMKGD <- MK_init %>%
    summarise(PopAvgMK = mean(ind_mean_kinship),
              PopAvgGD = 1 - mean(ind_mean_kinship))

  MK_final <- NULL
  if (subset) {
    highest <- character(0)                       # nothing removed yet
    while (nrow(dat1) > targetN) {
      dat1    <- dat1 %>% filter(!(Coral_ID %in% highest))
      kin     <- kinshipCalcs(determineAllAlleleMatchesOthers(dat1))
      highest <- indMeanKinship(kin) %>% slice(1) %>% pull(Coral_ID)
    }
    MK_final <- indMeanKinship(kin)
  }

  list(PopAvgMKGD = PopAvgMKGD, MK_init = MK_init, MK_final = MK_final)
}

# Internal helper (not exported): individual mean kinship, one row per colony,
# highest first. Imports are declared on kinshipCalcsNoInvar() above.
indMeanKinship <- function(kin) {
  kin %>%
    pivot_longer(cols = c(coral1, coral2), values_to = "Coral_ID") %>%
    select(Coral_ID, avg_kinship) %>%
    group_by(Coral_ID) %>%
    summarise(ind_mean_kinship = mean(avg_kinship), .groups = "drop") %>%
    arrange(desc(ind_mean_kinship))
}
