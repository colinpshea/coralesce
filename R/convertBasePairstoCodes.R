#' Convert base-pair genotypes to single-letter IUPAC codes
#'
#' @description Keeps only `Coral_ID` and columns whose values are all valid
#'   allele pairs (per the `IUPAC` table), translates the paired alleles (e.g.,
#'   `G:G`) to single-letter IUPAC codes (e.g., `G`), and returns the result in
#'   wide format (one column per locus). Missing markers (`"?"`) are converted
#'   to `NA`.
#' @param initdata A data frame with a `Coral_ID` column and one or more columns
#'   of paired-allele SNP data. Non-conforming columns (e.g., site name, an
#'   invalid base pair, or the bookkeeping `MatchMaker_Index`) are dropped here;
#'   they are reported separately by [handleError_ProhibitedData()].
#' @returns A wide data frame with `Coral_ID` and one single-letter-code column
#'   per valid locus.
#' @importFrom dplyr select left_join
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom magrittr %>%
#' @export
convertBasePairstoCodes <- function(initdata) {
  # Keep Coral_ID plus any column whose values are all valid allele pairs.
  allowable <- vapply(initdata,
                      function(col) all(checkforAllowableData(col)),
                      logical(1))
  keep <- names(initdata)[allowable | names(initdata) == "Coral_ID"]
  raw  <- initdata[, keep, drop = FALSE]

  raw %>%
    pivot_longer(!Coral_ID, names_to = "Locus", values_to = "Allelepairs") %>%
    left_join(IUPAC, by = "Allelepairs") %>%
    select(-Allelepairs) %>%
    pivot_wider(names_from = "Locus", values_from = "IUPACallele") -> wideclean

  wideclean[wideclean == "?"] <- NA   # unscored alleles become NA for genet calcs
  wideclean
}
