#' Separate colonies with no valid SNP data from the rest
#'
#' @description Splits the input into (1) colonies that have at least one scored
#'   locus and (2) colonies that are `NA` at every locus. The latter are
#'   identified by [handleError_allZeros()] and returned with `pctNull = 100` for
#'   later appending to the genet-assignment output.
#' @param dataset A wide data frame with a unique `Coral_ID` in column 1 and
#'   single-letter `IUPAC` allele codes in the remaining columns (as produced by
#'   [convertBasePairstoCodes()]).
#' @returns A list of two elements: `[[1]]` the data with all-`NA` colonies
#'   removed, and `[[2]]` a data frame of the all-`NA` colonies (columns
#'   `Coral_ID`, `genet`, `pctNull`, `AdequateData`), or `NULL` if there are none.
#' @importFrom dplyr mutate select filter
#' @importFrom magrittr %>%
#' @export
isolateAllNAColonies <- function(dataset) {
  allNA <- handleError_allZeros(dataset)
  if (is.null(allNA)) return(list(dataset, NULL))

  allNA <- allNA %>%
    mutate(pctNull = 100) %>%
    select(Coral_ID, genet, pctNull, AdequateData)

  dataset <- dataset %>% filter(!(Coral_ID %in% allNA$Coral_ID))
  list(dataset, allNA)
}
