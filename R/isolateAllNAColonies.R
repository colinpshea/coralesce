# isolate all NA colonies
#' Convert to long format and then convert base pairs to IUPAC codes
#' @description This function replaces the base pairs with IUPAC codes and
#' converts the data to wide format so that there is a column per code.
#' # Make sure that the data we're using are ONLY either the
#' Coral_ID OR one of the character strings included in OnlyDataAllowed. If
#' there are no such columns then nothing happens. There's probably a more
#' elegant way to do this but it works fine.
#' @param initdata A dataset containing a column "Coral_ID" with unique value per row, as well as at least one column containing base pair data with at least one value matching base pair data in IUPAC.
#' @importFrom dplyr select left_join
#' @importFrom tidyr pivot_longer pivot_wider
#' @export
isolateAllNAColonies <- function(dataset){
  allNA <- handleError_allZeros(dataset)
  if (is.null(allNA)==TRUE){return(list(dataset, NULL))}
  if (is.null(allNA)==FALSE){
    allNA <- allNA %>% mutate(pctNull = 100) %>% select(Coral_ID, genet, pctNull, AdequateData)
    dataset <- dataset %>% filter(!(Coral_ID %in% allNA$Coral_ID))
    return(list(dataset, allNA))
  }
}
