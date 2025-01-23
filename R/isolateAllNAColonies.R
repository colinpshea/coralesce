#' Isolate all NA colonies
#' @description This function determines whether any colonies have NA values at all loci, in which case they are isolated in a new data frame called allNA. The remaining data frame contains colonies with at an allele at at least one locus. 
#' @param dataset A dataset containing a column "Coral_ID" base pair data at each of N loci.
#' @returns A data frame called allNA and version of the input data frame that excludes the allNA colonies. If there are no allNA colonies then the original data frame is returned and the allNA data frame is set to NULL. 
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
