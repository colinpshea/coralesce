#' Isolate all NA colonies
#' @description This function determines whether any colonies have `NA` values at all loci, in which case they are isolated in a new data frame called `allNA`. The remaining data frame contains colonies with SNP data at, at the very least, one locus. 
#' @param dataset A data frame containing a `Coral_ID` column for each of N colonies and valid `IUPAC-allele` codes (i.e., G, not G:G) for SNP data at each of N loci. This function is used in conjunction with `convertBasePairstoCodes()` during initial processing. 
#' @returns This function returns two objects: (1) if such individuals are present, a data frame called `allNA` is returned with the identify of `allNA` colonies, and (2) a version of the input data frame that excludes any `allNA` colonies. If there are no `allNA` colonies, then the original data frame is returned and the `allNA` data frame is a `NULL` object. Additionally, if `allNA` individuals are present, their identifies will be listed in a warning message issues by the `handleError_allZeros()` function. 
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
