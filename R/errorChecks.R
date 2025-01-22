#' Notify user that Coral_ID had to be changed
#' @export
handleError_CoralID <- function(dataset) {
  if (!("Coral_ID" %in% names(dataset))) {message("The colony identifier field was manually renamed Coral_ID prior to genet assignment and kinship/gene diversity calculations")}
}

#' Notify user that non-conforming fields were omitted from the data set
#' @export
handleError_ProhibitedData <- function(dataset, acceptableData) {
  identifyProhibitedData <- function (x) {x %in% acceptableData$Allelepairs}
  if (length(names(dataset)[colSums(apply(dataset[,2:ncol(dataset)], 2, identifyProhibitedData))==0])>1) {
    message("Columns other than Coral_ID that do not adhere to the required base pair format (e.g., site name) were removed from this file prior to genet assignment and kinship/gene diversity calculations")}
 }

#' find colonies with all NA values and report a warning
#' @description We import our initial data that's got paired alleles (e.g., C:G,
#' A:T, etc.). This function names the first column "Coral_ID", makes sure that
#' R interprets "T" as a character rather than the logical "TRUE", which would be problematic.
#' The function also omits any potential whitespace (e.g., " C:T"). Lastly,
#' "NA" values are converted to "?" because we need to keep track of allele
#' combos that couldn't be compared, and NA could be problematic when
#' converting allele pairs to single letters with our cipher.
#' @export
handleError_allZeros <- function(dataset){
  testData <- dataset
  testData$allNA <- rowSums(is.na(testData[,2:ncol(testData)]))==length(2:ncol(testData))
  allNA <- testData %>% filter(allNA==TRUE) %>% mutate(genet = NA, AdequateData = "No") 
  if (nrow(allNA) > 0) {message("At least one colony has NA values at all loci and has been assigned genet = NA and adequateData = `No`")
    }
  return(allNA)
}
