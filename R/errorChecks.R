#' Notify user that `Coral_ID` field name had to be changed.
#' @description This function notifies users that the `Coral_ID` column name was changed. 
#' @param dataset A data frame with some sort of coral colony identifier (ideally named `Coral_ID`) and SNP genetic data (and perhaps other fields) in the remaining columns.
#' @returns This function returns a warning message notifying users that the colony identifier field was changed to `Coral_ID` prior to any genet assignment or kinship calculations. If the column is already called `Coral_ID` then nothing happens and no warning message is issued.  
#' @export
handleError_CoralID <- function(dataset) {
  if (!("Coral_ID" %in% names(dataset))) {message("The colony identifier field was manually renamed Coral_ID prior to genet assignment and/or kinship calculations.
                                                  ")}
}

#' Notify user that colonies with no valid SNP data were found. 
#' @description This function finds colonies with all `NA` values, report with a message, and separates the `allNA` individuals from the rest of the data. These `allNA` individuals are appended to the genet classification file and assigned `genet = XXXX_NA`, `pctNULL = 100`, and `AdequateData = No`.
#' @param dataset Takes a file with SNP data with a `Coral_ID` column as the first column and SNP data in the remaining columns. 
#' @returns This function returns a data frame of colonies with all NA values, which is then evaluated by `isolateAllNAColonies()`. If this is `NULL`, then nothing happens, but if a data frame with at least one observation is returned, then `isolateAllNAColonies()` will remove the offending individuals and append their data to the genet assignment table. If any `allNA` colonies are identified, this function also issues a warning message notifying users of the identity of the offending colonies. 
#' @export
handleError_allZeros <- function(dataset){
  testData <- dataset
  testData$allNA <- rowSums(is.na(testData[,2:ncol(testData)]))==length(2:ncol(testData))
  allNA <- testData %>% filter(allNA==TRUE) %>% mutate(genet = NA, AdequateData = "No") 
  if (nrow(allNA) > 0) {message("At least one colony has NA values at all loci. These colonies have been added to the genetAssignment data frame and have been assigned genet = XXXX_NA, pctNull = 100, and AdequateData = No.
                                ")
  }
  if (nrow(allNA) > 0) {message(cat("The offending colonies are:", allNA$Coral_ID, sep = "\n"))
  }
  return(allNA)
}

#' Notify user that non-conforming fields were found and omitted. 
#' @description This function notifies users that non-conforming fields/columns such as site name or, more importantly, allele pairs that are not included in the `IUPAC` reference file, were omitted from the data set. This function omits the offending columns and reports a message that they were removed along with their names.
#' @param dataset A data frame to be tested with a `Coral_ID` column and SNP data, along with other fields such as site name etc. 
#' @param acceptableData A data frame with acceptable data that is allowable; here, these are `IUPAC` SNP data. 
#' @returns This function returns a warning message notifying users of the presence and identity of columns in their data that do not adhere to formatting requirements. This is used in conjunction with `checkforAllowableData()`, which uses this function and also omits the offending columns.
#' @export
handleError_ProhibitedData <- function(dataset, acceptableData) {
  if (sum(colSums(apply(dataset[,2:ncol(dataset)], 2, checkforAllowableData)) < nrow(dataset)) > 0) {
    message("Columns other than Coral_ID that do not adhere to the required base pair format (e.g., a site name column or an invalid base pair) were removed prior to genet assignment and/or kinship calculations.
            ")
    message(cat("The offending columns are:", names(dataset[2:ncol(dataset)])[colSums(apply(dataset[2:ncol(dataset)], 2, checkforAllowableData)) != nrow(dataset)], sep = "\n"))
  }
}
