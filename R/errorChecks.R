#' Notify user that Coral_ID field name had to be changed.
#' @description This function notifies users that the Coral_ID column name was changed. 
#' @param dataset A data frame with some sort of coral colony identifier (ideally named Coral_ID) 
#' @export
handleError_CoralID <- function(dataset) {
  if (!("Coral_ID" %in% names(dataset))) {message("The colony identifier field was manually renamed Coral_ID prior to genet assignment and kinship calculations.
                                                  ")}
}

#' Notify user that colonies with no valid SNP data were found. 
#' @description This function finds colonies with all NA values, report with a message, and separates the `allNA` individuals from the rest of the data. These `all NA`individuals are appended to the genet classification file and assigned genet = NA, pctNULL = 100, and adequateData = No.
#' @param dataset Takes a file with SNP data with a Coral_ID column as the first column and SNP data in the remaining columns. 
#' @returns Returns a data frame of colonies with all NA values. If this is NULL, then nothing happens, but if this is a data frame with at least one observation, then isolateAllNAColonies() will remove the offending individuals and append to the genet assigment table. 
#' @export
handleError_allZeros <- function(dataset){
  testData <- dataset
  testData$allNA <- rowSums(is.na(testData[,2:ncol(testData)]))==length(2:ncol(testData))
  allNA <- testData %>% filter(allNA==TRUE) %>% mutate(genet = NA, AdequateData = "No") 
  if (nrow(allNA) > 0) {message("At least one colony has NA values at all loci. These colonies have been added to the genetAssignment data frame and have been assigned genet = NA, pctNull = 100, and adequateData = No.
                                ")
  }
  if (nrow(allNA) > 0) {message(cat("The offending colonies are:", allNA$Coral_ID, sep = "\n"))
  }
  return(allNA)
}

#' Notify user that non-conforming fields were found and omitted. 
#' @description This function notifies users that non-conforming fields/columns such as site name or, more importantly, allele pairs that are not included in the IUPAC reference file, were omitted from the data set. This function omits the offending columns and reports a message that they were removed along with their names.
#' @param dataset A data frame to be tested with a Coral_ID column and SNP data, along with other fields such as site name etc. 
#' @param acceptableData A data frame with acceptable data that is allowable; here, these are IUPAC SNP data. 
#' @export
handleError_ProhibitedData <- function(dataset, acceptableData) {
  if (sum(colSums(apply(dataset[,2:ncol(dataset)], 2, checkforAllowableData)) < nrow(dataset)) > 0) {
    message("Columns other than Coral_ID that do not adhere to the required base pair format (e.g., a site name column or an invalid base pair) were removed from this file prior to genet assignment and kinship calculations.
            ")
    message(cat("The offending column(s) are:", names(dataset[2:ncol(dataset)])[colSums(apply(dataset[2:ncol(dataset)], 2, checkforAllowableData)) != nrow(dataset)], sep = "\n"))
  }
}
