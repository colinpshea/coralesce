#' Notify user that Coral_ID had to be changed. This is just a message and the presence of a name other than Coral_ID doesn't cause any problems. 
#' @export
handleError_CoralID <- function(dataset) {
  if (!("Coral_ID" %in% names(dataset))) {message("The colony identifier field was manually renamed Coral_ID prior to genet assignment and kinship/gene diversity calculations.
                                                  ")}
}

#' Find colonies with all NA values, report with a message, and separate those individuals from the rest of the data.  with at least some genetic data (adequate or not). These "all NA" individuals are appended to the genet classification file and have genet = NA, pctNULL = 100, and adequateData = No.
#' @description This function finds colonies with all NA values, report with a message, and separate those individuals from the rest of the data.  with at least some genetic data (adequate or not). These "all NA" individuals are appended to the genet classification file and have genet = NA, pctNULL = 100, and adequateData = No.
#' @export
handleError_allZeros <- function(dataset){
  testData <- dataset
  testData$allNA <- rowSums(is.na(testData[,2:ncol(testData)]))==length(2:ncol(testData))
  allNA <- testData %>% filter(allNA==TRUE) %>% mutate(genet = NA, AdequateData = "No") 
  if (nrow(allNA) > 0) {message("At least one colony has NA values at all loci. These colonies have been added to the genetAssignment data frame and have been genet = NA, pctNull = 100, and adequateData = `No`.
                                ")
  }
  if (nrow(allNA) > 0) {message(cat("The offending colonies are:", allNA$Coral_ID, sep = "\n"))
    }
}

#' Notify user that non-conforming fields such as site name or - more importantly - allele pairs that are not included in the IUPAC reference file, were omitted from the data set. This is just a message as extraneous columns (i.e., those that aren't either locus data or Coral_ID) are removed for the genet classification and kinship calculation; however, the input file should be checked carefully because invalid allele pairs will cause entire loci to be dropped from the file. 
#' @export
handleError_ProhibitedData <- function(dataset, acceptableData) {
  if (sum(colSums(apply(dataset[,2:ncol(dataset)], 2, checkforAllowableData)) < nrow(dataset)) > 0) {
    message("Columns other than Coral_ID that do not adhere to the required base pair format (e.g., a site name column or an invalid base pair) were removed from this file prior to genet assignment and/or kinship/gene diversity calculations.
            ")
    message(cat("The offending column(s) are:", names(dataset[2:ncol(dataset)])[colSums(apply(dataset[2:ncol(dataset)], 2, checkforAllowableData)) != nrow(dataset)], sep = "\n"))
  }
}
