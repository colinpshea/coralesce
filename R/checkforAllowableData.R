#' Filter valid base pair combinations
#' @description This function is used to only select columns that contain allowable characters; specifically, the allele pairs listed in the `IUPAC` table.
#' @param dataset A data frame, where the first column named `Coral_ID` uniquely identifies the colony in each row, and the second through last column each contain allele data for a single locus. The dataset may also contain additional variables such as site name etc. that are ultimately removed; the use is notified of such columns being removed by `handleError_ProhibitedData()` in `readGeneticData()`.
#' @returns The supplied dataset with only SNP data and a column named `Coral_ID`. 
#' @export
checkforAllowableData <- function(dataset) {dataset %in% IUPAC$Allelepairs}
