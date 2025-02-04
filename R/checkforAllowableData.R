#' Filter valid base pair combinations
#' @description This function is used to only select columns that contain allowable characters; specifically, the allele pairs listed in the `IUPAC` table.
#' @param dataset A data frame, where the first column named `Coral_ID` uniquely identifies the colony in each row, and the second through last column each contain allele data for a single locus.
#' @returns description
#' @export
checkforAllowableData <- function(dataset) {dataset %in% IUPAC$Allelepairs}
