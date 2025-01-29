#' Filter valid base pair combinations
#' @description This function is used to select only only columns that contain allowable characters (those listed in IUPAC)
#' @param dataset A data frame, where the first column titled "Coral_ID" uniquely identifies each row, and the second through last column each contain allele data for a single locus.
#' @export
checkforAllowableData <- function(dataset) {dataset %in% IUPAC$Allelepairs}
