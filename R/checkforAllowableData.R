#' Filter valid base pair combinations
#' @description This function is used to select only only columns that contain allowable characters (those listed in IUPAC)
#' by the convertBasePairstoAlleles function. DOESNT THIS JUST RETURN THE ALLELEPAIRS?
#' @param x the dataset to use
#' @export
checkforAllowableData <- function(x) {x %in% IUPAC$Allelepairs}
