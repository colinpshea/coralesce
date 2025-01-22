#' Filter valid base pair combinations
#' @description This function removes invalid base pair combinations and is used
#' by the convertBasePairstoAlleles function. DOESNT THIS JUST RETURN THE ALLELEPAIRS?
#' @param x the dataset to use
#' @export
filterOnlyDataAllowed <- function(x) {IUPAC$Allelepairs %in% x}
