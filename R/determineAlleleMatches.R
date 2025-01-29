#' Assess matching of alleles, comparing to selves and other colonies.
#' @description For each pairwise combination of individuals, classify allele pairs at each locus as either matching or not. This works one column (locus) at a time, by looping through all loci. This function included matches of colonies with themselves and with other colonies. 
#' @param dataset A data frame, where the first column titled "Coral_ID" uniquely identifies each row, and the second through last column each contain allele data for a single locus.
#' @export
determineAllAlleleMatches <- function(dataset){
  A <- list()
  for (j in 2:ncol(dataset)){
    A[[j]] <- classifyAllelePairs(
      dataset = dataset,
      locus = j
    )
  }
  do.call(rbind.data.frame, A)
}
