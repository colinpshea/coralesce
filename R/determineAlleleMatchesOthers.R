#' Assess matching of alleles, comparing only to other colonies.
#' @description For each pairwise combination of individuals, classify allele pairs at each locus as either matching or not. This works one column (locus) at a time, by looping through all loci. This function ignores matches of colonies with themselves i.e., it only compares to other colonies. 
#' @param dataset A data frame, where the first column titled "Coral_ID" uniquely identifies each row, and the second through last column each contain allele data for a single locus.
#' @export
determineAllAlleleMatchesOthers <- function(dataset){
  A <- list()
  for (j in 2:ncol(dataset)){
    A[[j]] <- classifyAllelePairsOthers(
      dataset = dataset,
      locus = j
    )
  }
  do.call(rbind.data.frame, A)
}
