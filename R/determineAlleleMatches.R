#' Assess matching of alleles, comparing to selves and other colonies.
#' @description For each pairwise combination of individuals, this function classifies allele pairs at each locus as either matching or not. This function applies the classifyAllelePairs function, one column or locus at a time, by looping sequentially through all of N loci. This function includes allele matches at each locus for colonies compared to themselves and compared to other colonies. 
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
