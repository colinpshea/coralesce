#' Assess matching of alleles, comparing only to other colonies.
#' @description For each pairwise combination of individuals, this function classifies allele pairs at each locus as either matching or not. This function applies the classifyAllelePairs function, one column or locus at a time, by looping sequentially through all of `n:ncol(dataset)` loci. This function ignores allele matches of colonies with themselves i.e., it only compares and assess allele matches at an individual's loci to those of other colonies, and is specifically used for kinship calculations using `kinshipCalcsNoInvar()`, which requires only comparisons of colonies with other colonies, not themselves.
#' @param dataset A data frame, where the first column named `Coral_ID` that uniquely identifies each of `nrow(dataset)` colonies, and the second through last column each contain single-letter `IUPAC` allele data at `2:ncol(dataset)` loci.
#' @returns This function, like `determineAlleleMatches()` which uses `classifyAllelePairs()` for the genet assignment procedure, uses `classifyAllelePairsOthers()` to create and return a data frame with all possible comparisons of alleles at loci for all colonies; here, however, comparisons of colonies with themselves are excluded as those comparisons are not relevant for kinship calculations. 
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
