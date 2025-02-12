#' Assess matching of alleles, comparing to selves and other colonies.
#' @description For each pairwise combination of individuals, this function classifies allele pairs at each locus as either matching or not. This function applies `classifyAllelePairs()`, one column or locus at a time, by looping sequentially through all of `2:ncol(dataset)` loci. This function includes allele matches at each locus for colonies compared to themselves and compared to other colonies and is specifically used for genet assignment, which requires both types of comparisons. 
#' @param dataset A data frame with `Coral_ID` as the first column and single-letter `IUPAC` allele data for `nrow(dataset)` colonies at each of `2:ncol(dataset)` loci. 
#' @returns This function returns a data frame with six columns: `coral1` (name of coral 1), `coral2` (name of coral 2), `allele1` (allele for coral 1 at that locus), `allele2` (allele for coral 2 at that locus), `locus` (the locus name), and match (are the `IUPAC` alleles the same: `TRUE` or `FALSE`). All possible pairwise comparisons with selves and others are included here. The resulting file is used by `groupByGenets()`, which assigns colonies to genets. 
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
