#' Classify allele pairs at all loci for all pairwise combinations of individuals, but exclude comparisons with selves
#' @description For each pairwise combination of individuals, classify allele pairs at each locus as either matching or not. This works one column (locus) at a time, by looping through all loci
#' @param dataset A dataframe, where the first column titled "Coral_ID" uniquely identifies each row, and the second through last column each contain allele data for a single locus.
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
