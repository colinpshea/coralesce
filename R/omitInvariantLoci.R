#' Calculate pairwise kinship for all pairwise comparisons across all individuals (only others) and loci (only invariant)
#'
#' @description Calculates pairwise kinship across all individuals and loci and but excludes invariant loci
#' @param dataset uses result of convertBasePairstoCodes, which is the product of this code is used in classifyAllelePairsOthers and then determineAllAlleleMatchesOthers
#' @importFrom matrixStats rowProds
#' @return the supplied data frame with invariant loci (i.e., columns removed)
#' @export
omitInvariantLoci <- function(dataset){
locus_names <- colnames(dataset[-1])
variant_invariant_loci <- as.data.frame(t(dataset)[-1,]) %>% rowwise() %>% mutate(numAlleles =  n_distinct(c_across(everything()), na.rm = TRUE)) %>% ungroup() %>% mutate(locus = locus_names) %>% select(locus, numAlleles) 
invariant_loci <- variant_invariant_loci %>% filter(numAlleles==1) %>% pull(locus)
#datNoInvar <- dater %>% select(-all_of(invariant_loci))
dataset %>% select(-all_of(invariant_loci))
}
