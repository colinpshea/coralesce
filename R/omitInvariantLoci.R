#' Identify and omit invariant loci for subsequent kinship calculations.
#'
#' @description This function omits invariant loci for subsequent population- and individual-level kinship calculations.
#' @param dataset A data frame with `1:nrow(dataset)` individual colonies as rows and loci with SNP data in `2:ncol(dataset)` columns. 
#' @importFrom matrixStats rowProds
#' @returns This function returns the supplied data frame with invariant loci (i.e., columns) removed. This is done for subsequent kinship calculations that must exclude invariant loci (see `runKinship()` and `kinshipCalcsNoInvar()`). 
#' @export
omitInvariantLoci <- function(dataset){
  locus_names <- colnames(dataset[-1])
  variant_invariant_loci <- as.data.frame(t(dataset)[-1,]) %>% rowwise() %>% mutate(numAlleles =  n_distinct(c_across(everything()), na.rm = TRUE)) %>% ungroup() %>% mutate(locus = locus_names) %>% select(locus, numAlleles) 
  invariant_loci <- variant_invariant_loci %>% filter(numAlleles==1) %>% pull(locus)
  dataset %>% select(-all_of(invariant_loci))
}
