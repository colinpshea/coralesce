#' Classify allele pairs
#' @description This function creates a data frame of all possible pairwise combinations of individuals and determines, for each  single locus, whether two individuals have the same allele or not. This function only returns comparisons for colonies with other colonies and not with themselves (`classifyAllelePairs()` does both for the genet assignment procedure). This function works on one locus (i.e., column) at a time. 
#' @param dataset A data frame, where the first is named `Coral_ID` that uniquely identifies each row (colony), and the second through last column each contain allele data for a single locus.
#' @param locus The column number containing allele data for a single locus.
#' @returns This function returns a data frame with six columns: `coral1` (name of coral 1), `coral2` (name of coral 2), `allele1` (allele for coral 1 at that locus), `allele2` (allele for coral 2 at that locus), `locus` (the locus name), and match (are the `IUPAC` alleles the same: `TRUE` or `FALSE`). All possible pairwise comparisons others (not selves) are included here, and the resulting file is used by `kinshipCalcsNoInvar()` which calculates kinship and gene diversity for all colonies and subsets of colonies. 
#' @importFrom magrittr %>% %<>%
#' @export
classifyAllelePairsOthers <- function(dataset, locus){
  locus_data <- data.frame(
    dataset[,c(1,locus)]
  )
  possible_combos_others <- data.frame(
    Coral_ID = t(combn(dataset$Coral_ID,2))
  )
  possible_combos_others <- possible_combos_others %>%
    merge(locus_data, by.x="Coral_ID.1", by.y="Coral_ID") %>%
    mutate(locus=names(.)[3]) %>%
    merge(locus_data, by.x="Coral_ID.2", by.y="Coral_ID")
  names(possible_combos_others) <-c(
    "coral2",
    "coral1",
    "allele1",
    "locus",
    "allele2"
  )
  possible_combos_others %<>%
    select(coral1, coral2, allele1, allele2, locus) %>%
    mutate(
      match = ifelse(
        allele1==allele2,
        TRUE,
        FALSE
      )
    )
  return(possible_combos_others)
}
