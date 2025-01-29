#' Classify allele pairs at a single locus for all pairwise combinations of individuals, but excluding pairings with selves
#' @description Create a data frame of all possible pairwise combinations of individuals and determine for a single locus whether each pair of individuals have the same allele or not. This function is only for comparisons with other colonies, not comparisons with others AND themelves (that's what classifyAllelePairs is for).
#' @param dataset A data frame, where the first column titled "Coral_ID" uniquely identifies each row, and the second through last column each contain allele data for a single locus.
#' @param locus The column number containing allele data for a single locus
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
