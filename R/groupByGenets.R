#' Determine if individuals come from the same genet.
#' @description This function compares alleles at loci among colonies and calculates the percentage match and percentage not null and uses this information to determine if the pairwise comparison indicates that the two individuals are from the same genet. Uses `returnGenetIdentity()` for the genet classification.
#' 
#' @param CoralAlleleData Raw coral data (rows = coral colonies and columns = SNP data at loci; with one-letter codes for allele combinations and `NA` for `?` or blank values). This is used to append the colony-specific `pctNull` values to the genet assignment file created by this function.
#' @param AlleleMatchResults data set with results of pairwise matching of alleles at each locus for all individuals and loci
#' @param PctMatchThreshold Defaults to `NULL` but must be specified (can be specified in `runGenets()` wrapper function)
#' @param PctNotNullThreshold Defaults to `NULL` but must be specified (can be specified in `runGenets()` wrapper function)
#' @param getPairwiseAlleleMatches Do you want the function to return all pairwise comparisons of colonies? The default is `FALSE`.
#' @returns This function returns up to two objects: one is a data frame with final genet assignments, and the other is an optional data frame with all possible pairwise allele matches (returned when `getPairwiseAlleleMatches = TRUE`).
#' @importFrom dplyr arrange if_else n rename add_row distinct
#' @export
groupByGenets <- function(CoralAlleleData, AlleleMatchResults, PctMatchThreshold = NULL, PctNotNullThreshold = NULL, getPairwiseAlleleMatches = FALSE) {
  
  CoralAlleleData$pctNull <- 100 - apply(subset(CoralAlleleData, select = -Coral_ID), 1, calcPercentNotNull)
  
  CoralAlleleData <- CoralAlleleData %>% select(Coral_ID, pctNull)
  temp <- AlleleMatchResults %>% mutate(CoralPair = interaction(coral1, coral2)) %>% select(CoralPair, coral1, coral2, locus, match) %>% arrange(CoralPair, locus) %>% pivot_wider(names_from = locus, values_from = match)
  
  temp$pctMatch = rowMeans(temp[, 4:(ncol(temp))], na.rm = T)*100
  
  temp$pctNotNull <- apply(subset(temp, select = -c(CoralPair, coral1, coral2, pctMatch)), 1, calcPercentNotNull)
  
  temp %<>% mutate(PartOfGenet = ifelse(pctMatch >= PctMatchThreshold, "Yes", "No")) %>%
    select(coral1, coral2, CoralPair, pctMatch, pctNotNull, PartOfGenet)
  
  PartOfGenet_No <- temp %>% filter(PartOfGenet == "No")
  
  PartOfGenet_Yes <- temp %>% filter(PartOfGenet == "Yes") %>% 
    mutate(flag = if_else(coral1 != coral2 & pctNotNull < PctNotNullThreshold, "drop", "keep")) %>%
    filter(flag == "keep") %>% select(coral1, coral2, CoralPair, pctMatch, pctNotNull, PartOfGenet) %>%
    mutate(AdequateData = if_else(coral1 == coral2 & pctNotNull < PctNotNullThreshold, "No", "Yes"))
  
  finalYesClonesAdequateYes <- PartOfGenet_Yes %>% filter(AdequateData == "Yes") %>% mutate(obs = 1:n())
  
  finalYesClonesAdequateNo <- PartOfGenet_Yes %>% filter(AdequateData == "No") %>%
    mutate(genet = NA, pctNull = 100 - pctNotNull) %>%
    select(coral1, genet, pctNull, AdequateData) %>%
    rename(Coral_ID = coral1)
  
  groupedGenets <- returnGenetIdentity(finalYesClonesAdequateYes)
  
  genetAssignment <- finalYesClonesAdequateYes %>%
    left_join(groupedGenets, by = "obs") %>% 
    select(coral1, coral2, genet, AdequateData) %>%
    pivot_longer(-c(genet, AdequateData), names_to = NULL, values_to = "Coral_ID") %>% 
    select(Coral_ID, genet, AdequateData) %>%
    distinct(.) %>%  
    arrange(genet) %>% 
    left_join(CoralAlleleData, by = "Coral_ID") %>%
    select(Coral_ID, genet, pctNull, AdequateData) %>% 
    add_row(finalYesClonesAdequateNo)
  if (getPairwiseAlleleMatches==TRUE) {return(list(genetAssignment = genetAssignment, pairwiseAlleleMatches = temp))}
  if (getPairwiseAlleleMatches==FALSE) {return(list(genetAssignment = genetAssignment, pairwiseAlleleMatches = NULL))}
  }
