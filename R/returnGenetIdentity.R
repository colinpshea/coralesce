#' Determine the genet identity of each coral colony based on pairwise comparisons of SNPs
#'
#' @description This function assigns coral colonies to genets based on all possible comparisons (i.e., comparing colonies to themselves and other colonies) genetic similarity of SNPs at each of N loci.
#' @param obs a data frame created by `groupByGenets()` that has eight columns: `coral1`, `coral2`, `CoralPair`, `pctMatch`, `pctNotNull`, `PartOfGenet`, `AdequateData`, and `obs`. `coral1` identifies the first coral in the pair, `coral2` identifies the second coral in the pair, and `CoralPair` identifies the pair, taking the format of the first coral identity concatenated together with the second coral identity, separated by a period. `pctMatch` and `pctNotNull` summarize the alleles for the comparison being made, `PartOfGenet` means `pctMatch` exceeds a threshold for identifying individuals as a match or clone, `AdequateData` characterizes the frequency of `NULL` or `NA` values, and `obs` is simply a sequential observation number that is used to group genets. 
#' @importFrom igraph graph_from_adjacency_matrix components
#' @importFrom Matrix sparseMatrix tcrossprod
#' @returns This function returns a data frame with each colony assigned to a numeric genet number. Multiple colonies can belong to the same genet, but all corals are compared to themselves; hence, each colony is classified as a clone of itself by default, provided it had adequate SNP data. and they therefore occupy their own genet and are unrelated to any other colonies. Colonies previously determined to have insufficient data for this classification are not included in this analysis and are instead added to the genet assignment data frame at a later step with `genet = XXXX_NA`, where `XXXX` is the 4-letter species code.
#' @export
returnGenetIdentity <- function(obs) {
  ## Create a list of names for each individual in the data set
  indList <- as.character(sort(unique(c(obs$coral1, obs$coral2))))
  
  ## Create empty list for storing results
  grpList <- vector("list", length(indList))
  
  ## Create another empty list for storing results
  matchList <- list()
  
  ## Loop through each individual and determine all pairwise comparisons for which it was determined to be a clone
  for (i in 1:length(indList)){
    matchList <- grep(paste0("([[:punct:]]|^)", indList[i], "([[:punct:]]|$)"), obs$CoralPair)
    grpList[[i]] <- sort(c(matchList, grpList[[i]]))
  }
  
  ## Some network/graph analysis magic happens here: result is assignment of individuals to genets based on relatedness (i.e., all clones are identified and placed in unique groups called genets)
  ii <- rep(1:length(grpList), lengths(grpList))
  jj <- factor(unlist(grpList))
  tab <- sparseMatrix(
    i = ii,
    j = as.integer(jj),
    x = TRUE,
    dimnames = list(NULL, levels(jj))
  )
  connects <- tcrossprod(tab, boolArith = TRUE)
  group <- components(graph_from_adjacency_matrix(as(connects, "lMatrix")))$membership
  results <- tapply(grpList, group, function(x) sort(unique(unlist(x))))
  grpResults <- list()
  for (i in 1:length(results)){
    grpResults[[i]] <- data.frame(
      obs = unlist(results[i]),
      genet = as.numeric(names(results)[i])
    )
  }
  ## Collapse list of data frames to a single data frame
  do.call(rbind.data.frame, grpResults) %>% arrange(obs)
}
