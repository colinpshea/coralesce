#' For each coral pair observation, return the genet identity
#'
#' @description For each coral pair, return the identity of the genet to which the pair belongs..
#' @param obs a dataframe of three columns named coral1, coral2, and CoralPair. coral1 identifies the first coral in the pair, coral2 identifies the second coral in the pair, and CoralPair identifies the pair, taking the format of the first coral identity concatenated together with the second coral identity, separated by a period.
#' @examples
#' example <- data.frame(
#'   coral1=c(1,2,3,3,4,5,6),
#'   coral2=c(1,2,2,3,4,5,6),
#'   CoralPair=c("1.1","2.2","3.2","3.3","4.4","5.5","6.6")
#' )
#' returnGenetIdentity(example)
#' @importFrom igraph graph_from_adjacency_matrix components
#' @importFrom Matrix sparseMatrix tcrossprod
#' @export
returnGenetIdentity <- function(obs) {
  ## Create a list of names for each individual in the data set
  indList <- as.character(sort(unique(c(
    obs$coral1,
    obs$coral2
  ))))
  ## Create empty list for storing results
  grpList <- vector("list", length(indList))
  ## Create another empty list for storing results
  matchList <- list()
  ## Loop through each individual and determine all pairwise comparisons for which it was determined to be a clone
  for (i in 1:length(indList)){
    matchList <- grep(paste0("([[:punct:]]|^)", indList[i], "([[:punct:]]|$)"), obs$CoralPair)
    grpList[[i]] <- sort(c(matchList, grpList[[i]]))
  }
  ## Some network/graph analysis magic happens here: result is assigment of individuals to genets based on relatedness (i.e., all clones are identified and placed in unique groups).
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

#' Determine if the two individuals in the pairwise comparison come from the
#' same genet
#'
#'@param finalResultsWide the results to return in "wide" format
#'@param PctMatchThreshold Defaults to 99.
#'
#' @description The function calculates the percentage match and percentage not
#' null and uses this information to determine if the pairwise comparison
#' indicates that the two individuals are from the same genet.
#' @importFrom dplyr arrange if_else n rename add_row distinct
#' @export
groupByGenets <- function(AlleleMatchResults, PctMatchThreshold = NULL, PctNotNullThreshold = NULL) {
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
    mutate(genet = NA) %>%
    select(coral1, genet, AdequateData) %>%
    rename(Coral_ID = coral1)
  groupedGenets <- returnGenetIdentity(finalYesClonesAdequateYes)
  genetAssignment <- finalYesClonesAdequateYes %>%
    left_join(groupedGenets, by = "obs") %>% 
    select(coral1, coral2, genet, pctNotNull, AdequateData) %>%
    pivot_longer(-c(genet, AdequateData), names_to = NULL, values_to = "Coral_ID") %>% 
    select(Coral_ID, genet, AdequateData) %>% distinct(.) 
    arrange(genet) %>% add_row(finalYesClonesAdequateNo)
    return(genetAssignment)
}
