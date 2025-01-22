#' For each coral pair observation, return the genet identity
#'
#' @description For each coral pair, return the identity of the genet to which the pair belongs..
#' @param obs a dataframe of three columns named coral1, coral2, and CoralPair. coral1 identifies the first coral in the pair, coral2 identifies the second coral in the pair, and CoralPair identifies the pair, taking the format of the first coral identity concatenated together with the second coral identity, separated by a period.
#' @examples None yet
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
