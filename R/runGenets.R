#' Run all genet assignments in a single function
#'
#' @description Assigns colonies to genets and calculates pairwise kinship across all individuals and loci. Average kinship is calculated at the individual and population level, and both calculations exclude invariant loci. Colonies for which data are inadequate (too many NULL observations for SNP loci or none a all) are classified as such and assigned to genet = NA; these colonies can be re-assigned to existing or new colonies should additional data become available.  
#' @param PctMatchThreshold The desired threshold for percent match of alleles across all loci between individuals for identifying pairings as matches or clones.
#' @param PctNotNullThreshold The desired threshold for percent match of alleles across all loci between individuals for identifying pairings as matches or clones. 
#' @param getPairwiseAlleleMatches Do you want to return a data frame with all pairwise comparisons 
#' @importFrom stringr str_detect str_pad
#' @export
runGenets <- function(PctMatchThreshold = NULL, PctNotNullThreshold = NULL, getPairwiseAlleleMatches = FALSE){
  #### Determine folder paths - just set the working directory to the right place and this will work fine: we're just looking for the names of all the folders in the working directory here. 
  folderPaths <- list.dirs(path = paste0(getwd()), full.names = TRUE, recursive = F)
  
  #### Specify locations of data and results folders: changed this is bit so we're not hard-wiring folderPaths[1] and folderPaths[2] but instead matching "Data" and "Results" strings. This is case sensitive so it's got to be Data and Results. Note that you don't want to have any other folders with "Data" or "Results" in the name in your working directory because it will just muddle things up (and the script will crash). 
  dataLocation <- folderPaths[which(str_detect(folderPaths, "Data")==TRUE)]
  resultsLocation <- folderPaths[which(str_detect(folderPaths, "Results")==TRUE)]
  
  #### Create a list of file names to be processed for genet assignment/kinship etc. If there is more than one file, they will each be processed separately.  
  fileList <- as.list(list.files(path = dataLocation, pattern = "\\.csv$"))
  
  #### Loop through all available data files
  for (i in 1:length(fileList)){
    a <- readGeneticData(fileloc = paste0(dataLocation,"/", fileList[[i]])) 
    b1 <- isolateAllNAColonies(convertBasePairstoCodes(initdata = a))[[1]]
    b2 <- isolateAllNAColonies(convertBasePairstoCodes(initdata = a))[[2]]
    c <- determineAllAlleleMatches(dataset = b1)
    d1 <- groupByGenets(CoralAlleleData = b1, AlleleMatchResults = c, PctMatchThreshold = PctMatchThreshold, PctNotNullThreshold = PctNotNullThreshold, getPairwiseAlleleMatches = getPairwiseAlleleMatches)
    d2 <- d1$genetAssignment %>% add_row(b2) %>% arrange(genet, Coral_ID) %>% mutate(genet = paste0(substr(fileList[[i]], start = 1, stop = 4), "_", str_pad(genet, 5, side = "left", pad = 0)))
    write.csv(d2, paste0(resultsLocation,"/","genetAssignment_", PctMatchThreshold,"_", PctNotNullThreshold, "_", nrow(b1), "_", paste0(fileList[[i]])), row.names = F)
    if (getPairwiseAlleleMatches==TRUE){write.csv(d1$pairwiseAlleleMatches, paste0(resultsLocation,"/","pairwiseAlleleMatches_", PctMatchThreshold,"_", PctNotNullThreshold, "_", nrow(b1), "_", paste0(fileList[[i]])), row.names = F)
    }
  }
  if (getPairwiseAlleleMatches==TRUE){ return(list(genetAssignments = d2, pairwiseAlleleMatches = d1))}
  if (getPairwiseAlleleMatches==FALSE){ return(list(genetAssignments = d2))}
}
