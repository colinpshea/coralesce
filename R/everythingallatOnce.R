#' Run all genet assignments in a single function
#'
#' @description Assigns colonies to genets and calculates pairwise kinship across all individuals and loci. Average kinship is calculated at the individual and population level, and both calculations exclude invariant loci. Colonies for which data are inadequate (too many NULL observations for SNP loci or none a all) are classified as such and assigned to genet = NA; these colonies can be re-assigned to existing or new colonies should additional data become available.  
#' @param PctMatchThreshold The desired threshold for percent match of alleles across all loci between individuals for identifying pairings as matches or clones. 
#' @param PctNotNullThreshold The desired threshold for percent match of alleles across all loci between individuals for identifying pairings as matches or clones. 
#' @param targetN The desired number of individuals over which population-average kinship and gene diversity are calculated. This is ignored if left as NULL or if its value is greater than or equal to nrow(dataset) i.e., the number of individuals in a data set, nothing happens. If this value is smaller than nrow(), then kinship is recalculated repeatedly, removing the individual with the highest average kinship, one-by-one, until targetN individuals with the lowest averge kinship remain.
#' @importFrom stringr str_detect
#' @export
runGenets <- function(PctMatchThreshold = NULL, PctNotNullThreshold = NULL){
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
    d1 <- groupByGenets(AlleleMatchResults = c, PctMatchThreshold = PctMatchThreshold, PctNotNullThreshold = PctNotNullThreshold)
    d2 <- d1 %>% add_row(b2) %>% arrange(genet, Coral_ID)
    write.csv(d2, paste0(resultsLocation,"/","genetAssignment_", PctMatchThreshold,"_", PctNotNullThreshold, "_", nrow(b1), "_", paste0(fileList[[i]])), row.names = F)
  }
  return(list(genetAssignments = d2))
}


#' Run all kinship and  gene diversity calculations in a single function
#'
#' @description Calculates pairwise kinship across all individuals and loci. Average kinship is calculated at the individual and population level, and both calculations exclude invariant loci. Colonies for which data are inadequate (too many NULL observations for SNP loci or none a all) are classified as such and assigned to genet = NA; these colonies can be re-assigned to existing or new colonies should additional data become available.  
#' @param PctMatchThreshold The desired threshold for percent match of alleles across all loci between individuals for identifying pairings as matches or clones. 
#' @param PctNotNullThreshold The desired threshold for percent match of alleles across all loci between individuals for identifying pairings as matches or clones. 
#' @param targetN The desired number of individuals over which population-average kinship and gene diversity are calculated. This is ignored if left as NULL or if its value is greater than or equal to nrow(dataset) i.e., the number of individuals in a data set, nothing happens. If this value is smaller than nrow(), then kinship is recalculated repeatedly, removing the individual with the highest average kinship, one-by-one, until targetN individuals with the lowest averge kinship remain.
#' @importFrom stringr str_detect
#' @export
runKinship <- function(targetN = NULL){
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
    b <- isolateAllNAColonies(convertBasePairstoCodes(initdata = a))[[1]]
    c <- omitInvariantLoci(b)
    d <- kinshipCalcsNoInvar(dataset = c, targetN = targetN)
    write.csv(c[[1]], paste0(resultsLocation,"/","popAvgMKGD_", nrow(b), "_", paste0(fileList[[i]])), row.names = F)
    write.csv(c[[2]], paste0(resultsLocation,"/","mnKinshipALL_", nrow(b), "_", paste0(fileList[[i]])), row.names = F)
    write.csv(c[[3]], paste0(resultsLocation,"/","mnKinshipTargetN_", targetN, "_", paste0(fileList[[i]])), row.names = F)
  }
  return(list(kinshipCalculations = d))
}

