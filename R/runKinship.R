#' Run all kinship and  gene diversity calculations in a single function
#'
#' @description This function calculates pairwise kinship across all individuals and loci. Average kinship is calculated at the individual and population level, and both calculations exclude invariant loci. 
#' @param subset Do you want to subset the data and calculate kinship for the `targetN` least-related colonies? Default is FALSE; if changed to TRUE, then you must specify `targetN` individuals (i.e., some value less than or equal to the total number of colonies in the data set).
#' @param targetN The desired number of individuals over which population-average kinship and gene diversity are calculated. This is ignored if `subset = FALSE`. If `subset = TRUE` then `targetN` must specified and can be ≤ the number of individuals in a data set. If `targetN` is equal to the number of individuals in the data set, then MK_init and MK_final will be identical; if this value is smaller than the number of individuals in the data set, MK_final will differ because kinship is recalculated repeatedly, removing the individual with the highest average kinship, one-by-one, until `targetN` individuals remain.
#' @return This function returns up to three objects depending on user inputs: 
#' 
#' The first object, `PopAvgMKGD` is a data frame with a single row and two values, mean population-level kinship and mean population-level gene diversity (1 - mean population level kinship). 
#' 
#' The second object, `MK_init`, is a data frame with a row for each coral colony and an average kinship column that's an average of pairwise kinship for that individual across all other individuals and all loci.
#' 
#' The third object, `MK_final` is simular to `MK_init` except it only has targetN rows. 
#'  
#' @importFrom stringr str_detect
#' @export
runKinship <- function(subset = FALSE, targetN = NULL){
  #### Determine folder paths - just set the working directory to the right place and this will work fine: we're just looking for the names of all the folders in the working directory here. 
  folderPaths <- list.dirs(path = paste0(getwd()), full.names = TRUE, recursive = F)
  
  #### Specify locations of data and results folders
  dataLocation <- folderPaths[which(str_detect(folderPaths, "Data")==TRUE)]
  resultsLocation <- folderPaths[which(str_detect(folderPaths, "Results")==TRUE)]
  
  #### Create a list of file names to be processed for genet assignment/kinship etc. If there is more than one file, they will each be processed separately.  
  fileList <- as.list(list.files(path = dataLocation, pattern = "\\.csv$"))
  
  #### Loop through all available data files
  for (i in 1:length(fileList)){
    a <- readGeneticData(fileloc = paste0(dataLocation,"/", fileList[[i]])) 
    b <- isolateAllNAColonies(convertBasePairstoCodes(initdata = a))[[1]]
    c <- omitInvariantLoci(b)
    d <- kinshipCalcsNoInvar(dataset = c, targetN = targetN, subset = subset)
    if (subset==FALSE){
    write.csv(d$PopAvgMKGD, paste0(resultsLocation,"/","popAvgMKGD_", nrow(b), "_", paste0(fileList[[i]])), row.names = F)
    write.csv(d$MK_init, paste0(resultsLocation,"/","kinshipInit_", nrow(b), "_", paste0(fileList[[i]])), row.names = F)
    }
    if (subset==TRUE){
    write.csv(d$PopAvgMKGD, paste0(resultsLocation,"/","popAvgMKGD_", nrow(b), "_", paste0(fileList[[i]])), row.names = F)
    write.csv(d$MK_init, paste0(resultsLocation,"/","kinshipInit_", nrow(b), "_", paste0(fileList[[i]])), row.names = F)
    write.csv(d$MK_final, paste0(resultsLocation,"/","kinshipTargetN_", targetN, "_", paste0(fileList[[i]])), row.names = F)
    }
  }
  return(list(kinshipCalculations = d))
}
