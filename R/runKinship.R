#' Run all kinship calculations in a single function
#'
#' @description This wrapper function calculates pairwise kinship across all individuals and loci. Average kinship is calculated at the individual and population level, and both calculations exclude invariant loci. For this function to work properly, you **MUST** have folders called `Data` and `Results` in your working directory. This function will look for a genetics data file in `Data` and save results to the `Results` folder. To read in the data, this function passes the `Data` path to `readGeneticData()` as the `fileloc` argument. If there is more than one file in `Data`, the function will cycle through each file and save results to the `Results` folder; the file name for each file in the `Results` folder will include the name of the original data file. 
#' @param subset Do you want to subset the data and calculate kinship for the `targetN` least related colonies in the data set? The default value is FALSE, in which case `targetN` is ignored and defaults to NULL. If `subset = TRUE` then you MUST enter a value for `targetN`; otherwise the function will quit and print an error message.
#' @param targetN The desired number of individuals over which population-average kinship and gene diversity are calculated. This is ignored if left as `NULL` or if its value is greater than or equal to `nrow(dataset)` i.e., the number of individuals in a data set. If this value is `< nrow(dataset)`, then kinship is recalculated, removing the individual with the highest average kinship, one individual at a time, until `targetN` individuals with the lowest average kinship remain.
#' @returns This function returns up to three objects depending on user inputs: 
#' 
#' The first object, `PopAvgMKGD` is a data frame with a single row and two values, population-level mean kinship and gene diversity (1 - population-level mean kinship). 
#' 
#' The second object, `MK_init`, is a data frame with a row for each coral colony and a kinship column that's an average of pairwise kinship for that individual across all other individuals and all loci (i.e., individual-level mean kinship).
#' 
#' The third object, `MK_final` is simular to `MK_init` except it only has `targetN` rows/individuals. 
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
    write.csv(d$PopAvgMKGD, paste0(resultsLocation,"/","popAvgMKGD_", paste0(fileList[[i]])), row.names = F)
    write.csv(d$MK_init, paste0(resultsLocation,"/","kinship_Init_", paste0(fileList[[i]])), row.names = F)
    }
    if (subset==TRUE){
    write.csv(d$PopAvgMKGD, paste0(resultsLocation,"/","popAvgMKGD_", paste0(fileList[[i]])), row.names = F)
    write.csv(d$MK_init, paste0(resultsLocation,"/","kinship_Init_", paste0(fileList[[i]])), row.names = F)
    write.csv(d$MK_final, paste0(resultsLocation,"/","kinship_targetN_", paste0(fileList[[i]])), row.names = F)
    }
  }
  return(list(kinshipCalculations = d))
}
