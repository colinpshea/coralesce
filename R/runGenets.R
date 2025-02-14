#' Calculate genet assignments and pariwise comparisons in a single function
#'
#' @description This wrapper function assigns colonies to genets and (optionally) summarizes pairwise comparisons of alleles for all possible pairwise comparisons. Note that for this function to work properly (or at all), you **MUST** have folders called `Data` and `Results` in your working directory. This function will look for a genetics data file in `Data` and save results to the `Results` folder. To read in the data, this function passes the `Data` path to `readGeneticData()` as the `fileloc` argument. If there is more than one file in `Data`, this function will cycle through each file and save results for each data file to the `Results` folder; the file name for each file in the `Results` folder will include the name of the original data file. 
#' 
#' @param PctMatchThreshold The desired threshold for percent match of alleles across all loci between individuals for identifying pairings as matches or clones.
#' @param PctNotNullThreshold The desired threshold for the percentage of valid SNP data across all loci. Given two individuals are determined to be a match (i.e., percent match of alleles ≥ `PctMatchThreshold`), `PctNotNullThreshold` is the minimum allowable percentage of loci with valid (i.e., not `NULL` or `NA`) allele data for making a genet assignment. Individuals with values `≥ PctNotNullThreshold` that are determined to be matches (with themselves or others) are classified as `AdequateData = No` and appended to the genet assignment file with `genet = XXXX_NA`, where `XXXX` is the 4-letter species code.
#' @param getPairwiseAlleleMatches Set to `TRUE` if you want to return a data frame with all pairwise comparisons and their corresponding percent match and percent not null values. The default value is `FALSE`.
#' @returns This function returns up to two objects depending on user inputs: 
#' 
#' The first object, `genetAssignment` is a data frame with a single row for each colony along with their genet number, percent null values across all of their loci, and whether or not the data were adequate for assigning them to a genet. Data adequacy is defined by the user-defined `PctMatchThreshold` and `PctNotNullThreshold` values. 
#' 
#' The second object, `pairwiseAlleleMatches`, containing ALL possible pairwise comparisons (each colony with itself and other colonies) at each locus, calculating percent match and percent not null. These are typically very large files and are only saved to the working directory if `getPairwiseAlleleMatches = TRUE`. 
#'  
#' @importFrom stringr str_detect str_pad
#' @export
runGenets <- function(PctMatchThreshold = NULL, PctNotNullThreshold = NULL, getPairwiseAlleleMatches = FALSE){
  #### Determine folder paths - just set the working directory to the right place and this will work fine: we're just looking for the names of all the folders in the working directory here. 
  folderPaths <- list.dirs(path = paste0(getwd()), full.names = TRUE, recursive = F)
  
  #### Specify locations of data and results folders 
  dataLocation <- folderPaths[which(str_detect(folderPaths, "Data")==TRUE)]
  resultsLocation <- folderPaths[which(str_detect(folderPaths, "Results")==TRUE)]
  
  #### Create a list of file names to be processed for genet assignment/kinship etc. If there is more than one file, they will each be processed separately, and if there's only one file this does nothing special.  
  fileList <- as.list(list.files(path = dataLocation, pattern = "\\.csv$"))
  
  #### Loop through all available data files
  for (i in 1:length(fileList)){
    a <- readGeneticData(fileloc = paste0(dataLocation,"/", fileList[[i]]))
    b <- isolateAllNAColonies(convertBasePairstoCodes(initdata = a))
    b1 <- b[[1]] # data frame with colonies that DO NOT have NA values at all loci
    b2 <- b[[2]] # data frame with colonies that DO have NA values at all loci
    c <- determineAllAlleleMatches(dataset = b1)
    d1 <- groupByGenets(CoralAlleleData = b1, AlleleMatchResults = c, PctMatchThreshold = PctMatchThreshold, PctNotNullThreshold = PctNotNullThreshold, getPairwiseAlleleMatches = getPairwiseAlleleMatches)
    d2 <- d1$genetAssignment %>% add_row(b2) %>% arrange(genet, Coral_ID) %>% mutate(genet = paste0(substr(fileList[[i]], start = 1, stop = 4), "_", str_pad(genet, 5, side = "left", pad = 0)))
    write.csv(d2, paste0(resultsLocation,"/","genetAssignment_", paste0(fileList[[i]])), row.names = F)
    if (getPairwiseAlleleMatches==TRUE){write.csv(d1$pairwiseAlleleMatches, paste0(resultsLocation,"/","pairwiseAlleleMatches_", paste0(fileList[[i]])), row.names = F)
    }
  }
  if (getPairwiseAlleleMatches==TRUE){ return(list(genetAssignments = d2, pairwiseAlleleMatches = d1))}
  if (getPairwiseAlleleMatches==FALSE){ return(list(genetAssignments = d2))}
}
