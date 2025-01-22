#' Read in Genetic Data into R and Format Properly
#'
#' @description With the fileloc specified, one can import the initial data that's got paired alleles (e.g., C:G, A:T, etc.) along with a Coral_ID columns and any other columns like site name. This function names the first column "Coral_ID", makes sure that R interprets "T" as a character rather than the logical "TRUE", which would be problematic. The function also omits any potential whitespace (e.g., " C:T"). Lastly, "NA" values are converted to "?" because we need to keep track of allele combos that couldn't be compared, and NA could be problematic when converting allele pairs to single letters with our IUPAC cipher.
#' 
#' This fileloc path is built into `runGets` as follows: 
#' folderPaths <- list.dirs(path = paste0(getwd()), full.names = TRUE, recursive = F)
#' dataLocation <- folderPaths[which(str_detect(folderPaths, "Data") == TRUE)]
#' resultsLocation <- folderPaths[which(str_detect(folderPaths, "Results") == TRUE)]
#' fileList <- as.list(list.files(path = dataLocation, pattern = "\\.csv$"))
#' for (i in 1:length(fileList)) {
#' a <- readGeneticData(fileloc = paste0(dataLocation, "/", fileList[[i]]))
#' @param fileloc The location of the data file, specified as a path. 
#' @importFrom stringi stri_replace_all_regex
#' @importFrom purrr map
#' @importFrom magrittr %>%
#' @export
readGeneticData <- function(fileloc) {
  raw <- read.csv(
    fileloc,
    stringsAsFactors = FALSE,
    colClasses = c("character"),
    na.strings = c("", " ")
  )
  raw <- map(raw, stri_replace_all_regex," ", "") %>%
    as.data.frame
  handleError_CoralID(raw)
  names(raw)[1] <- "Coral_ID"
  raw[is.na(raw)] <- "?"
  handleError_ProhibitedData(raw, acceptableData = IUPAC)
  return(raw)
}
