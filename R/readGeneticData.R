#' Read in genetic SNP data format for analysis
#'
#' @description With the `fileloc` specified, one can import the initial data that's got paired alleles (e.g., C:G, A:T, etc.) along with a `Coral_ID` columns and any other columns like site name. This function names the first column `Coral_ID` and ensures that R interprets `T` as a character rather than the logical `TRUE`. The function also omits any potential white space (e.g., " C:T" or "C: T"). Lastly, `NA` values are converted to `?` because we need to keep track of allele combos that couldn't be compared, and `NA` could be problematic when converting allele pairs to single letters using the `IUPAC` table.
#' @param fileloc The location of the data file, specified as a path (see documentation for `runGenets()` and `runKinship()` for information about `fileloc`).
#' 
#' Note that if using this function on its own outside of the `runGenets()` or `runKinship()` wrapper functions, you must specify the path to the genetics data file manually e.g.,
#'  
#' `readGeneticData(fileloc = paste0(getwd(),"/Data/", paste0(list.files(paste0(getwd(),"/Data/")))))`
#' 
#'  if the genetics data are located (as they should be) in a working directory folder called `Data`.
#' 
#' @returns This function returns a data frame with only: `Coral_ID` as column 1 and the remaining columns as loci with valid SNP data (here, valid SNP data are those allele combinations listed in the `IUPAC` table); columns with prohibited data are removed via `handleError_ProhibitedData()`, which identifies and removes such columns, if present, and reports their names in a warning message.
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
