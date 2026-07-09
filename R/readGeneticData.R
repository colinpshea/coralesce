#' Read and validate a genetic SNP data file
#'
#' @description Reads a genotype CSV of paired alleles (e.g., `C:G`, `A:T`) with
#'   `Coral_ID` and `MatchMaker_Index` as the first two columns. Values are read
#'   as character (so `T` is not treated as `TRUE`), internal whitespace is
#'   stripped, the input-column contract is enforced by
#'   [handleError_ColumnContract()], and `MatchMaker_Index` is converted to
#'   integer. Missing SNP values are recorded as `"?"` so unscored allele combos
#'   can be tracked through the `IUPAC` translation.
#' @param fileloc Path to the CSV data file. When calling this function outside
#'   [runGenets()] / [runKinship()], supply the full path yourself, e.g.
#'   `readGeneticData(file.path(getwd(), "Data", "myfile.csv"))`.
#' @returns A list of two data frames: `[[1]]` the cleaned genotype data
#'   (`Coral_ID`, `MatchMaker_Index`, and SNP columns), and `[[2]]` an index of
#'   `Coral_ID` and integer `MatchMaker_Index`.
#' @importFrom stringi stri_replace_all_regex
#' @importFrom dplyr select
#' @importFrom magrittr %>%
#' @importFrom utils read.csv
#' @export
readGeneticData <- function(fileloc) {
  raw <- read.csv(fileloc, stringsAsFactors = FALSE,
                  colClasses = "character", na.strings = c("", " "))

  # Strip stray whitespace (e.g., " C:T", "C: T") from every value.
  raw[] <- lapply(raw, stri_replace_all_regex, " ", "")

  # Enforce column contract BEFORE anything is renamed or coerced.
  handleError_ColumnContract(raw)
  raw$MatchMaker_Index <- as.integer(as.numeric(raw$MatchMaker_Index))

  # Record missing SNP values as "?" (SNP columns only; keys are already clean).
  snpCols <- setdiff(names(raw), c("Coral_ID", "MatchMaker_Index"))
  if (length(snpCols) > 0) {
    sub <- raw[snpCols]
    sub[is.na(sub)] <- "?"
    raw[snpCols] <- sub
  }

  handleError_ProhibitedData(raw, acceptableData = IUPAC)
  find_dups(raw)

  index <- raw %>% select(Coral_ID, MatchMaker_Index)
  list(raw, index)
}
