#' Assign genets for every data file in the working directory
#'
#' @description Wrapper that reads each genotype CSV in a `Data` folder, assigns
#'   colonies to genets, and writes a `genetAssignment_<file>.csv` to a `Results`
#'   folder (optionally also a `pairwiseAlleleMatches_<file>.csv`). Requires
#'   `Data` and `Results` folders in the working directory. All CSV files in
#'   `Data` are processed; results are keyed by input file name.
#' @param PctMatchThreshold Minimum percent allele match across loci for two
#'   colonies to be called a match/clone. Required.
#' @param PctNotNullThreshold Minimum percent of loci with valid data required to
#'   trust a comparison. Required.
#' @param getPairwiseAlleleMatches Logical; if `TRUE`, also write and return the
#'   full pairwise comparison table for each file. Default `FALSE`.
#' @returns Invisibly, a list with `genetAssignments` (a named list of
#'   per-file genet-assignment data frames) and, when
#'   `getPairwiseAlleleMatches = TRUE`, `pairwiseAlleleMatches` (a named list of
#'   per-file pairwise tables). Each genet-assignment data frame has columns
#'   `Coral_ID`, `MatchMaker_Index`, `genet`, `pctNull`, `AdequateData`.
#' @importFrom dplyr bind_rows left_join arrange mutate select
#' @importFrom stringr str_detect str_pad
#' @importFrom magrittr %>%
#' @importFrom utils write.csv
#' @export
runGenets <- function(PctMatchThreshold = NULL, PctNotNullThreshold = NULL,
                      getPairwiseAlleleMatches = FALSE) {

  loc             <- locateDataResults()
  dataLocation    <- loc$data
  resultsLocation <- loc$results

  fileList <- list.files(path = dataLocation, pattern = "\\.csv$")
  if (length(fileList) == 0) {
    warning("No .csv files found in the Data folder: ", dataLocation)
    return(invisible(list(genetAssignments = list())))
  }

  genetResults    <- vector("list", length(fileList))
  pairwiseResults <- vector("list", length(fileList))
  names(genetResults)    <- fileList
  names(pairwiseResults) <- fileList

  for (f in fileList) {
    dat <- readGeneticData(fileloc = file.path(dataLocation, f))
    snp   <- dat[[1]]     # cleaned genotype data
    index <- dat[[2]]     # Coral_ID + MatchMaker_Index

    split   <- isolateAllNAColonies(convertBasePairstoCodes(initdata = snp))
    withData <- split[[1]]
    noData   <- split[[2]]

    matches <- determineAllAlleleMatches(dataset = withData)
    gg <- groupByGenets(
      CoralAlleleData          = withData,
      AlleleMatchResults       = matches,
      PctMatchThreshold        = PctMatchThreshold,
      PctNotNullThreshold      = PctNotNullThreshold,
      getPairwiseAlleleMatches = getPairwiseAlleleMatches
    )

    speciesCode <- substr(f, 1, 4)
    assignment <- gg$genetAssignment %>%
      bind_rows(noData) %>%
      left_join(index, by = "Coral_ID") %>%
      arrange(MatchMaker_Index) %>%
      mutate(genet = paste0(speciesCode, "_", str_pad(genet, 5, "left", "0"))) %>%
      select(Coral_ID, MatchMaker_Index, genet, pctNull, AdequateData)

    write.csv(assignment,
              file.path(resultsLocation, paste0("genetAssignment_", f)),
              row.names = FALSE)
    genetResults[[f]] <- assignment

    if (isTRUE(getPairwiseAlleleMatches)) {
      write.csv(gg$pairwiseAlleleMatches,
                file.path(resultsLocation, paste0("pairwiseAlleleMatches_", f)),
                row.names = FALSE)
      pairwiseResults[[f]] <- gg$pairwiseAlleleMatches
    }
  }

  out <- list(genetAssignments = genetResults)
  if (isTRUE(getPairwiseAlleleMatches)) out$pairwiseAlleleMatches <- pairwiseResults
  invisible(out)
}
