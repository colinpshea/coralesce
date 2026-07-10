#' Run kinship calculations for every data file in the working directory
#'
#' @description Wrapper that reads each genotype CSV in a `Data` folder and
#'   computes individual- and population-level mean kinship and gene diversity
#'   (excluding invariant loci), writing results to a `Results` folder. Requires
#'   `Data` and `Results` folders in the working directory. All CSV files in
#'   `Data` are processed; results are keyed by input file name.
#' @param subset Logical; if `TRUE`, also compute the `targetN`-colony summary.
#'   Requires `targetN`. Default `FALSE`.
#' @param targetN Number of least-related colonies to retain when
#'   `subset = TRUE`. Must be >= 2.
#' @returns Invisibly, a list with `PopAvgMKGD` and `kinship_init` (each a named
#'   list of per-file data frames) and, when `subset = TRUE`, `kinship_targetN`.
#'   `kinship_*` data frames have columns `Coral_ID`, `MatchMaker_Index`,
#'   `ind_mean_kinship`.
#' @details
#' Kinship and gene diversity are computed across all colonies in the input,
#' including multiple ramets (physical fragments) of the same genet. This is
#' intentional: every colony is treated as an individually eligible candidate
#' (e.g., for breeding or outplanting), so clonal replicates are retained rather
#' than collapsed to one representative per genet. Consequently, `PopAvgGD`
#' describes the genetic diversity of the *eligible colony pool*, not a
#' clone-corrected, genet-level population estimate. When `subset = TRUE`, the
#' routine preferentially sheds genetically redundant colonies but keeps a
#' representative of each clonal group. If a clone-corrected population estimate
#' is needed instead, collapse colonies to one representative per genet (e.g.,
#' using the genet assignments from [runGenets()]) before running kinship.
#' @importFrom dplyr left_join arrange select
#' @importFrom stringr str_detect
#' @importFrom magrittr %>%
#' @importFrom utils write.csv
#' @export
runKinship <- function(subset = FALSE, targetN = NULL) {

  loc             <- locateDataResults()
  dataLocation    <- loc$data
  resultsLocation <- loc$results

  fileList <- list.files(path = dataLocation, pattern = "\\.csv$")
  if (length(fileList) == 0) {
    warning("No .csv files found in the Data folder: ", dataLocation)
    return(invisible(list(PopAvgMKGD = list(), kinship_init = list())))
  }

  popResults    <- vector("list", length(fileList))
  initResults   <- vector("list", length(fileList))
  targetResults <- vector("list", length(fileList))
  names(popResults)    <- fileList
  names(initResults)   <- fileList
  names(targetResults) <- fileList

  for (f in fileList) {
    dat <- readGeneticData(fileloc = file.path(dataLocation, f))
    snp   <- dat[[1]]
    index <- dat[[2]]

    withData <- isolateAllNAColonies(convertBasePairstoCodes(initdata = snp))[[1]]
    reduced  <- omitInvariantLoci(withData)
    k <- kinshipCalcsNoInvar(dataset = reduced, targetN = targetN, subset = subset)

    initTbl <- k$MK_init %>%
      left_join(index, by = "Coral_ID") %>%
      arrange(MatchMaker_Index) %>%
      select(Coral_ID, MatchMaker_Index, ind_mean_kinship)

    write.csv(k$PopAvgMKGD,
              file.path(resultsLocation, paste0("popAvgMKGD_", f)), row.names = FALSE)
    write.csv(initTbl,
              file.path(resultsLocation, paste0("kinship_Init_", f)), row.names = FALSE)
    popResults[[f]]  <- k$PopAvgMKGD
    initResults[[f]] <- initTbl

    if (isTRUE(subset)) {
      finalTbl <- k$MK_final %>%
        left_join(index, by = "Coral_ID") %>%
        arrange(MatchMaker_Index) %>%
        select(Coral_ID, MatchMaker_Index, ind_mean_kinship)
      write.csv(finalTbl,
                file.path(resultsLocation, paste0("kinship_targetN_", f)),
                row.names = FALSE)
      targetResults[[f]] <- finalTbl
    }
  }

  out <- list(PopAvgMKGD = popResults, kinship_init = initResults)
  if (isTRUE(subset)) out$kinship_targetN <- targetResults
  invisible(out)
}

# Internal helper (not exported): locate the single Data and Results folders in
# the working directory, erroring clearly if they are missing or ambiguous.
locateDataResults <- function() {
  folderPaths <- list.dirs(path = getwd(), full.names = TRUE, recursive = FALSE)
  dataLocation    <- folderPaths[stringr::str_detect(basename(folderPaths), "Data")]
  resultsLocation <- folderPaths[stringr::str_detect(basename(folderPaths), "Results")]

  if (length(dataLocation) != 1) {
    stop("Expected exactly one 'Data' folder in the working directory; found ",
         length(dataLocation), ".", call. = FALSE)
  }
  if (length(resultsLocation) != 1) {
    stop("Expected exactly one 'Results' folder in the working directory; found ",
         length(resultsLocation), ".", call. = FALSE)
  }
  list(data = dataLocation, results = resultsLocation)
}
