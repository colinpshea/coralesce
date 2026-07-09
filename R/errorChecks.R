#' Validate the required input-column contract
#'
#' @description Enforces the input format the rest of MatchMakeR relies on:
#'   column 1 must be named `Coral_ID`, column 2 must be named
#'   `MatchMaker_Index`, and every `MatchMaker_Index` value must be a whole
#'   number (no missing values). If any condition fails, the function stops with
#'   an informative message rather than silently renaming or coercing columns.
#'   `MatchMaker_Index` is used only for bookkeeping (ordering output back to the
#'   database); it is not used in any genetic calculation.
#' @param dataset A data frame read from the input CSV, before any processing.
#' @returns Invisibly returns `dataset` if the contract is satisfied; otherwise
#'   throws an error.
#' @export
handleError_ColumnContract <- function(dataset) {
  nm <- names(dataset)

  if (length(nm) < 2) {
    stop("The input must have at least two columns: 'Coral_ID' (column 1) and ",
         "'MatchMaker_Index' (column 2).", call. = FALSE)
  }
  if (nm[1] != "Coral_ID") {
    stop("Column 1 must be named 'Coral_ID' but is named '", nm[1], "'. ",
         "Rename the first column to 'Coral_ID' and re-run.", call. = FALSE)
  }
  if (nm[2] != "MatchMaker_Index") {
    stop("Column 2 must be named 'MatchMaker_Index' but is named '", nm[2], "'. ",
         "This column must be an integer whose values match groupings present in ",
         "'Coral_ID'. Make 'MatchMaker_Index' the second column and re-run.",
         call. = FALSE)
  }

  idx <- dataset[[2]]
  num <- suppressWarnings(as.numeric(idx))
  ok  <- !is.na(num) & num == floor(num)
  if (!all(ok)) {
    bad <- unique(idx[!ok])
    stop("'MatchMaker_Index' must contain only whole numbers with no missing ",
         "values. Offending entries: ", paste(bad, collapse = ", "), ". ",
         "Clean 'MatchMaker_Index' and re-run.", call. = FALSE)
  }

  invisible(dataset)
}

#' Identify colonies with no valid SNP data
#'
#' @description Finds colonies that are `NA` at every locus. Such colonies
#'   cannot be assigned a genet from data and are handled separately: they are
#'   appended to the genet-assignment output with `genet = XXXX_NA`,
#'   `pctNull = 100`, and `AdequateData = No`. Emits a message listing the
#'   affected colonies when any are found.
#' @param dataset A data frame with `Coral_ID` in column 1 and single-letter
#'   allele codes in the remaining columns.
#' @returns A data frame of the all-`NA` colonies (with `genet` and
#'   `AdequateData` columns added), or `NULL` if there are none.
#' @export
handleError_allZeros <- function(dataset) {
  snp        <- dataset[, -1, drop = FALSE]
  allNArows  <- rowSums(is.na(snp)) == ncol(snp)
  allNA      <- dataset[allNArows, , drop = FALSE]

  if (nrow(allNA) == 0) return(NULL)

  allNA$genet        <- NA_integer_
  allNA$AdequateData <- "No"

  message("At least one colony has NA at all loci. These colonies were added to ",
          "the genet assignment with genet = XXXX_NA, pctNull = 100, and ",
          "AdequateData = No.")
  message("Offending colonies:\n", paste(allNA$Coral_ID, collapse = "\n"))
  allNA
}

#' Report columns that do not hold valid SNP data
#'
#' @description Checks every column other than the `Coral_ID` and
#'   `MatchMaker_Index` keys and reports (via message) any whose values are not
#'   all valid `IUPAC` allele pairs (e.g., a site-name column or an invalid base
#'   pair). Reporting only; the offending columns are actually dropped by
#'   [convertBasePairstoCodes()].
#' @param dataset A data frame with `Coral_ID`, `MatchMaker_Index`, SNP data, and
#'   possibly other fields.
#' @param acceptableData A reference table of allowable allele pairs (defaults to
#'   the package `IUPAC` table); must contain an `Allelepairs` column.
#' @returns Invisibly `NULL`; called for its message side effect.
#' @export
handleError_ProhibitedData <- function(dataset, acceptableData = IUPAC) {
  candidates <- setdiff(names(dataset), c("Coral_ID", "MatchMaker_Index"))
  if (length(candidates) == 0) return(invisible(NULL))

  allowable <- vapply(dataset[candidates],
                      function(col) all(col %in% acceptableData$Allelepairs),
                      logical(1))
  offenders <- candidates[!allowable]

  if (length(offenders) > 0) {
    message("Columns that do not match the required IUPAC base-pair format ",
            "(e.g., a site-name column or an invalid base pair) were removed ",
            "prior to genet and/or kinship calculations.")
    message("Offending columns:\n", paste(offenders, collapse = "\n"))
  }
  invisible(NULL)
}

#' Stop if any colony has multiple rows of SNP data
#'
#' @description Downstream reshaping requires one row per colony. This check
#'   stops with an informative error, listing the duplicated `Coral_ID`s, if any
#'   colony appears more than once.
#' @param df A data frame with a `Coral_ID` column.
#' @returns Invisibly `df` if all `Coral_ID`s are unique; otherwise an error.
#' @importFrom dplyr group_by summarise filter pull n
#' @importFrom magrittr %>%
#' @export
find_dups <- function(df) {
  dups <- df %>%
    group_by(Coral_ID) %>%
    summarise(N = n(), .groups = "drop") %>%
    filter(N > 1) %>%
    pull(Coral_ID)

  if (length(dups) > 0) {
    stop("At least one colony has multiple rows of SNP data, which is not ",
         "permitted. Contact the coral section to check whether these are ramets ",
         "or re-samples that can be renamed. Offending colonies:\n",
         paste(dups, collapse = "\n"), call. = FALSE)
  }
  invisible(df)
}
