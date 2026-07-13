#' Flag possible replicate (re-run) samples
#'
#' @description Detects `Coral_ID`s that look like repeat observations of the
#'   same physical colony, and emits an informative message so they can be
#'   reviewed before analysis. Replicate/re-run samples are typically named by
#'   appending a suffix to the original's ID (e.g., `ACER_001` and
#'   `ACER_001_Rerun`, or `ACER_001` and `ACER_001_1`). This function flags any
#'   `Coral_ID` that contains another `Coral_ID` as a punctuation-bounded prefix,
#'   which catches those cases regardless of the exact suffix wording (`_Rerun`,
#'   `_1`, `_rep`, etc.).
#'
#'   This is a **non-blocking** heads-up, not an error. Nothing is wrong with the
#'   computation: genet assignment handles any identifier correctly. The reason
#'   to review is that leaving replicates in the data inflates the genet count
#'   and biases diversity/kinship, because a replicate is a second observation of
#'   a colony that is already represented. Decide whether to remove replicates
#'   (or keep the higher-coverage record per colony) before interpreting results.
#' @param dataset A data frame with a `Coral_ID` column.
#' @returns Invisibly, a data frame of the flagged pairs (`replicate` and the
#'   `original` it appears to repeat), or `NULL` if none are found. Called mainly
#'   for its message side effect.
#' @export
flagPossibleReplicates <- function(dataset) {
  if (!("Coral_ID" %in% names(dataset))) return(invisible(NULL))

  ids <- unique(as.character(dataset$Coral_ID))
  if (length(ids) < 2) return(invisible(NULL))

  # For each id, does another id form a punctuation-bounded prefix of it?
  # e.g. base = "ACER_001" is a prefix of "ACER_001_Rerun", and the next
  # character after the base ("_") is punctuation -> a replicate-style name.
  hits <- lapply(ids, function(id) {
    others <- ids[ids != id]
    is_base <- vapply(others, function(o) {
      startsWith(id, o) &&
        nchar(id) > nchar(o) &&
        grepl("[[:punct:]]", substr(id, nchar(o) + 1, nchar(o) + 1))
    }, logical(1))
    if (any(is_base)) {
      # if several ids are prefixes, attribute to the longest (closest) one
      base <- others[is_base]
      base <- base[which.max(nchar(base))]
      data.frame(replicate = id, original = base, stringsAsFactors = FALSE)
    } else {
      NULL
    }
  })
  flagged <- do.call(rbind, hits)

  if (is.null(flagged) || nrow(flagged) == 0) return(invisible(NULL))

  message("Possible replicate/re-run samples detected: ", nrow(flagged),
          " colony ID(s) appear to repeat another colony (one ID contains ",
          "another as a suffix). Leaving replicates in the data inflates the ",
          "genet count and biases diversity/kinship. Review these and decide ",
          "whether to remove them before interpreting results.")
  message(paste0("  ", flagged$replicate, "  (repeats ", flagged$original, ")",
                 collapse = "\n"))

  invisible(flagged)
}
