#' Assign colonies to genets from clonal pairings
#'
#' @description Given the set of colony pairs judged to be clones (matching above
#'   threshold, with adequate data), assigns each colony to a genet. Colonies are
#'   nodes and clonal pairings are edges; each connected component is one genet.
#'   A colony that clones with no other colony forms its own single-member genet.
#'   The colony graph is built directly from the `coral1`/`coral2` columns, so
#'   colony identifiers may contain any characters (including periods).
#' @param clonePairs A data frame of clone pairings with (at least) `coral1` and
#'   `coral2` columns, as produced within [groupByGenets()]. Self-pairs
#'   (`coral1 == coral2`) are used to register singleton colonies.
#' @returns A data frame with columns `Coral_ID` and `genet` (an integer genet
#'   label), one row per colony.
#' @importFrom igraph graph_from_data_frame components
#' @export
returnGenetIdentity <- function(clonePairs) {
  colonies <- sort(unique(c(clonePairs$coral1, clonePairs$coral2)))
  if (length(colonies) == 0) {
    return(data.frame(Coral_ID = character(0), genet = integer(0),
                      stringsAsFactors = FALSE))
  }

  edges <- clonePairs[clonePairs$coral1 != clonePairs$coral2,
                      c("coral1", "coral2"), drop = FALSE]

  g <- graph_from_data_frame(
    d        = edges,
    directed = FALSE,
    vertices = data.frame(name = colonies, stringsAsFactors = FALSE)
  )
  memb <- components(g)$membership

  data.frame(
    Coral_ID = names(memb),
    genet    = as.integer(unname(memb)),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}
