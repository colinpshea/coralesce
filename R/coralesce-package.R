#' @keywords internal
"_PACKAGE"

# Quiet R CMD check NOTEs about "no visible binding for global variable ...".
# These names are column names referenced inside dplyr/tidyr non-standard
# evaluation, plus the lazy-loaded `IUPAC` data object. They are not true
# globals; declaring them here is the standard companion to importing dplyr/
# tidyr rather than attaching them via Depends.
#' @importFrom utils globalVariables
utils::globalVariables(c(
  "IUPAC", "Allelepairs", "IUPACallele",
  "Coral_ID", "MatchMaker_Index",
  "coral1", "coral2", "allele1", "allele2", "locus", "match",
  "pctMatch", "pctNotNull", "pctNull", "PartOfGenet", "AdequateData",
  "genet", "flag", "avg_kinship", "ind_mean_kinship", "N", "NumAlleles",
  "multProbSum", "totalProbLoci", "totalScorable"
))
