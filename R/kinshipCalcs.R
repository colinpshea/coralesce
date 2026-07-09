#' Calculate pairwise kinship across all comparisons
#'
#' @description Calculates pairwise kinship for each colony pair from their
#'   allele dosages at every fully scorable locus. This function does *not*
#'   exclude invariant loci; [kinshipCalcsNoInvar()] handles that.
#' @param dataset A pairwise allele-match table (columns `coral1`, `coral2`,
#'   `allele1`, `allele2`, `locus`, `match`) from
#'   [determineAllAlleleMatchesOthers()].
#' @returns A data frame with `coral1`, `coral2`, `totalProbLoci`,
#'   `totalScorable`, and `avg_kinship` for each colony pair.
#' @importFrom dplyr filter group_by summarise mutate n
#' @importFrom magrittr %>%
#' @export
kinshipCalcs <- function(dataset) {
  # Per-base allele dosage for each single-letter IUPAC code
  # (homozygous = 1; heterozygous = 0.5 for each constituent base).
  dosage <- rbind(
    A = c(A = 1,   C = 0,   G = 0,   T = 0),
    C = c(A = 0,   C = 1,   G = 0,   T = 0),
    G = c(A = 0,   C = 0,   G = 1,   T = 0),
    T = c(A = 0,   C = 0,   G = 0,   T = 1),
    M = c(A = 0.5, C = 0.5, G = 0,   T = 0),
    R = c(A = 0.5, C = 0,   G = 0.5, T = 0),
    W = c(A = 0.5, C = 0,   G = 0,   T = 0.5),
    S = c(A = 0,   C = 0.5, G = 0.5, T = 0),
    Y = c(A = 0,   C = 0.5, G = 0,   T = 0.5),
    K = c(A = 0,   C = 0,   G = 0.5, T = 0.5)
  )

  idx1 <- match(dataset$allele1, rownames(dosage))   # unknown code / NA -> NA
  idx2 <- match(dataset$allele2, rownames(dosage))
  p1 <- dosage[idx1, , drop = FALSE]                 # NA index -> row of NAs
  p2 <- dosage[idx2, , drop = FALSE]

  # A locus counts for a pair only if BOTH alleles are recognised IUPAC codes;
  # NA and any unrecognised code are treated as non-scorable. For clean data
  # (all present alleles valid) this equals the original NA-only rule.
  dataset$NumAlleles  <- 2 * (!is.na(idx1)) + 2 * (!is.na(idx2))
  dataset$multProbSum <- rowSums(p1 * p2)   # sum over A/C/G/T of the dosage products

  dataset %>%
    filter(NumAlleles == 4) %>%             # both colonies scored at this locus
    group_by(coral1, coral2) %>%
    summarise(totalProbLoci = sum(multProbSum),
              totalScorable = n(),
              .groups = "drop") %>%
    mutate(avg_kinship = totalProbLoci / totalScorable)
}
