# Regression: kinshipCalcs used a strict matrix lookup that errored with
# "subscript out of bounds" on missing (NA) alleles. It must instead treat
# NA / unknown codes as non-scorable and return without error.

make_pair <- function(allele1, allele2) {
  # one row per locus for a single coral pair
  data.frame(
    coral1  = "c1",
    coral2  = "c2",
    allele1 = allele1,
    allele2 = allele2,
    locus   = paste0("L", seq_along(allele1)),
    match   = allele1 == allele2,
    stringsAsFactors = FALSE
  )
}

test_that("kinshipCalcs runs when an allele is missing (NA)", {
  dat <- make_pair(c("A", "C", NA), c("A", "C", "G"))
  expect_no_error(res <- kinshipCalcs(dat))
  # only the two fully-scored loci count
  expect_equal(res$totalScorable, 2)
})

test_that("kinshipCalcs tolerates an unexpected code without erroring", {
  dat <- make_pair(c("A", "Z"), c("A", "C"))   # "Z" is not a valid IUPAC code
  expect_no_error(res <- kinshipCalcs(dat))
  expect_equal(res$totalScorable, 1)           # the "Z" locus drops out
})

test_that("kinshipCalcs computes the expected kinship for known genotypes", {
  # L1: A vs A  -> 1*1                = 1
  # L2: S(C/G) vs S(C/G) -> .5*.5 + .5*.5 = 0.5
  # avg_kinship = (1 + 0.5) / 2 = 0.75
  dat <- make_pair(c("A", "S"), c("A", "S"))
  res <- kinshipCalcs(dat)
  expect_equal(res$avg_kinship, 0.75)
})

test_that("identical homozygotes give kinship 1, opposite homozygotes give 0", {
  expect_equal(kinshipCalcs(make_pair("A", "A"))$avg_kinship, 1)
  expect_equal(kinshipCalcs(make_pair("A", "C"))$avg_kinship, 0)
})
