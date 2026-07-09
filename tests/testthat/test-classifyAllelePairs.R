dat <- data.frame(
  Coral_ID = c("a", "b", "c"),
  L1 = c("A", "A", "C"),
  stringsAsFactors = FALSE
)

test_that("include_self = TRUE adds self-comparisons", {
  res <- classifyAllelePairs(dat, locus = 2, include_self = TRUE)
  # 3 cross pairs + 3 self pairs
  expect_equal(nrow(res), 6)
  expect_true(any(res$coral1 == res$coral2))
})

test_that("classifyAllelePairsOthers omits self-comparisons", {
  res <- classifyAllelePairsOthers(dat, locus = 2)
  expect_equal(nrow(res), 3)           # 3 cross pairs only
  expect_false(any(res$coral1 == res$coral2))
})

test_that("match flags equal alleles and is NA when an allele is missing", {
  d2 <- data.frame(Coral_ID = c("a", "b"), L1 = c("A", NA),
                   stringsAsFactors = FALSE)
  res <- classifyAllelePairsOthers(d2, locus = 2)
  expect_true(is.na(res$match))        # a vs NA -> NA, not FALSE
})

test_that("determineAllAlleleMatches errors on a single-column input", {
  expect_error(determineAllAlleleMatches(dat[, 1, drop = FALSE]),
               "at least one locus")
})
