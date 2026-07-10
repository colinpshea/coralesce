# computeKinship is the data-frame entry point. It should accept raw paired-
# allele genotypes and match what the manual convert -> isolate -> omit ->
# kinship chain produces.

raw <- data.frame(
  Coral_ID = c("c1", "c2", "c3", "c4"),
  MatchMaker_Index = 1:4,
  L1 = c("A:A", "A:A", "C:C", "G:G"),
  L2 = c("C:G", "C:G", "A:T", "G:T"),
  L3 = c("A:G", "A:G", "A:G", "C:T"),
  stringsAsFactors = FALSE
)

test_that("computeKinship runs on a raw data frame and returns the expected shape", {
  k <- computeKinship(raw)
  expect_named(k, c("PopAvgMKGD", "MK_init", "MK_final"))
  expect_true(all(c("PopAvgMK", "PopAvgGD") %in% names(k$PopAvgMKGD)))
  expect_setequal(k$MK_init$Coral_ID, raw$Coral_ID)
  expect_null(k$MK_final)                      # subset = FALSE
})

test_that("computeKinship matches the manual pipeline", {
  manual_coded <- convertBasePairstoCodes(raw)
  manual_wd    <- isolateAllNAColonies(manual_coded)[[1]]
  manual <- kinshipCalcsNoInvar(manual_wd, subset = FALSE)
  auto   <- computeKinship(raw)
  expect_equal(auto$PopAvgMKGD, manual$PopAvgMKGD)
})

test_that("MatchMaker_Index is optional", {
  no_idx <- raw[, setdiff(names(raw), "MatchMaker_Index")]
  expect_no_error(computeKinship(no_idx))
})

test_that("missing Coral_ID errors informatively", {
  expect_error(computeKinship(raw[, -1]), "Coral_ID")
})

test_that("single-letter coded input gives a helpful error, not a silent empty result", {
  coded <- data.frame(Coral_ID = c("c1", "c2"), L1 = c("A", "R"),
                      stringsAsFactors = FALSE)
  expect_error(computeKinship(coded), "paired")
})
