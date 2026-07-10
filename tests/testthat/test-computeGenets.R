# computeGenets is the data-frame entry point for genet assignment. It should
# group clones the same way runGenets does, carry MatchMaker_Index through when
# present, and produce genet labels that collapseToGenets understands.

raw <- data.frame(
  Coral_ID = c("c1", "c2", "c3", "c4"),
  MatchMaker_Index = c(2, 1, 4, 3),
  L1 = c("A:A", "A:A", "C:C", "G:G"),
  L2 = c("C:G", "C:G", "A:T", "G:T"),
  L3 = c("T:T", "T:T", "T:T", "A:A"),
  L4 = c("A:G", "A:G", "A:G", "C:T"),
  stringsAsFactors = FALSE
)

test_that("clones share a genet; distinct colonies do not", {
  ga <- computeGenets(raw, PctMatchThreshold = 95, PctNotNullThreshold = 25)
  g  <- setNames(ga$genet, ga$Coral_ID)
  expect_equal(unname(g["c1"]), unname(g["c2"]))     # c1/c2 identical -> same genet
  expect_false(unname(g["c1"]) == unname(g["c3"]))
})

test_that("output carries MatchMaker_Index and is ordered by it", {
  ga <- computeGenets(raw, PctMatchThreshold = 95, PctNotNullThreshold = 25)
  expect_true("MatchMaker_Index" %in% names(ga))
  expect_false(is.unsorted(ga$MatchMaker_Index))
})

test_that("MatchMaker_Index is optional", {
  no_idx <- raw[, setdiff(names(raw), "MatchMaker_Index")]
  ga <- computeGenets(no_idx, PctMatchThreshold = 95, PctNotNullThreshold = 25)
  expect_false("MatchMaker_Index" %in% names(ga))
  expect_true(all(c("Coral_ID", "genet", "pctNull", "AdequateData") %in% names(ga)))
})

test_that("speciesCode prefixes the genet label", {
  ga <- computeGenets(raw, PctMatchThreshold = 95, PctNotNullThreshold = 25,
                      speciesCode = "ACER")
  expect_true(all(grepl("^ACER_", ga$genet)))
})

test_that("output feeds collapseToGenets to give one row per genet", {
  ga  <- computeGenets(raw, PctMatchThreshold = 95, PctNotNullThreshold = 25)
  red <- collapseToGenets(raw, ga)
  # c1/c2 are one genet, c3 and c4 distinct -> 3 representatives
  expect_equal(nrow(red), 3)
})

test_that("getPairwiseAlleleMatches returns a list with the pairwise table", {
  out <- computeGenets(raw, PctMatchThreshold = 95, PctNotNullThreshold = 25,
                       getPairwiseAlleleMatches = TRUE)
  expect_named(out, c("genetAssignment", "pairwiseAlleleMatches"))
  expect_s3_class(out$pairwiseAlleleMatches, "data.frame")
})

test_that("missing Coral_ID and missing thresholds error informatively", {
  expect_error(computeGenets(raw[, -1], PctMatchThreshold = 95, PctNotNullThreshold = 25),
               "Coral_ID")
  expect_error(computeGenets(raw), "Threshold")
})

test_that("single-letter coded input gives a helpful error", {
  coded <- data.frame(Coral_ID = c("c1", "c2"), L1 = c("A", "R"),
                      stringsAsFactors = FALSE)
  expect_error(computeGenets(coded, PctMatchThreshold = 95, PctNotNullThreshold = 25),
               "paired")
})
