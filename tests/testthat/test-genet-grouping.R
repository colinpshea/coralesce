# returnGenetIdentity assigns colonies to genets by connected components of the
# clone graph. The hardened version builds the graph directly from coral1/coral2
# so identifiers may contain ANY characters (periods, etc.).

test_that("clones are grouped and singletons get their own genet", {
  # c1==c2 are clones; c3 stands alone. Self-pairs register every colony.
  pairs <- data.frame(
    coral1 = c("c1", "c2", "c3", "c1"),
    coral2 = c("c1", "c2", "c3", "c2"),
    stringsAsFactors = FALSE
  )
  res <- returnGenetIdentity(pairs)
  expect_setequal(res$Coral_ID, c("c1", "c2", "c3"))
  # c1 and c2 share a genet; c3 differs
  g <- setNames(res$genet, res$Coral_ID)
  expect_equal(unname(g["c1"]), unname(g["c2"]))
  expect_false(unname(g["c1"]) == unname(g["c3"]))
})

test_that("colony IDs containing periods still group correctly", {
  # This is the exact case the old regex/separate logic broke on.
  pairs <- data.frame(
    coral1 = c("ACER.01", "ACER.02", "ACER.01"),
    coral2 = c("ACER.01", "ACER.02", "ACER.02"),
    stringsAsFactors = FALSE
  )
  res <- returnGenetIdentity(pairs)
  g <- setNames(res$genet, res$Coral_ID)
  expect_equal(unname(g["ACER.01"]), unname(g["ACER.02"]))
})

test_that("empty input yields an empty result, not an error", {
  empty <- data.frame(coral1 = character(0), coral2 = character(0),
                      stringsAsFactors = FALSE)
  expect_no_error(res <- returnGenetIdentity(empty))
  expect_equal(nrow(res), 0)
})
