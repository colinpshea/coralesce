test_that("calcPercentNotNull returns percent of non-NA values", {
  expect_equal(calcPercentNotNull(c(1, 2, NA, 4)), 75)
  expect_equal(calcPercentNotNull(c(NA, NA)), 0)
  expect_equal(calcPercentNotNull(1:10), 100)
})

test_that("omitInvariantLoci drops single-allele loci and keeps variable ones", {
  dat <- data.frame(
    Coral_ID = c("a", "b", "c"),
    Lvar = c("A", "C", "G"),   # variable -> keep
    Linv = c("T", "T", "T"),   # invariant -> drop
    stringsAsFactors = FALSE
  )
  res <- omitInvariantLoci(dat)
  expect_true("Lvar" %in% names(res))
  expect_false("Linv" %in% names(res))
  expect_true("Coral_ID" %in% names(res))
})

test_that("find_dups stops when a Coral_ID is duplicated", {
  dup <- data.frame(Coral_ID = c("a", "a", "b"), L1 = c("A", "C", "G"),
                    stringsAsFactors = FALSE)
  expect_error(find_dups(dup), "multiple rows")
})

test_that("find_dups passes silently when all Coral_IDs are unique", {
  uniq <- data.frame(Coral_ID = c("a", "b"), L1 = c("A", "C"),
                     stringsAsFactors = FALSE)
  expect_silent(find_dups(uniq))
})
