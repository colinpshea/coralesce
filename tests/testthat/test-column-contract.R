test_that("valid input passes the column contract", {
  ok <- data.frame(Coral_ID = c("a", "b"),
                   MatchMaker_Index = c("1", "2"),
                   L1 = c("A:A", "C:C"),
                   stringsAsFactors = FALSE)
  expect_silent(handleError_ColumnContract(ok))
  expect_identical(handleError_ColumnContract(ok), ok)  # returns input invisibly
})

test_that("wrong column 1 name is rejected", {
  bad <- data.frame(SampleID = "a", MatchMaker_Index = "1",
                    stringsAsFactors = FALSE)
  expect_error(handleError_ColumnContract(bad), "Coral_ID")
})

test_that("wrong column 2 name is rejected", {
  bad <- data.frame(Coral_ID = "a", Longitude = "1",
                    stringsAsFactors = FALSE)
  expect_error(handleError_ColumnContract(bad), "MatchMaker_Index")
})

test_that("non-integer MatchMaker_Index is rejected", {
  bad <- data.frame(Coral_ID = c("a", "b"),
                    MatchMaker_Index = c("1", "2.5"),
                    stringsAsFactors = FALSE)
  expect_error(handleError_ColumnContract(bad), "whole numbers")
})

test_that("missing MatchMaker_Index value is rejected", {
  bad <- data.frame(Coral_ID = c("a", "b"),
                    MatchMaker_Index = c("1", NA),
                    stringsAsFactors = FALSE)
  expect_error(handleError_ColumnContract(bad))
})

test_that("fewer than two columns is rejected", {
  bad <- data.frame(Coral_ID = "a", stringsAsFactors = FALSE)
  expect_error(handleError_ColumnContract(bad), "at least two columns")
})
