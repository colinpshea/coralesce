# These exercise the IUPAC-dependent path. They adapt to whatever the package's
# IUPAC table actually contains rather than hard-coding specific pairs.

test_that("convertBasePairstoCodes keeps Coral_ID, drops non-allele columns, and maps '?' to NA", {
  skip_if_not(exists("IUPAC"), "IUPAC data not available")
  # pick a real allele pair -> code from the package table (not the '?' row)
  realrow <- which(IUPAC$IUPACallele != "?")[1]
  pair <- IUPAC$Allelepairs[realrow]
  code <- IUPAC$IUPACallele[realrow]

  d <- data.frame(
    Coral_ID = c("a", "b"),
    L1 = c(pair, pair),
    L2 = c(pair, "?"),          # one unscored value
    Site = c("north", "south"), # non-allele column -> dropped
    stringsAsFactors = FALSE
  )
  suppressMessages(res <- convertBasePairstoCodes(d))

  expect_true("Coral_ID" %in% names(res))
  expect_false("Site" %in% names(res))
  expect_equal(as.character(res$L1), c(code, code))
  expect_true(is.na(res$L2[2]))          # "?" became NA
})

test_that("checkforAllowableData recognises valid pairs and rejects junk", {
  skip_if_not(exists("IUPAC"), "IUPAC data not available")
  valid <- IUPAC$Allelepairs[1]
  expect_true(checkforAllowableData(valid))
  expect_false(checkforAllowableData("not_an_allele"))
})
