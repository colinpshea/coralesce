# The subset loop reuses the pairwise kinship table (filtering it) instead of
# recomputing it from genotypes each round. That shortcut is valid only because
# omitInvariantLoci() runs ONCE, before the loop, so the locus set is fixed.
#
# These tests pin the behaviour: the subset routine must still remove the
# most-related colony each round, retain exactly targetN colonies, and produce
# results consistent with a from-scratch recomputation on the retained set.

make_test_data <- function(n = 30, L = 25, seed = 123) {
  set.seed(seed)
  codes <- c("A", "R", "G", "C", "Y", "T")
  pm <- c(A = "A:A", R = "A:G", G = "G:G", C = "C:C", Y = "C:T", T = "T:T")
  mat <- replicate(L, sample(codes, n, replace = TRUE))
  raw <- data.frame(Coral_ID = paste0("c", seq_len(n)),
                    matrix(pm[mat], nrow = n),
                    stringsAsFactors = FALSE)
  names(raw)[-1] <- paste0("L", seq_len(L))
  convertBasePairstoCodes(raw)
}

test_that("subset retains exactly targetN colonies", {
  d <- make_test_data()
  k <- kinshipCalcsNoInvar(d, subset = TRUE, targetN = 20)
  expect_equal(nrow(k$MK_final), 20)
})

test_that("subset removes the most-related colonies, keeping the least-related", {
  d <- make_test_data()
  k <- kinshipCalcsNoInvar(d, subset = TRUE, targetN = 20)
  # the colony with the highest initial mean kinship should not survive a
  # subset that removes 10 colonies
  worst <- k$MK_init$Coral_ID[1]              # MK_init is sorted desc
  expect_false(worst %in% k$MK_final$Coral_ID)
})

test_that("retained colonies' kinship matches a from-scratch recomputation", {
  # This is the guard on the optimisation: filtering the pairwise table must
  # give the same numbers as recomputing kinship on just the retained colonies.
  d <- make_test_data()
  k <- kinshipCalcsNoInvar(d, subset = TRUE, targetN = 20)

  kept    <- k$MK_final$Coral_ID
  d_kept  <- d[d$Coral_ID %in% kept, , drop = FALSE]
  scratch <- kinshipCalcsNoInvar(d_kept, subset = FALSE)$MK_init

  a <- k$MK_final[order(k$MK_final$Coral_ID), ]
  b <- scratch[order(scratch$Coral_ID), ]
  expect_equal(a$ind_mean_kinship, b$ind_mean_kinship, tolerance = 1e-10)
})

test_that("subset = FALSE returns no MK_final and full MK_init", {
  d <- make_test_data()
  k <- kinshipCalcsNoInvar(d, subset = FALSE)
  expect_null(k$MK_final)
  expect_equal(nrow(k$MK_init), nrow(d))
})

test_that("population averages are unaffected by subsetting", {
  d <- make_test_data()
  a <- kinshipCalcsNoInvar(d, subset = FALSE)$PopAvgMKGD
  b <- kinshipCalcsNoInvar(d, subset = TRUE, targetN = 20)$PopAvgMKGD
  # PopAvgMKGD is computed on the FULL set in both cases
  expect_equal(a, b)
})
