# collapseToGenets keeps one representative per genet and matches rows on
# Coral_ID, so it must handle clones, singletons, and unassigned colonies.

geno <- data.frame(
  Coral_ID = c("a", "b", "c", "d", "e"),
  MatchMaker_Index = 1:5,
  L1 = c("A:A", "A:A", "C:C", "G:G", "A:A"),
  stringsAsFactors = FALSE
)

# a & b are one genet; c and d are singletons; e is unassigned (inadequate data)
ga <- data.frame(
  Coral_ID = c("a", "b", "c", "d", "e"),
  genet    = c("ACER_00001", "ACER_00001", "ACER_00002", "ACER_00003", "ACER_000NA"),
  pctNull  = c(10, 5, 20, 0, 100),
  AdequateData = c("Yes", "Yes", "Yes", "Yes", "No"),
  stringsAsFactors = FALSE
)

test_that("one representative is kept per real genet", {
  res <- collapseToGenets(geno, ga, drop_unassigned = TRUE)
  # 3 real genets -> 3 rows
  expect_equal(nrow(res), 3)
  # the a/b genet keeps b (lower pctNull = more data)
  expect_true("b" %in% res$Coral_ID)
  expect_false("a" %in% res$Coral_ID)
})

test_that("representative = 'first' keeps the first-listed colony", {
  res <- collapseToGenets(geno, ga, representative = "first", drop_unassigned = TRUE)
  expect_true("a" %in% res$Coral_ID)   # a is listed before b
  expect_false("b" %in% res$Coral_ID)
})

test_that("unassigned colonies are retained by default, dropped on request", {
  kept    <- collapseToGenets(geno, ga)                        # default keeps e
  dropped <- collapseToGenets(geno, ga, drop_unassigned = TRUE)
  expect_true("e" %in% kept$Coral_ID)
  expect_false("e" %in% dropped$Coral_ID)
})

test_that("unassigned colonies are not collapsed together", {
  ga2 <- ga
  ga2$Coral_ID <- c("a", "b", "c", "d", "e")
  ga2$genet <- c("G1", "G1", NA, "G_000NA", "G_000NA")  # c NA, d & e placeholder
  res <- collapseToGenets(geno, ga2)
  # one rep for G1 (a or b) + c + d + e all retained as individuals = 4
  expect_equal(nrow(res), 4)
})

test_that("column structure is preserved and only rows are filtered", {
  res <- collapseToGenets(geno, ga, drop_unassigned = TRUE)
  expect_identical(names(res), names(geno))
})

test_that("missing required columns raise informative errors", {
  expect_error(collapseToGenets(geno[, -1], ga), "Coral_ID")
  expect_error(collapseToGenets(geno, ga[, "Coral_ID", drop = FALSE]), "genet")
})
