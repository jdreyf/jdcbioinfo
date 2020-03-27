context("revcomp")

test_that("single string", {
  x <- "ATCGATGC"
  y <- "GCATCGAT"
  z1 <- revcomp(x)
  z2 <- revcomp(y)
  expect_equal(y, z1)
  expect_equal(x, z2)
})

test_that("vector of strings", {
  x <- c("ATCGATGC", "ttaaggcc")
  y <- c("GCATCGAT", "ggccttaa")
  z1 <- revcomp(x)
  z2 <- revcomp(y)
  expect_equal(y, z1)
  expect_equal(x, z2)
})
