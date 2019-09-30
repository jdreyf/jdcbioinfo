context("pval2z")

test_that("Numeric input", {
  z <- pval2z(pval, direction=logfc)
  expect_equal(z[1:3], c(-0.48967977, -0.25665501, -0.48650188))
})

test_that("Character input", {
  z <- pval2z(pval, direction=direction)
  expect_equal(z[1:3], c(-0.48967977, 0.25665501, -0.48650188))
})

test_that("Invalid p-value input", {
  pval_invalid <- pval
  pval_invalid[1] <- 2
  expect_error(pval2z(pval_invalid, direction=direction))
})
