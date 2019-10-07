context("combine_directional_pvalues")

test_that("only.p", {
  pval <- combine_directional_pvalues(mtt, only.p=TRUE)
  expect_equal(as.vector(pval[1:3]), c(0.000000e+00, 8.974161e-09, 1.341179e-07))
})

test_that("opposite direction", {
  res <- combine_directional_pvalues(mtt)
  expect_equal(as.vector(res[1:3, 3]), c(0.000000e+00, 4.487080e-07, 4.470598e-06))
})
