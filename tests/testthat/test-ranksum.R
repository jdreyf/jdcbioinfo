context("ranksum")

test_that("rank sum", {
  mat_z <- abs(multi_pval2z(mtt))
  rs <- ranksum(mat_z, nsim = 1e5-1)
  expect_equal(rownames(rs)[1:3], c("gene3", "gene41", "gene2"))
})


