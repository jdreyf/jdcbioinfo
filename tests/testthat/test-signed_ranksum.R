context("signed_ranksum")

test_that("signed rank sum", {
  mat_z <- multi_pval2z(mtt)
  rs <- signed_ranksum(mat_z, nsim = 1e5-1)
  expect_equal(rownames(rs)[1:3], c("gene3", "gene41", "gene2"))
})


